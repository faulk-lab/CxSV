#!/usr/bin/env python3
"""Choose the graded privacy rule on the panel (roadmap item 10).

Input: loci table from collect_scores.py (every panel genome's own CxSV loci
scored by score_cxsv.py against the panel without its individual / family)
and the chance-match null of null_shift.py.

A candidate rule only changes the panel_privacy points of shared_junctions
loci:  points = shared + (public - shared) x r(f), where f is the matched
junction fraction and r a monotone weight in [0, 1]:
  ramp(f_lo, f_hi)   0 below f_lo, 1 above f_hi, linear between
  power(k)           f ** k
  complement(k)      1 - (1 - f) ** k        (the "(1 - f)^k" family)
Every other feature is unchanged (nonprivacy_points from the scorer), and the
catalogue feature is recomputed exactly as cxsv_scoring.f_catalogue_frequency
does with the new panel points (privacy = min(panel, catalogue) when the
locus is contained in a catalogue).

Two privacy views are evaluated:
  panel_only  total = nonprivacy + panel points  (no catalogue term). This is
              the clean leave-out view and the one used for choosing: the
              catalogues share individuals with the panel (HPRC, HGSVC3 and
              CoLoRSdb contain 1kGP genomes), so a panel genome's own loci are
              often "contained" in a catalogue because the same person is in it.
  deployed    total as the scorer computes it (panel + catalogue). Reported,
              not used for choosing (contaminated by that overlap).

Guard against penalising real events (chance matches). A truly private event
can share junctions with panel sites by coincidence; null_shift.py measures how
often, by moving every panel locus 1-5 Mb and matching it again. The proxy for
"true events" is the set of Tier-1 loci of healthy genomes that are PRIVATE
under LOFO (no junction in the panel): for each, the chance-match outcomes of
shifted loci with the same number of junctions (capped at 20) give the
probability that a rule would drop it below Tier 1 by chance (demotion risk;
the current rule already has one, because a chance match removes the private
bonus).

Selection rule:
  1. constraint: demotion risk at most --max-risk-increase (0.01 = 1 percentage
     point) above the current rule's;
  2. among feasible rules, the lowest mean Tier-1 calls per healthy genome
     (LOFO, panel_only, all genomes);
  3. parsimony: among feasible rules within --tolerance (2%) of that minimum,
     the one with the smallest demotion risk; ties -> ramp before power, then
     the larger f_lo;
  4. fraction: the scorer's matched fraction (Panel_Junctions_Seen /
     Panel_Junctions, f_seen) unless the best-individual fraction
     (Panel_Best_Junction_Fraction, f_best) gives > 2% fewer Tier-1 calls.
The first pre-registered constraint (at most 1% of shifted loci losing >= 0.10
points) was infeasible for every f_seen rule (2% of shifted loci have f_seen >
0.8) and allowed only ramps starting at f_best >= 0.7; it is still reported
(null_share_loss_ge) next to the risk. The unconstrained minimum is reported too.
Stability: the selection is repeated on random halves of the families (both
directions) and leaving out each superpopulation, and the chosen rule's
held-out Tier-1 counts are reported.
"""
import argparse
import itertools
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gp_common  # noqa: E402,F401  (puts workflow/lib on sys.path)
import cxsv_scoring as sc  # noqa: E402

T1, T2, T3 = 0.99, 0.50, 0.0
CATS = ("CoLoRSdb", "HGSVC3", "HPRC", "gnomAD-SV")


# ─── rules ───────────────────────────────────────────────────────────────────
def rule_grid():
    rules = [{"shape": "none"}]
    lo = np.round(np.arange(0.0, 0.751, 0.05), 2)
    hi = np.round(np.arange(0.05, 1.001, 0.05), 2)
    for a, b in itertools.product(lo, hi):
        if b > a:
            rules.append({"shape": "ramp", "f_lo": float(a), "f_hi": float(b)})
    for k in (0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0):
        rules.append({"shape": "power", "k": k})
    for k in (0.5, 1.0, 2.0, 3.0, 4.0):
        rules.append({"shape": "complement", "k": k})
    return rules


def rule_name(g):
    if g["shape"] == "none":
        return "current (shared_junctions 0)"
    if g["shape"] == "ramp":
        return f"ramp f_lo={g['f_lo']:.2f} f_hi={g['f_hi']:.2f}"
    if g["shape"] == "power":
        return f"power k={g['k']:g}"
    return f"complement k={g['k']:g}"


def weight(f, g):
    f = np.clip(f, 0, 1)
    if g["shape"] == "none":
        return np.zeros_like(f)
    if g["shape"] == "ramp":
        r = (f - g["f_lo"]) / (g["f_hi"] - g["f_lo"])
        r = np.where(f <= g["f_lo"] + 1e-12, 0.0, np.where(f >= g["f_hi"] - 1e-12, 1.0, r))
        return np.clip(r, 0, 1)
    if g["shape"] == "power":
        return f ** g["k"]
    return 1 - (1 - f) ** g["k"]


def check_against_scorer(g):
    """The vectorised weight equals cxsv_scoring.graded_weight (ramp / power)."""
    if g["shape"] not in ("ramp", "power"):
        return
    fs = np.linspace(0, 1, 41)
    ours = weight(fs, g)
    theirs = np.array([sc.graded_weight(x, g) for x in fs])
    assert np.allclose(ours, theirs), g


# ─── privacy re-scoring ──────────────────────────────────────────────────────
class Loci:
    def __init__(self, df, W):
        self.df = df
        pw = W["panel_privacy"]
        self.base, self.pub = float(pw["shared_junctions"]), float(pw["public"])
        self.shared = (df["Panel_Status"] == "shared_junctions").to_numpy()
        self.nonpriv = df["nonprivacy_points"].to_numpy(float)
        self.pp_now = df["panel_privacy_points"].to_numpy(float)
        self.cat_now = df["catalogue_frequency_points"].to_numpy(float)
        # catalogue: most severe containing catalogue (f_catalogue_frequency)
        cw = W["catalogue_frequency"]
        worst = np.zeros(len(df))
        anyc = np.zeros(len(df), dtype=bool)
        for c in CATS:
            cont = (df[f"{c}_Contains"] == "yes").to_numpy()
            af = pd.to_numeric(df[f"{c}_Max_AF"], errors="coerce").to_numpy(float)
            pen = np.where(af >= cw["common_af"], cw["common"],
                           np.where(af >= cw["rare_af"], cw["rare"], cw["seen"]))
            pen = np.where(np.isnan(af), cw["seen"], pen)
            worst = np.where(cont, np.minimum(worst, pen), worst)
            anyc |= cont
        self.cat_worst, self.cat_any = worst, anyc
        self.samples, self.sample_code = np.unique(df["sample_key"].to_numpy(), return_inverse=True)

    def catalogue_points(self, pp):
        combined = np.where(pp < 0, np.minimum(self.cat_worst, pp), self.cat_worst)
        return np.where(self.cat_any, np.round(combined - pp, 4), 0.0)

    def panel_points(self, f, g):
        pp = self.pp_now.copy()
        if g["shape"] != "none":
            r = weight(f, g)
            pp = np.where(self.shared, np.round(self.base + (self.pub - self.base) * r, 4), pp)
        return pp

    def totals(self, f, g, view):
        pp = self.panel_points(f, g)
        if view == "panel_only":
            return self.nonpriv + pp
        return self.nonpriv + pp + self.catalogue_points(pp)


def per_sample_counts(L, total, mask=None):
    m = np.ones(len(total), dtype=bool) if mask is None else mask
    n = len(L.samples)
    t1 = np.bincount(L.sample_code[m & (total >= T1 - 1e-9)], minlength=n)
    t2 = np.bincount(L.sample_code[m & (total >= T2 - 1e-9) & (total < T1 - 1e-9)], minlength=n)
    return t1, t2


def null_cost(null_f, null_shared, g, base, pub):
    """(share of null loci losing >= x points, mean loss) for shared null loci;
    null loci that are public by chance lose the public penalty under any rule
    and are counted separately (rule-independent)."""
    if g["shape"] == "none":
        loss = np.zeros(len(null_f))
    else:
        loss = np.where(null_shared, -(pub - base) * weight(null_f, g), 0.0)
    return loss


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--loci", required=True)
    ap.add_argument("--null", required=True)
    ap.add_argument("--scoring-yaml", required=True)
    ap.add_argument("--max-null-share", type=float, default=0.01, help="reported only")
    ap.add_argument("--max-risk-increase", type=float, default=0.01)
    ap.add_argument("--null-points", type=float, default=0.10)
    ap.add_argument("--tolerance", type=float, default=0.02)
    ap.add_argument("--low-coverage", type=float, default=5.0)
    ap.add_argument("--seed", type=int, default=20261103)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    import yaml
    W = yaml.safe_load(open(a.scoring_yaml))["weights"]

    df = pd.read_csv(a.loci, sep="\t", low_memory=False,
                     dtype={"Panel_Callable": str, "Catalogue_Max_AF": str})
    df["sample_key"] = df["cohort"] + "/" + df["sample_id"]
    df["f_seen"] = df["Panel_Matched_Fraction"].astype(float)
    df["f_best"] = df["Panel_Best_Junction_Fraction"].astype(float)
    df["platform_group"] = np.where(df["platform"] == "PACBIO", "PacBio HiFi",
                                    np.where(df["cohort"] == "1kGP_ONT_1019_CHM13", "ONT (1kGP-ONT 1019)",
                                             "ONT LRSC " + df["chemistry"].fillna("").astype(str)))
    meta_cols = ["sample_key", "cohort", "sample_id", "individual_id", "family_id", "n_panel_family",
                 "platform_group", "coverage_x", "superpopulation"]
    meta = df[df["mode"] == "LOFO"].drop_duplicates("sample_key")[meta_cols].set_index("sample_key").sort_index()
    null = pd.read_csv(a.null, sep="\t", low_memory=False)
    null = null[null["null_status"] != "unshiftable"]
    null_f = null["null_matched_fraction"].astype(float).to_numpy()
    null_best = null["null_best_fraction"].astype(float).to_numpy()
    null_shared = (null["null_status"] == "shared_junctions").to_numpy()
    null_public = (null["null_status"] == "public").to_numpy()

    null_n = np.minimum(null["n_junctions"].astype(int).to_numpy(), 20)
    null_pub_by_n = {b: int((null_public & (null_n == b)).sum()) for b in np.unique(null_n)}
    null_N_by_n = {b: int((null_n == b).sum()) for b in np.unique(null_n)}
    null_f_by_n = {fd: {b: arr[null_shared & (null_n == b)] for b in np.unique(null_n)}
                   for fd, arr in (("f_seen", null_f), ("f_best", null_best))}

    def demotion_risk(g, fdef, E_n, E_np, base, pub):
        """Mean over proxy true events of P(drop below Tier 1 | chance-match outcome of
        shifted loci with the same junction count) under rule g (panel_only view)."""
        out = np.zeros(len(E_n))
        for b in np.unique(E_n):
            idx = np.flatnonzero(E_n == b)
            fb_ = null_f_by_n[fdef].get(b, np.array([]))
            pen = np.sort(-(pub - base) * (weight(fb_, g) if g["shape"] != "none" else np.zeros(len(fb_))))
            N = null_N_by_n.get(b, 0)
            if not N:
                continue
            margin = E_np[idx] - T1          # private bonus lost on any chance match; demoted if penalty > margin
            n_shared_dem = len(pen) - np.searchsorted(pen, margin + 1e-9, side="right")
            n_pub_dem = np.where(E_np[idx] + pub < T1 - 1e-9, null_pub_by_n.get(b, 0), 0)
            out[idx] = (n_shared_dem + n_pub_dem) / N
        return out

    rules = rule_grid()
    for g in rules:
        check_against_scorer(g)
    out_curve, per_ind = [], {}
    Ls = {}
    for mode in ("LOFO", "LOIO"):
        d = df[df["mode"] == mode].reset_index(drop=True)
        L = Loci(d, W)
        assert np.array_equal(L.samples, meta.index.to_numpy()), "sample sets differ between modes"
        # the re-scoring reproduces the scorer exactly under the current rule
        tot_dep = L.totals(d["f_seen"].to_numpy(), {"shape": "none"}, "deployed")
        err = np.abs(np.round(tot_dep, 4) - d["Total_Score"].to_numpy(float)).max()
        cat_err = np.abs(L.catalogue_points(L.pp_now) - L.cat_now).max()
        print(f"{mode}: {len(d)} loci; re-scoring vs scorer max |total diff| {err:.2e}, catalogue {cat_err:.2e}",
              file=sys.stderr)
        assert err < 1e-6 and cat_err < 1e-6
        Ls[mode] = (L, d)
    lowcov = (meta["coverage_x"].astype(float) < a.low_coverage).to_numpy()
    strata = {"all": np.ones(len(meta), dtype=bool), f"coverage >= {a.low_coverage:g}x": ~lowcov,
              f"coverage < {a.low_coverage:g}x": lowcov}
    for pg in sorted(meta["platform_group"].unique()):
        strata[pg] = (meta["platform_group"] == pg).to_numpy()
    for sp in sorted(meta["superpopulation"].unique()):
        strata[f"superpopulation {sp}"] = (meta["superpopulation"] == sp).to_numpy()

    E_sets = {}
    for mode in ("LOFO", "LOIO"):
        L, d = Ls[mode]
        tot0 = L.totals(d["f_seen"].to_numpy(float), {"shape": "none"}, "panel_only")
        E = (d["Panel_Status"] == "private").to_numpy() & (tot0 >= T1 - 1e-9)
        E_sets[mode] = (E, np.minimum(d.loc[E, "Panel_Junctions"].astype(int).to_numpy(), 20), L.nonpriv[E],
                        L.sample_code[E])
        print(f"{mode}: {int(E.sum())} proxy true events (private Tier-1 loci of healthy genomes)", file=sys.stderr)
    for mode in ("LOFO", "LOIO"):
        L, d = Ls[mode]
        E, E_n, E_np, E_s = E_sets[mode]
        for fdef in ("f_seen", "f_best"):
            f = d[fdef].to_numpy(float)
            nf = null_f if fdef == "f_seen" else null_best
            for g in rules:
                loss = null_cost(nf, null_shared, g, L.base, L.pub)
                risk = demotion_risk(g, fdef, E_n, E_np, L.base, L.pub)
                row = {"mode": mode, "fraction": fdef, "rule": rule_name(g), **{k: g.get(k) for k in ("shape", "f_lo", "f_hi", "k")},
                       "null_share_loss_ge": float(np.mean(loss >= a.null_points - 1e-9)),
                       "null_mean_loss": float(loss.mean()),
                       "null_share_loss_ge_0.5": float(np.mean(loss >= 0.5 - 1e-9)),
                       "demotion_risk": float(risk.mean()), "_risk": risk}
                for view in ("panel_only", "deployed"):
                    tot = L.totals(f, g, view)
                    t1, t2 = per_sample_counts(L, tot)
                    for sname, sm in strata.items():
                        row[f"{view}|{sname}|T1_mean"] = float(t1[sm].mean()) if sm.any() else np.nan
                        row[f"{view}|{sname}|T2_mean"] = float(t2[sm].mean()) if sm.any() else np.nan
                    if view == "panel_only":
                        row["_t1"] = t1
                out_curve.append(row)
    curve = pd.DataFrame(out_curve)

    # ─── selection (LOFO, f_seen and f_best compared; panel_only, all genomes) ──
    def select(c, mask_samples=None, risk_mask=None):
        c = c.copy()
        if mask_samples is not None:
            c["obj"] = [t[mask_samples].mean() for t in c["_t1"]]
        else:
            c["obj"] = c["panel_only|all|T1_mean"]
        if risk_mask is not None:
            c["risk"] = [r[risk_mask].mean() for r in c["_risk"]]
        else:
            c["risk"] = c["demotion_risk"]
        r0 = c[c["shape"] == "none"]["risk"].iloc[0]
        feas = c[(c["risk"] <= r0 + a.max_risk_increase + 1e-12) & (c["shape"] != "none")]
        best = feas["obj"].min()
        near = feas[feas["obj"] <= best * (1 + a.tolerance) + 1e-12].copy()
        near["shape_rank"] = near["shape"].map({"ramp": 0, "power": 1, "complement": 2})
        near = near.sort_values(["risk", "shape_rank", "f_lo"], ascending=[True, True, False])
        return near.iloc[0], best, feas

    lofo = curve[curve["mode"] == "LOFO"]
    choice = {}
    for fdef in ("f_seen", "f_best"):
        ch, best, feas = select(lofo[lofo["fraction"] == fdef])
        choice[fdef] = (ch, best, len(feas))
        print(f"{fdef}: chosen {ch['rule']}: LOFO T1/genome {ch['panel_only|all|T1_mean']:.3f} "
              f"(feasible min {best:.3f}, {len(feas)} feasible); demotion risk {ch['demotion_risk']:.4f}; "
              f"null share >= {a.null_points}: {ch['null_share_loss_ge']:.4f}", file=sys.stderr)
    # f definition: the scorer's matched fraction (f_seen) unless f_best is clearly better (> tolerance)
    fsel = "f_seen"
    if choice["f_best"][0]["panel_only|all|T1_mean"] < choice["f_seen"][0]["panel_only|all|T1_mean"] * (1 - a.tolerance):
        fsel = "f_best"
    chosen = choice[fsel][0]
    g = {k: chosen[k] for k in ("shape", "f_lo", "f_hi", "k") if pd.notna(chosen[k])}

    # ─── stability ───────────────────────────────────────────────────────────
    rng = np.random.default_rng(a.seed)
    fams = meta["family_id"].unique()
    rng.shuffle(fams)
    half = set(fams[: len(fams) // 2])
    inA = meta["family_id"].isin(half).to_numpy()
    lf = lofo[lofo["fraction"] == fsel]
    stab = []
    splits = [("families half A -> B", inA, ~inA), ("families half B -> A", ~inA, inA)]
    for sp in sorted(meta["superpopulation"].unique()):
        m = (meta["superpopulation"] == sp).to_numpy()
        splits.append((f"leave out {sp}", ~m, m))
    none_t1 = lf[lf["shape"] == "none"]["_t1"].iloc[0]
    E_s = E_sets["LOFO"][3]
    for name, tr, te in splits:
        ch, _, _ = select(lf, tr, risk_mask=tr[E_s])
        t1 = ch["_t1"]
        chosen_t1 = lf[lf["rule"] == chosen["rule"]]["_t1"].iloc[0]
        stab.append({"split": name, "n_train": int(tr.sum()), "n_test": int(te.sum()), "rule_chosen_on_train": ch["rule"],
                     "test_T1_before": float(none_t1[te].mean()), "test_T1_rule_from_train": float(t1[te].mean()),
                     "test_T1_global_choice": float(chosen_t1[te].mean())})
    stab = pd.DataFrame(stab)

    # ─── per-individual before / after ───────────────────────────────────────
    rows = meta.copy()
    resid = []
    for mode in ("LOFO", "LOIO"):
        L, d = Ls[mode]
        f = d[fsel].to_numpy(float)
        for view in ("panel_only", "deployed"):
            for tag, rule in (("before", {"shape": "none"}), ("after", g)):
                tot = L.totals(f, rule, view)
                t1, t2 = per_sample_counts(L, tot)
                rows[f"{mode}_{view}_T1_{tag}"] = t1
                rows[f"{mode}_{view}_T2_{tag}"] = t2
                if mode == "LOFO":
                    sel = tot >= T1 - 1e-9
                    x = d.loc[sel, ["sample_key", "Locus_ID", "Panel_Status", "Cmplx_Class", "Panel_Junctions",
                                     "Panel_Junctions_Seen", "f_seen", "f_best", "Top_Features"]].copy()
                    x["view"], x["rule"], x["total"] = view, tag, tot[sel]
                    resid.append(x)
        rows[f"{mode}_n_loci"] = np.bincount(L.sample_code, minlength=len(L.samples))
        rows[f"{mode}_n_shared"] = np.bincount(L.sample_code[L.shared], minlength=len(L.samples))
        rows[f"{mode}_n_private"] = np.bincount(L.sample_code[(d["Panel_Status"] == "private").to_numpy()],
                                                minlength=len(L.samples))
    rows["low_coverage_flag"] = np.where(lowcov, "yes", "no")
    rows.reset_index().to_csv(f"{a.outdir}/per_individual_leaveout.tsv", sep="\t", index=False)
    resid = pd.concat(resid)
    resid.to_csv(f"{a.outdir}/tier1_loci_LOFO.tsv.gz", sep="\t", index=False)
    curve.drop(columns=["_t1", "_risk"]).to_csv(f"{a.outdir}/sweep_curve.tsv", sep="\t", index=False)
    stab.to_csv(f"{a.outdir}/stability.tsv", sep="\t", index=False)

    # matched fraction distributions (LOFO) by class and status
    L, d = Ls["LOFO"]
    dist = d[d["Panel_Status"] == "shared_junctions"].groupby("Cmplx_Class")[["f_seen", "f_best"]].describe()
    dist.to_csv(f"{a.outdir}/matched_fraction_by_class_LOFO.tsv", sep="\t")
    status_counts = d.groupby("Panel_Status").size()
    json.dump({
        "fraction": fsel, "graded": {"enabled": True, **{k: (float(v) if k != "shape" else v) for k, v in g.items()}},
        "rule": chosen["rule"],
        "selection": {"max_risk_increase": a.max_risk_increase, "max_null_share_reported": a.max_null_share,
                      "null_points": a.null_points, "tolerance": a.tolerance,
                      "feasible_min_T1_per_genome": float(choice[fsel][1]), "n_feasible": int(choice[fsel][2])},
        "LOFO_panel_only_T1_per_genome": {"before": float(none_t1.mean()), "after": float(chosen["panel_only|all|T1_mean"])},
        "alternative_fraction": {k: {"rule": v[0]["rule"], "T1_per_genome": float(v[0]["panel_only|all|T1_mean"]),
                                     "demotion_risk": float(v[0]["demotion_risk"])} for k, v in choice.items()},
        "demotion_risk": {"current": float(lf[lf["shape"] == "none"]["demotion_risk"].iloc[0]),
                          "chosen": float(chosen["demotion_risk"]),
                          "n_proxy_true_events": int(E_sets["LOFO"][0].sum())},
        "unconstrained_minimum": {k: {"rule": r["rule"], "T1_per_genome": float(r["panel_only|all|T1_mean"]),
                                      "demotion_risk": float(r["demotion_risk"])}
                                  for k, r in ((fd, lofo[(lofo["fraction"] == fd) & (lofo["shape"] != "none")]
                                                .sort_values(["panel_only|all|T1_mean", "demotion_risk"]).iloc[0])
                                               for fd in ("f_seen", "f_best"))},
        "null": {"n_loci": int(len(null)), "share_shared": float(null_shared.mean()),
                 "share_public_by_chance": float(null_public.mean()),
                 "chosen_share_loss_ge": float(chosen["null_share_loss_ge"]),
                 "chosen_mean_loss": float(chosen["null_mean_loss"])},
        "LOFO_status_counts": {k: int(v) for k, v in status_counts.items()},
        "n_genomes": int(len(meta)), "n_low_coverage": int(lowcov.sum()),
    }, open(f"{a.outdir}/chosen_parameters.json", "w"), indent=2)
    print(open(f"{a.outdir}/chosen_parameters.json").read(), file=sys.stderr)


if __name__ == "__main__":
    main()
