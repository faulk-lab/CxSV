#!/usr/bin/env python3
"""Figures of the item-13 note (docs/OVERSIZED_CLUSTERS.md) from the collected tables.

  fig_blind_span.png        genome span inside oversized clusters per panel genome,
                            panel-v1.0 rule vs item 13 (all 1,619 genomes, re-clustered)
  fig_chaining_records.png  what makes clusters oversized: the records > 500 kb
                            (type, callers, length, recurrence, read depth)
  fig_pos_control.png       positive-control events before/after
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def groups(manifest, registry):
    man = pd.read_csv(manifest, sep="\t")
    reg = pd.read_csv(registry, sep="\t", low_memory=False)
    reg = reg[reg.status == "verified"]
    chem = {(r.cohort, r.individual_id): r.chemistry for r in reg.itertuples()}

    def g(r):
        if r.cohort == "1kGP_ONT_1019_CHM13":
            return "ONT 1kGP-ONT"
        if r.cohort == "1kGP_PacBio_CHM13":
            return "PacBio HiFi"
        return f"ONT LRSC {chem.get((r.cohort, r.individual_id), '')}"
    man["group"] = man.apply(g, axis=1)
    return man.set_index(["cohort", "sample_id"]).group


ORDER = ["ONT 1kGP-ONT", "ONT LRSC R9", "ONT LRSC R10", "PacBio HiFi"]
BEFORE, AFTER = "#9aa0a6", "#1f5aa6"


AFTER_VARIANT = "S500_exclude"      # the adopted rule (cxsv_params.long_span_rule exclude)


def fig_blind_span(gv, grp, out, style):
    gv = gv[gv.variant.isin(["V0_current", AFTER_VARIANT])].replace({AFTER_VARIANT: "S500_ends"}).copy()
    gv["group"] = [grp.get((c, s), "?") for c, s in zip(gv.cohort, gv["sample"])]
    gv["Mb"] = gv.blind_span / 1e6
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw={"width_ratios": [1.15, 1]})
    ax = axes[0]
    floor = 0.01
    for v, col, lab in (("V0_current", BEFORE, "panel-v1.0 rule"), ("S500_ends", AFTER, "item 13")):
        x = np.sort(np.maximum(gv[gv.variant == v].Mb.values, floor))
        ax.step(x, np.arange(1, len(x) + 1) / len(x), where="post", color=col, lw=1.6)
        med = np.median(gv[gv.variant == v].Mb)
        if v == "S500_ends":
            ax.text(0.015, 0.62, f"{lab}\nmedian {med:.3g} Mb", ha="left", va="center", color=col, fontsize=style[1])
        else:
            ax.text(med * 1.6, 0.35, f"{lab}\nmedian {med:.0f} Mb", ha="left", va="center", color=col,
                    fontsize=style[1])
    ax.set_xscale("log")
    ax.set_xlim(floor / 1.3, 3000)
    ax.set_xticks([0.01, 0.1, 1, 10, 100, 1000])
    ax.set_xticklabels(["0", "0.1", "1", "10", "100", "1,000"])
    ax.set_xlabel("genome span in oversized clusters (Mb, log)")
    ax.set_ylabel("fraction of panel genomes")
    ax.set_title(f"Blind span falls 200-fold (n = {int((gv.variant == 'V0_current').sum()):,} genomes)", loc="left")
    ax = axes[1]
    for k, g in enumerate(ORDER):
        sub = gv[gv.group == g]
        b = sub[sub.variant == "V0_current"].Mb.median()
        a = sub[sub.variant == "S500_ends"].Mb.median()
        n = len(sub[sub.variant == "V0_current"])
        ax.plot([a, b], [k, k], color="#cccccc", lw=1.2, zorder=1)
        ax.scatter([b], [k], color=BEFORE, s=28, zorder=2)
        ax.scatter([a], [k], color=AFTER, s=28, zorder=2)
        ax.text(b * 1.25, k, f"{b:.0f}", va="center", fontsize=style[2], color=BEFORE)
        ax.text(a / 1.25, k, f"{a:.2g}", va="center", ha="right", fontsize=style[2], color=AFTER)
        ax.text(1.02, k, f"n={n}", transform=ax.get_yaxis_transform(), va="center", fontsize=style[2])
    ax.set_yticks(range(len(ORDER)))
    ax.set_yticklabels(ORDER)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlim(0.12, 2500)
    ax.set_xticks([0.1, 1, 10, 100, 1000]); ax.set_xticklabels(["0.1", "1", "10", "100", "1,000"])
    ax.set_xlabel("median span per genome (Mb, log)")
    ax.set_title("Every platform: median before (grey), after (blue)", loc="left")
    ax.margins(y=0.12)
    for a_, l in zip(axes, "ab"):
        a_.text(-0.12, 1.06, l, transform=a_.transAxes, fontweight="bold", fontsize=style[0] + 2)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    return fig


def fig_chaining(gr, sites, dt, dr, out, style):
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.4))
    ax = axes[0, 0]
    gr = gr.copy()
    gr["origin"] = np.where(gr.bnd_alt == 1, "breakend pair\n(BND notation)", gr.svtype)
    sv = gr.supp_vec.astype(str).str.zfill(3)
    gr["callers"] = sv.map({"001": "Sniffles2 only", "010": "SVIM only", "011": "SVIM + Sniffles2"}).fillna("other")
    orig = ["breakend pair\n(BND notation)", "DEL", "DUP", "INV"]
    cal = ["Sniffles2 only", "SVIM only", "SVIM + Sniffles2"]
    cols = {"Sniffles2 only": "#7fa7d9", "SVIM only": "#c7d7ee", "SVIM + Sniffles2": "#1f5aa6"}
    ct = pd.crosstab(gr.origin, gr.callers).reindex(index=orig, columns=cal).fillna(0) / 1000
    left = np.zeros(len(orig))
    for c in cal:
        ax.barh(range(len(orig)), ct[c].values, left=left, color=cols[c], label=c, edgecolor="white", lw=0.5)
        left += ct[c].values
    for k, v in enumerate(left):
        ax.text(v + 0.4, k, f"{v:.1f}k", va="center", fontsize=style[2])
    ax.set_yticks(range(len(orig))); ax.set_yticklabels(orig); ax.invert_yaxis()
    ax.set_xlabel("records > 500 kb in oversized clusters (thousands)")
    ax.set_xlim(0, left.max() * 1.18)
    ax.legend(frameon=False, fontsize=style[1], loc="lower right")
    single = (gr.n_callers == 1).mean()
    ax.set_title(f"{single:.0%} single-caller; cuteSV calls none", loc="left")
    ax = axes[0, 1]
    bins = np.logspace(np.log10(5e5), np.log10(2.5e8), 40)
    ax.hist(gr.len, bins=bins, color="#7fa7d9", edgecolor="white", lw=0.3)
    ax.set_xscale("log")
    ax.set_xticks([1e6, 1e7, 1e8]); ax.set_xticklabels(["1 Mb", "10 Mb", "100 Mb"])
    ax.set_xlabel("record length (log)")
    ax.set_ylabel("records")
    ax.set_title(f"Median {np.median(gr.len) / 1e6:.1f} Mb, up to {gr.len.max() / 1e6:.0f} Mb", loc="left")
    ax = axes[1, 0]
    labs = ["1", "2-9", "10-99", ">=100"]
    rec = pd.cut(gr.site_genomes, [0, 1, 9, 99, 10 ** 6], labels=labs).value_counts(normalize=True).reindex(labs) * 100
    ax.bar(range(4), rec.values, color="#7fa7d9", width=0.7)
    for k, v in enumerate(rec.values):
        ax.text(k, v + 1.5, f"{v:.0f}%", ha="center", fontsize=style[2])
    ax.set_xticks(range(4)); ax.set_xticklabels(labs)
    ax.set_xlabel("panel genomes with the same record (ends within 100 kb)")
    ax.set_ylabel("% of records")
    ax.set_ylim(0, rec.max() * 1.15)
    top = sites.sort_values("genomes", ascending=False).iloc[0]
    ax.set_title(f"Mostly recurrent: top site in {int(top.genomes):,} genomes", loc="left")
    ax = axes[1, 1]
    d = dr.merge(dt, on="test_id")
    rng = np.random.default_rng(1)
    for k, (t, nc) in enumerate([("DEL", "1 caller"), ("DEL", "2 callers"), ("DUP", "1 caller"), ("DUP", "2 callers")]):
        s = d[(d.svtype == t) & (d.ncall == nc)].cn_ratio.values
        ax.scatter(k + rng.uniform(-0.18, 0.18, len(s)), np.minimum(s, 2.4), s=9, color="#7fa7d9", lw=0)
        ax.plot([k - 0.25, k + 0.25], [np.median(s)] * 2, color="#1f2937", lw=1.6)
        sup = np.mean(s <= 0.75) if t == "DEL" else np.mean(s >= 1.25)
        ax.text(k, 2.55, f"{sup:.0%}\nn={len(s)}", ha="center", fontsize=style[2])
    ax.axhline(1.0, color="#888888", lw=0.8, ls=":")
    ax.axhline(0.5, color="#bbbbbb", lw=0.8, ls="--"); ax.axhline(1.5, color="#bbbbbb", lw=0.8, ls="--")
    ax.text(3.4, 0.53, "het. loss", fontsize=style[2], va="bottom", color="#777777")
    ax.text(3.4, 1.53, "het. gain", fontsize=style[2], va="bottom", color="#777777")
    ax.set_xticks(range(4)); ax.set_xticklabels(["DEL\n1 caller", "DEL\n2 callers", "DUP\n1 caller", "DUP\n2 callers"])
    ax.set_ylim(0, 2.95); ax.set_xlim(-0.5, 4.2)
    ax.set_ylabel("read-depth ratio vs control genome")
    ax.set_title("Few change copy number (% supported)", loc="left")
    for a_, l in zip(axes.flat, "abcd"):
        a_.text(-0.12, 1.06, l, transform=a_.transAxes, fontweight="bold", fontsize=style[0] + 2)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    return fig


def fig_pos_control(ev_before, ev_after, host, out, style):
    """ev_*: score_pos_control.py event tables (panel-v1.0 rule, item 13); host: DataFrame
    rule, tier1, tier2 = Tier-1/Tier-2 calls of the host genome outside the planted events."""
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw={"width_ratios": [2.2, 1]})
    ax = axes[0]
    e = ev_before[["event", "score", "tier"]].merge(ev_after[["event", "score", "tier"]], on="event",
                                                     suffixes=("_b", "_a"))
    num = lambda v: float(v) if v not in (".", "", "NA") else np.nan
    for k, r in e.iterrows():
        b, a = num(r.score_b), num(r.score_a)
        if not np.isnan(b):
            ax.scatter(k - 0.15, b, s=26, color=BEFORE, zorder=3)
        else:
            ax.scatter(k - 0.15, 0, marker="x", s=30, color=BEFORE, zorder=3)
        if not np.isnan(a):
            ax.scatter(k + 0.15, a, s=26, color=AFTER, zorder=3)
        elif r.event in ("E11", "E12"):
            pass
        if r.event in ("E11", "E12"):
            ax.text(k, 0.08, "not\nCxSV", ha="center", fontsize=style[2], color="#555555")
    k2 = int(np.where(e.event == "E02")[0][0])
    ax.annotate("E02 (NF1): not called\nunder panel-v1.0 rule", (k2 - 0.15, 0), xytext=(k2 + 0.6, -0.55),
                fontsize=style[2], color="#555555", arrowprops=dict(arrowstyle="-", color="#999999", lw=0.6))
    ax.axhline(1.0, color="#888888", lw=0.7, ls="--"); ax.axhline(0.5, color="#bbbbbb", lw=0.7, ls="--")
    ax.text(len(e) - 0.4, 1.03, "Tier 1", fontsize=style[2], color="#777777", ha="right", va="bottom")
    ax.text(len(e) - 0.4, 0.53, "Tier 2", fontsize=style[2], color="#777777", ha="right", va="bottom")
    ax.set_xticks(range(len(e))); ax.set_xticklabels(e.event, rotation=90)
    ax.set_ylim(-0.8, 3.0); ax.set_xlim(-0.6, len(e) - 0.4)
    ax.set_ylabel("CxSV score (points)")
    ax.scatter([], [], color=BEFORE, s=26, label="panel-v1.0 rule"); ax.scatter([], [], color=AFTER, s=26, label="item 13")
    ax.legend(frameon=False, fontsize=style[1], loc="upper left")
    ax.set_title("E02 recovered; other planted events unchanged", loc="left")
    ax = axes[1]
    for k, r in host.reset_index(drop=True).iterrows():
        col = BEFORE if r.rule.startswith("panel") else ("#c0392b" if "breakend" in r.rule else AFTER)
        ax.barh(k, r.tier1, color=col, height=0.6)
        ax.text(r.tier1 + 0.08, k, f"{int(r.tier1)} Tier 1 / {int(r.tier2)} Tier 2", va="center", fontsize=style[2])
    ax.set_yticks(range(len(host))); ax.set_yticklabels(host.rule)
    ax.invert_yaxis(); ax.set_xlim(0, max(4, host.tier1.max() + 2.5))
    ax.set_xlabel("host Tier-1 calls (planted events excluded)")
    ax.set_title("Host genome", loc="left")
    for a_, l in zip(axes, "ab"):
        a_.text(-0.12, 1.06, l, transform=a_.transAxes, fontweight="bold", fontsize=style[0] + 2)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    return fig


def fig_tier_rates(per, out, style):
    """per: one row per genome x rule (rule, group, tier1, tier2), genomes scored in every rule."""
    rules = ["panel-v1.0 rule", "breakends (not adopted)", "item 13 (exclude)"]
    cols = {rules[0]: BEFORE, rules[1]: "#c0392b", rules[2]: AFTER}
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.9), sharey=True)
    groups_ = ["all"] + ORDER
    for ax, tier in zip(axes, ("tier1", "tier2")):
        for k, g in enumerate(groups_):
            sub = per if g == "all" else per[per.group == g]
            n = sub[sub.rule == rules[0]].shape[0]
            for j, r in enumerate(rules):
                v = sub[sub.rule == r][tier].mean()
                ax.scatter(v, k + (j - 1) * 0.22, color=cols[r], s=22, zorder=3, label=r if k == 0 else None)
            ax.text(1.0, k, f"n={n}", transform=ax.get_yaxis_transform(), va="center", fontsize=style[2])
        ax.set_yticks(range(len(groups_))); ax.set_yticklabels(["all genomes"] + ORDER)
        ax.invert_yaxis(); ax.margins(x=0.15)
        ax.set_xlabel(("Tier-1" if tier == "tier1" else "Tier-2") + " calls per healthy genome (LOFO)")
    axes[0].legend(frameon=False, fontsize=style[1], loc="lower right")
    axes[0].set_title("Excluding long records keeps Tier 1 close to v1.0", loc="left")
    for a_, l in zip(axes, "ab"):
        a_.text(-0.12 if l == "a" else -0.05, 1.06, l, transform=a_.transAxes, fontweight="bold", fontsize=style[0] + 2)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    return fig


def main():
    """figures.py GENOME_VARIANTS GIANT_RECORDS_ANNOTATED DOCS_TABLES SAMPLE_MANIFEST REGISTRY OUT_DIR
    GENOME_VARIANTS = collect.py genome_variants.tsv (characterise_v2); GIANT_RECORDS_ANNOTATED
    = giant records with site_genomes (tables in docs/item13 come from the same run)."""
    import sys
    gvp, grp_, T, man, reg, out = sys.argv[1:7]
    import os
    style = (8, 7, 6)
    try:
        from figure_style import apply_figure_style  # optional house style
        apply_figure_style(sizes=style)
    except Exception:
        plt.rcParams.update({"font.size": 8})
    grp = groups(man, reg)
    gv = pd.read_csv(gvp, sep="\t")
    g = pd.read_csv(grp_, sep="\t")
    sites = pd.read_csv(f"{T}/giant_record_sites.tsv", sep="\t")
    dep = pd.read_csv(f"{T}/giant_record_depth_tests.tsv", sep="\t")
    fig_blind_span(gv, grp, os.path.join(out, "fig_blind_span.png"), style)
    fig_chaining(g, sites, dep.drop(columns=["cn_ratio"]), dep[["test_id", "cn_ratio"]],
                 os.path.join(out, "fig_chaining_records.png"), style)
    pc = pd.read_csv(f"{T}/pos_control_events_before_after.tsv", sep="\t")
    before = pc.rename(columns={"score_v1.0": "score", "tier_v1.0": "tier"})[["event", "score", "tier"]]
    after = pc.rename(columns={"score_item13": "score", "tier_item13": "tier"})[["event", "score", "tier"]]
    fig_pos_control(before.astype(str), after.astype(str), pd.read_csv(f"{T}/pos_control_host_tiers.tsv", sep="\t"),
                    os.path.join(out, "fig_pos_control.png"), style)
    ts = pd.read_csv(f"{T}/lofo_tiers_per_genome.tsv", sep="\t")
    ts = ts[ts.scored == 1]
    lab = {"before": "panel-v1.0 rule", "breakends": "breakends (not adopted)", "exclude": "item 13 (exclude)"}
    n_arms = ts.groupby(["cohort", "sample"]).arm.nunique()
    keep = n_arms[n_arms == 3].index
    t = ts.set_index(["cohort", "sample"]).loc[keep].reset_index()
    out3 = t[(t.arm == "before") & (t.tier1 >= 30)].set_index(["cohort", "sample"]).index   # QC outliers
    t = t[~t.set_index(["cohort", "sample"]).index.isin(out3)]
    per = pd.DataFrame({"rule": t.arm.map(lab), "group": [grp.get((c, s)) for c, s in zip(t.cohort, t["sample"])],
                        "tier1": t.tier1, "tier2": t.tier2})
    fig_tier_rates(per, os.path.join(out, "fig_tier_rates.png"), style)


if __name__ == "__main__":
    main()
