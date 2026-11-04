#!/usr/bin/env python3
"""Figures of roadmap item 11 (docs/LITERATURE_SET.md) from the tables in docs/item11.

  python docs/item11/make_figures.py --literature /mnt/sargasso1/CxSV/resources/literature/lit-v1.0 \\
      [--out docs/figures/item11]

  set_composition.png       events per source and author-described structure, junctions per
                            event, junction precision (literature_events/junctions.tsv)
  literature_scores.png     best locus score per junction-level event under four discovery
                            rules (literature_event_scores_modes.tsv, test A)
  detection_L01_L04.png     whole-genome positive control, L01-L04 + E02 (detection_*.tsv, test B)

Plain matplotlib + pandas; colours and layout as in the report.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.ticker import NullLocator  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({"font.size": 7, "axes.titlesize": 8, "axes.labelsize": 8, "xtick.labelsize": 6,
                     "ytick.labelsize": 6, "legend.fontsize": 6, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.linewidth": 0.8})
SRC_LAB = {"nazaryan2018": "Nazaryan-Petersen 2018", "collins2017": "Collins 2017", "eisfeldt2024": "Eisfeldt 2024",
           "bahrambeigi2019": "Bahrambeigi 2019", "grochowski2024": "Grochowski 2024",
           "clinvar_nstd102": "ClinVar (nstd102)"}
ORDER = ["nazaryan2018", "grochowski2024", "bahrambeigi2019", "eisfeldt2024", "collins2017", "clinvar_nstd102"]
MECHS = ["DUP-TRP/INV-DUP and higher-order gains", "DUP-NML-DUP", "DUP-INV and mixed CNV-inversion",
         "DEL-INV-DEL / deletions with inversion", "insertion / translocation",
         "chromothripsis / chromoanasynthesis", "inversion (ClinVar)", "complex substitution (ClinVar)"]
MCOL = dict(zip(MECHS, ["#1b6ca8", "#6aaed6", "#2a9d8f", "#8ab17d", "#e9c46a", "#d1495b", "#9e9e9e", "#cfcfcf"]))


def mechanism(t):
    """Group the authors' event-type wording into the eight structure classes of the figure."""
    t = t.lower()
    if "trp/inv" in t or "trip-quad" in t:
        return MECHS[0]
    if "chromothripsis" in t or "chromoanasynthesis" in t:
        return MECHS[5]
    if "insert" in t or "translocat" in t:
        return MECHS[4]
    if any(k in t for k in ("del-inv-del", "delinvdel", "deletions with inversion", "deletion with inversion")):
        return MECHS[3]
    if "inv" in t and "dup" in t:
        return MECHS[2]
    if "dup-nml-dup" in t or "tandem dup" in t:
        return MECHS[1]
    return MECHS[6] if t == "inversion" else MECHS[7]


def tiers(ax, horizontal=False):
    for x0, c in [(0.5, "#bdbdbd"), (0.99, "#7f7f7f")]:
        (ax.axhline if horizontal else ax.axvline)(x0, color=c, lw=0.8, ls="--", zorder=0)


def fig_composition(lit, out):
    E = pd.read_csv(f"{lit}/literature_events.tsv", sep="\t", comment="#")
    J = pd.read_csv(f"{lit}/literature_junctions.tsv", sep="\t", comment="#")
    E["mechanism"] = E.author_event_type.map(mechanism)
    ct = pd.crosstab(E.source_id, E.mechanism).reindex(index=ORDER, columns=MECHS, fill_value=0)
    fig = plt.figure(figsize=(7.2, 5.4))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.25, 1], hspace=0.45, wspace=0.12)
    ax = fig.add_subplot(gs[0, :])
    left = np.zeros(len(ORDER))
    for m in MECHS:
        ax.barh(range(len(ORDER)), ct[m].values, left=left, color=MCOL[m], label=m, edgecolor="white", lw=0.4)
        left += ct[m].values
    for i, t in enumerate(left):
        ax.text(t + 0.8, i, f"{int(t)}", va="center", fontsize=7)
    ax.set_yticks(range(len(ORDER)), [SRC_LAB[s] for s in ORDER])
    ax.invert_yaxis()
    ax.set_xlabel("events")
    ax.set_xlim(0, 100)
    ax.legend(frameon=False, fontsize=6, loc="upper right", ncol=2, title="author-described structure",
              title_fontsize=6)
    n_junction_level = int((E.source_id != "clinvar_nstd102").sum())
    ax.set_title(f"{len(E)} published pathogenic germline complex SVs; {n_junction_level} resolved to junctions",
                 loc="left")
    ax = fig.add_subplot(gs[1, 0])
    cur = E[E.source_id != "clinvar_nstd102"].copy()
    cur["nj"] = cur.event_id.map(J.groupby("event_id").size()).fillna(0)
    rng = np.random.default_rng(1)
    for i, s in enumerate(ORDER[:-1]):
        d = cur[cur.source_id == s]
        ax.scatter(d.nj, i + rng.uniform(-0.18, 0.18, len(d)), s=12, color=[MCOL[m] for m in d.mechanism],
                   edgecolor="k", lw=0.3, zorder=3)
        ax.plot([d.nj.median()] * 2, [i - 0.32, i + 0.32], color="k", lw=1.2)
    ax.set_xscale("log")
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xticks([1, 2, 5, 10, 20, 50], ["1", "2", "5", "10", "20", "50"])
    ax.set_yticks(range(len(ORDER) - 1), [SRC_LAB[s] for s in ORDER[:-1]])
    ax.invert_yaxis()
    ax.set_xlabel("junctions per event")
    ax.set_title(f"Curated events carry {int(cur.nj.min())}–{int(cur.nj.max())} junctions", loc="left")
    ax = fig.add_subplot(gs[1, 1])
    J["source_id"] = J.event_id.map(E.set_index("event_id").source_id)
    pc = pd.crosstab(J.source_id, J.precision).reindex(index=ORDER, columns=["bp", "approximate", "region"],
                                                       fill_value=0)
    lab = {"bp": "base pair", "approximate": "approximate (array / short read)", "region": "region (ClinVar)"}
    col = {"bp": "#264653", "approximate": "#8d99ae", "region": "#d9d9d9"}
    left = np.zeros(len(ORDER))
    for p in pc.columns:
        ax.barh(range(len(ORDER)), pc[p].values, left=left, color=col[p], label=lab[p], edgecolor="white", lw=0.4)
        left += pc[p].values
    ax.set_yticks(range(len(ORDER)), [])
    ax.invert_yaxis()
    ax.set_xlim(0, 125)
    ax.set_xlabel("junctions lifted to hs1")
    lifted = int((J.lift_status == "lifted").sum())
    ax.set_title(f"{lifted} of {len(J)} junctions lifted to hs1", loc="left")
    ax.legend(frameon=False, fontsize=6, loc="center right", bbox_to_anchor=(1.0, 0.55), title="junction precision",
              title_fontsize=6)
    fig.savefig(f"{out}/set_composition.png", dpi=300, bbox_inches="tight")


def score_of(x):
    return np.nan if x == "NC" or pd.isna(x) else float(str(x).split()[1])


def fig_scores(out):
    T = pd.read_csv(f"{HERE}/literature_event_scores_modes.tsv", sep="\t", index_col=0)
    srco = ["nazaryan2018", "eisfeldt2024", "bahrambeigi2019", "collins2017", "grochowski2024"]
    cols = {"v1.0_tier_score": "v1", "breakends_de23a0e_tier_score": "brk", "main_65ef970_tier_score": "main",
            "interval_rule_filter_off_tier_score": "inf"}
    for c, k in cols.items():
        T[k] = T[c].map(score_of)
    T["rank"] = T.source_id.map({s: i for i, s in enumerate(srco)})
    T["key"] = T[["main", "brk", "inf"]].max(axis=1).fillna(-1)
    T = T.sort_values(["rank", "key"], ascending=[True, False])
    n = len(T)
    y = np.arange(n)
    rc = {"called": "#1b6ca8", "long-span junctions excluded": "#e76f51",
          "split into non-complex clusters": "#9e9e9e", "oversized cluster (>500 kb)": "#6d597a"}
    fig, (ax0, ax) = plt.subplots(1, 2, figsize=(7.2, 9.6), sharey=True,
                                  gridspec_kw=dict(width_ratios=[0.06, 1], wspace=0.03))
    ax0.scatter(np.zeros(n), y, marker="s", s=18, c=[rc[x] for x in T.main_65ef970_outcome])
    labels = [f"{i.split('_', 1)[1]}  {g if len(str(g)) < 16 else str(g)[:14] + '…'}".replace("  .", "")
              for i, g in zip(T.index, T.genes_reported)]
    ax0.set_xticks([])
    ax0.set_yticks(y, labels, fontsize=6)
    for sp in ("top", "right", "bottom"):
        ax0.spines[sp].set_visible(False)
    tiers(ax)
    ax.scatter(T.inf, y, marker="x", s=12, color="#b0b0b0", lw=0.8, label="interval rule, size filter off")
    ax.scatter(T.brk, y, marker="^", s=16, facecolor="none", edgecolor="#8d99ae", lw=0.8,
               label="breakends rule (de23a0e, superseded)")
    ax.scatter(T.v1, y, marker="o", s=24, facecolor="none", edgecolor="#264653", lw=0.9,
               label="v1.0 (panel-v1.0, bundle 33c075a)")
    ax.scatter(T.main, y, marker="o", s=12, color="#1b6ca8", label="main 65ef970 (long spans excluded)")
    for s in srco:
        k = np.where(T.source_id.values == s)[0]
        if s != srco[0]:
            ax.axhline(k[0] - 0.5, color="#dddddd", lw=0.6)
        ax.text(-0.3, k.min() - 0.1, SRC_LAB[s], fontsize=6, va="top", color="#555", style="italic")
    ax.set_xlim(-0.35, 3.8)
    ax.set_ylim(n - 0.3, -1.6)
    ax.set_xlabel("best locus score (scoring v1, no phenotype)")
    ax.set_title(f"Main 65ef970 calls {int(T.main.notna().sum())} of {n} junction-level events; "
                 f"the breakends rule called {int(T.brk.notna().sum())}", loc="left", fontsize=8)
    l1 = ax.legend(frameon=False, fontsize=6, loc="lower right", bbox_to_anchor=(1.0, 0.12))
    ax.add_artist(l1)
    ax.legend(handles=[Line2D([], [], marker="s", ls="", color=v, label="main: " + k) for k, v in rc.items()],
              frameon=False, fontsize=6, loc="lower right")
    fig.savefig(f"{out}/literature_scores.png", dpi=300, bbox_inches="tight")


def fig_detection(out):
    D = pd.read_csv(f"{HERE}/detection_L01_L04.tsv", sep="\t")
    LL = pd.read_csv(f"{HERE}/detection_loci_with_L_breakends.tsv", sep="\t")
    rows = ["L01", "L02", "L03", "L04", "E02"]
    names = {"L01": "L01 PLP1 DUP-TRP/INV-DUP\n(Bahrambeigi BAB8964)", "L02": "L02 17p13.3 DEL-INV-DEL\n(Nazaryan P2109_123)",
             "L03": "L03 balanced chromothripsis\n(Collins UTR22)", "L04": "L04 Xq28 DEL + inv. ins.\n(Eisfeldt P6)",
             "E02": "E02 NF1 (planted;\nitem-13 reference)"}
    runs = {"bundle 33c075a": "v33c075a", "main 65ef970": "main"}
    fig = plt.figure(figsize=(7.2, 4.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1, 0.8], wspace=0.25)
    axa = fig.add_subplot(gs[0])
    axb = fig.add_subplot(gs[1], sharey=axa)
    axc = fig.add_subplot(gs[2])
    m = D[D.run == "main 65ef970"].set_index("event")
    for i, e in enumerate(rows):
        n, ex, lg = m.loc[e, "junctions"], m.loc[e, "exact"], m.loc[e, "junctions_over_500kb_excluded"]
        axa.barh(i, (ex - lg) / n, color="#264653", height=0.6)
        axa.barh(i, lg / n, left=(ex - lg) / n, color="#e76f51", height=0.6)
        axa.text(1.03, i, f"{ex}/{n}" + (f"\n{lg} >500 kb" if lg else ""), va="center", fontsize=6)
    axa.set_yticks(range(len(rows)), [names[e] for e in rows], fontsize=6)
    axa.invert_yaxis()
    axa.set_xlim(0, 1.45)
    axa.set_xticks([0, 0.5, 1], ["0", "50%", "100%"])
    axa.set_xlabel("junctions called exactly")
    axa.set_title("Every junction called", loc="left")
    axa.legend(handles=[Patch(color="#264653", label="used for discovery"),
                        Patch(color="#e76f51", label="excluded on main (>500 kb)")],
               frameon=False, fontsize=6, loc="upper center", bbox_to_anchor=(0.45, -0.13))
    tiers(axb)
    for i, e in enumerate(rows):
        for run, face, size, dy in (("bundle 33c075a", "none", 40, -0.13), ("main 65ef970", "#1b6ca8", 16, 0.13)):
            r = D[(D.run == run) & (D.event == e)].iloc[0]
            s = np.nan if r.best_score == "." else float(r.best_score)
            colr = "#264653" if face == "none" else "#1b6ca8"
            if np.isfinite(s):
                axb.scatter(s, i, s=size, facecolor=face, edgecolor=colr, lw=1.0, zorder=3)
            else:
                axb.scatter(-0.25, i + dy, marker="x", s=20, color=colr, zorder=3)
    axb.set_xlim(-0.45, 2.3)
    axb.tick_params(labelleft=False)
    axb.set_xticks([-0.25, 0, 0.5, 1, 1.5, 2], ["n.c.", "0", "0.5", "1", "1.5", "2"])
    axb.set_xlabel("best locus score (n.c. = not called)")
    axb.set_title("Tiers unchanged; E02 rescued", loc="left")
    axb.legend(handles=[Line2D([], [], marker="o", ls="", mfc="none", mec="#264653", label="bundle 33c075a (v1.0 discovery)"),
                        Line2D([], [], marker="o", ls="", color="#1b6ca8", ms=4, label="main 65ef970 (long spans excluded)")],
               frameon=False, fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.13))
    rng = np.random.default_rng(3)
    for j, (lab, k) in enumerate(runs.items()):
        v = LL[(LL.run == k) & (LL.event == "L03")].score.astype(float).values
        axc.scatter(j + rng.uniform(-0.15, 0.15, len(v)), v, s=14, edgecolor=["#264653", "#1b6ca8"][j],
                    facecolor=["none", "#1b6ca8"][j], lw=0.9, zorder=3)
        r = D[(D.run == lab) & (D.event == "L03")].iloc[0]
        axc.text(j, 2.25, f"{len(v)} loci\n{r.junctions_in_loci}/40\njunctions", ha="center", va="top", fontsize=6)
    tiers(axc, horizontal=True)
    axc.set_xticks([0, 1], ["bundle\n33c075a", "main\n65ef970"])
    axc.set_xlim(-0.6, 1.6)
    axc.set_ylim(-0.3, 2.3)
    axc.set_ylabel("score of loci with L03 breakends")
    axc.set_title("L03 found as fragments", loc="left")
    fig.savefig(f"{out}/detection_L01_L04.png", dpi=300, bbox_inches="tight")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--literature", required=True, help="literature set directory (lit-v1.0)")
    ap.add_argument("--out", default=os.path.join(HERE, "..", "figures", "item11"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    fig_composition(a.literature, a.out)
    fig_scores(a.out)
    fig_detection(a.out)


if __name__ == "__main__":
    main()
