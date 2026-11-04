#!/usr/bin/env python3
"""Report for the whole-genome positive control scored against a full panel release
(tests/lab/pos_control_full_genome.sh), with the 10-slice regression as the expectation.

Inputs (IN directory):
  score.events.tsv, score.junctions.tsv, score.offtarget.tsv   score_pos_control.py on the run
  cxsv_scores.tsv, cxsv_panel_status.tsv, cxsv_summary.tsv     the run's scoring/, private/, cxsv/
  locus_depth.tsv                                              the run's scoring/ (depth estimate)
  ref/regression10_item10.score.events.tsv                     10-slice regression (expectation)
  ref/slicecalls_vs_fullpanel_item10.score.events.tsv          slice calls re-scored vs the full panel
  standin/scoring__cxsv_scores.tsv, standin/scoring__locus_depth.tsv   optional second non-panel genome
  survey/oversized_clusters_per_genome.tsv, survey/oversized_clusters_truth_events.tsv
                                                               oversized-cluster survey of the panel
Healthy-genome reference rates: item-10 calibration, scorer as deployed, leave-one-family-out over
the 1,619 panel genomes (calibration/graded_privacy/results/scorer_exact_before_after.tsv and
strata_before_after.tsv, columns LOFO_deployed_T1_after / T2_per_genome_after).

Writes to OUT: pc_full_panel_events.tsv, pc_host_calls_summary.tsv, pc_host_tier12_calls.tsv,
pc_full_panel_figure.{png,pdf}.     Usage: full_panel_report.py IN OUT
"""
import re
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

TIERS = ["Tier1_high_priority", "Tier2", "Tier3", "Benign_likely"]
TIER_LABEL = {"Tier1_high_priority": "Tier 1", "Tier2": "Tier 2", "Tier3": "Tier 3", "Benign_likely": "Benign-likely"}
TIER_COLOR = {"Tier1_high_priority": "#b2182b", "Tier2": "#ef8a62", "Tier3": "#67a9cf", "Benign_likely": "#bdbdbd"}
# per healthy genome, deployed item-10 scorer, LOFO (None: not computed for that stratum)
HEALTHY = [("healthy panel genomes, all (n=1,619)", 0.4571, 4.3533),
           ("healthy, 20-30x (n=292)", 0.3596, None),
           ("healthy, >=30x (n=628)", 0.8678, None)]
KEY = ["Chr", "Start", "End"]


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


def depth_estimate(path):
    """Median mean depth (MAPQ>=20) over the flank intervals of locus_depth.tsv (bedcov -j)."""
    d = pd.read_csv(path, sep="\t", header=None)
    fl = d[d[4].str.startswith("flank")]
    return float((fl[6] / (fl[2] - fl[1])).median())


def breakpoint_offsets(junctions):
    """Largest |offset| (bp) between each truth junction and its closest matching call, per event."""
    def ends(s):
        m = re.match(r"(\w+):(\d+)-(\w+):(\d+)", s)
        return (m.group(1), int(m.group(2))), (m.group(3), int(m.group(4)))
    out = {}
    for _, r in junctions.iterrows():
        if r.status == "missed" or not isinstance(r.calls, str):
            continue
        a, b = ends(r.truth)
        best = np.inf
        for call in r.calls.split(";"):
            try:
                x, y = ends(call.split(":", 1)[1])
            except (AttributeError, IndexError):
                continue
            for p, q in ((x, y), (y, x)):
                if p[0] == a[0] and q[0] == b[0]:
                    best = min(best, max(abs(p[1] - a[1]), abs(q[1] - b[1])))
        if r.status == "exact" and np.isfinite(best):
            out.setdefault(r.event, []).append(best)
    return {e: int(max(v)) for e, v in out.items()}


def event_table(IN):
    ev = pd.read_csv(IN / "score.events.tsv", sep="\t", dtype=str)
    reg = pd.read_csv(IN / "ref/regression10_item10.score.events.tsv", sep="\t", dtype=str)
    resc = pd.read_csv(IN / "ref/slicecalls_vs_fullpanel_item10.score.events.tsv", sep="\t", dtype=str)
    scores = pd.read_csv(IN / "cxsv_scores.tsv", sep="\t")
    status = pd.read_csv(IN / "cxsv_panel_status.tsv", sep="\t")
    junc = pd.read_csv(IN / "score.junctions.tsv", sep="\t", dtype=str)
    summ = pd.read_csv(IN / "cxsv_summary.tsv", sep="\t")
    surv = pd.read_csv(IN / "survey/oversized_clusters_truth_events.tsv", sep="\t")
    sc = scores.merge(status[KEY + ["Panel_Junctions", "Panel_Junctions_Seen"]], on=KEY, how="left")
    sc["Rank"] = sc["Rank"].astype(str)
    t = ev[["event", "name", "gene", "junctions", "exact", "partial", "missed", "cxsv_loci", "classes", "ideal_class",
            "class_match", "expect_status", "panel_status", "carriers", "junctions_matched_bp", "score", "tier",
            "rank", "top_features", "verdict"]].copy()
    best = sc.set_index("Rank")
    t["matched_fraction"] = [
        f"{int(best.loc[r, 'Panel_Junctions_Seen'])}/{int(best.loc[r, 'Panel_Junctions'])}" if r in best.index else "."
        for r in t["rank"]]
    t["locus"] = [best.loc[r, "Locus_ID"] if r in best.index else "." for r in t["rank"]]
    offs = breakpoint_offsets(junc)
    t["max_breakpoint_offset_bp"] = [offs.get(e, ".") for e in t.event]
    # oversized clusters of THIS genome overlapping the event's truth span, and the panel-wide count
    over = summ[summ.CxSV_Criteria == "FILTERED_oversized_cluster"]
    notes = []
    for e in t.event:
        spans = surv[surv.event == e]
        hits = [f"{o.Chr}:{o.Start}-{o.End} ({o.Size / 1e6:.1f} Mb, {o.N_Records} records)"
                for _, s in spans.iterrows() for _, o in over.iterrows()
                if o.Chr == s.chrom and o.Start < s.end and o.End > s.start]
        notes.append("; ".join(sorted(set(hits))) or ".")
    t["oversized_cluster_here"] = notes
    pg = surv.groupby("event").panel_genomes_with_event_span_in_oversized_cluster.max()
    t["panel_genomes_span_in_oversized_cluster"] = [f"{int(pg.get(e, 0))}/1619" for e in t.event]
    for tag, d in (("reg10", reg), ("slicecalls_fullpanel", resc)):
        t = t.merge(d[["event", "panel_status", "carriers", "classes", "score", "tier", "verdict"]].rename(
            columns=lambda c: c if c == "event" else f"{tag}_{c}"), on="event", how="left")
    return t


def host_tables(IN):
    scores = pd.read_csv(IN / "cxsv_scores.tsv", sep="\t")
    off = pd.read_csv(IN / "score.offtarget.tsv", sep="\t")[KEY].drop_duplicates()
    host = scores.merge(off, on=KEY, how="inner")
    assert len(host) == len(off), "every off-target locus needs a score row"
    genomes = [("pos_control host 3559 (non-planted loci)", host, depth_estimate(IN / "locus_depth.tsv"))]
    sp = IN / "standin/scoring__cxsv_scores.tsv"
    if sp.exists():
        genomes.append(("HG01868 stand-in (all loci)", pd.read_csv(sp, sep="\t"),
                        depth_estimate(IN / "standin/scoring__locus_depth.tsv")))
    rows, top = [], []
    for name, d, dep in genomes:
        r = {"genome": name, "depth_x_median_flank_MAPQ20": round(dep, 1), "loci": len(d)}
        for s in ("private", "shared_junctions", "public"):
            r[f"status_{s}"] = int((d.Panel_Status == s).sum())
        for k in TIERS:
            r[k] = int((d.Tier == k).sum())
            r[f"{k}_private"] = int(((d.Tier == k) & (d.Panel_Status == "private")).sum())
        rows.append(r)
        x = d[d.Tier.isin(TIERS[:2])].sort_values("Total_Score", ascending=False).copy()
        x.insert(0, "genome", name)
        top.append(x)
    for name, t1, t2 in HEALTHY:
        rows.append({"genome": name + " (per genome, LOFO)", "Tier1_high_priority": t1,
                     "Tier2": t2 if t2 is not None else np.nan})
    keep = ["genome", "Rank", "Locus_ID", "Size", "Cmplx_Class", "SV_Types", "N_Junctions", "Panel_Status",
            "Panel_Carriers", "Panel_Callable", "Panel_Individuals", "Catalogues_Containing", "Genes_Affected",
            "Total_Score", "Tier", "Top_Features", "panel_privacy_value"]
    return host, pd.DataFrame(rows), pd.concat(top)[keep]


def figure(t, summary, IN, OUT):
    fig = plt.figure(figsize=(7.2, 2.9))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.6, 1.0, 1.0], wspace=0.75)
    # a: planted events
    ax = fig.add_subplot(gs[0])
    x = np.arange(len(t))
    for i, r in t.reset_index(drop=True).iterrows():
        s, s10 = _num(r.score), _num(r.reg10_score)
        if np.isfinite(s10):
            ax.scatter(i, s10, s=40, marker="_", color="k", lw=1.2, zorder=4)
        if np.isfinite(s):
            ax.plot([i, i], [0, s], color="0.6", lw=0.8, zorder=1)
            ax.scatter(i, s, s=30, color=TIER_COLOR[r.tier], edgecolor="k", lw=0.5, zorder=3)
        elif r.event == "E02":
            ax.scatter(i, 0, s=30, marker="x", color=TIER_COLOR["Tier1_high_priority"], lw=1.2, zorder=3)
            ax.annotate("E02 not called: inside a\n25.8 Mb oversized cluster", (i, 0), xytext=(i + 0.5, -0.33),
                        fontsize=5.5, ha="left", va="center", arrowprops=dict(arrowstyle="-", lw=0.5, color="0.3"))
    ev = list(t.event)
    i11, i12, i13 = ev.index("E11"), ev.index("E12"), ev.index("E13")
    ax.text((i11 + i12) / 2, 0.62, "simple SVs,\nnot CxSV\n(expected)", fontsize=5.5, ha="center",
            va="bottom", color="0.35")
    ax.text(i13 + 0.75, -0.8, "E13 private: carriers only in the\nregression panel's spliced genomes",
            fontsize=5.5, ha="right", va="center", color="0.2")
    for y, lab in ((0.99, "Tier 1"), (0.50, "Tier 2")):
        ax.axhline(y, ls="--", lw=0.6, color="0.45", zorder=0)
        ax.text(-0.6, y + 0.03, lab, fontsize=5.5, color="0.3", va="bottom")
    ax.axhline(0, lw=0.6, color="k", zorder=0)
    ax.set_xticks(x, t.event, rotation=90, fontsize=6)
    ax.set_ylim(-1.3, 3.0)
    ax.set_xlim(-0.8, len(t) - 0.2)
    ax.set_ylabel("CxSV score (points)")
    ax.set_title("11 of 13 planted events as expected", loc="left", fontsize=8, pad=10)
    hd = [mpl.lines.Line2D([], [], marker="o", ls="", color=TIER_COLOR[k], mec="k", mew=0.5, label=TIER_LABEL[k])
          for k in ("Tier1_high_priority", "Tier3")]
    hd.append(mpl.lines.Line2D([], [], marker="_", ls="", color="k", mew=1.2, ms=7, label="10-slice regression"))
    ax.legend(handles=hd, loc="upper left", fontsize=5.5, frameon=False, handletextpad=0.2, borderaxespad=0.2)

    # b: Tier-1 / Tier-2 calls per genome outside planted events vs healthy panel genomes
    ax2 = fig.add_subplot(gs[1])
    lab = {"pos_control host 3559 (non-planted loci)": "host 3559, 35x",
           "HG01868 stand-in (all loci)": "HG01868, 28x"}
    rows = list(summary.itertuples())
    short = {HEALTHY[0][0]: "healthy: all", HEALTHY[1][0]: "healthy: 20-30x", HEALTHY[2][0]: "healthy: >=30x"}
    names = [lab.get(r.genome, short.get(r.genome.replace(" (per genome, LOFO)", ""), r.genome)) for r in rows]
    y = np.arange(len(rows))[::-1]
    for yi, r in zip(y, rows):
        for k, dy, mk in (("Tier1_high_priority", 0.12, "o"), ("Tier2", -0.12, "s")):
            v = getattr(r, k)
            if pd.notna(v):
                ax2.scatter(v, yi + dy, s=22, marker=mk, color=TIER_COLOR[k], edgecolor="k", lw=0.5, zorder=3)
                ax2.text(v + 0.2, yi + dy, f"{v:g}" if v == int(v) else f"{v:.2f}", fontsize=5.5, va="center")
    ax2.axhline(len(rows) - 2.5, color="0.7", lw=0.5)
    ax2.set_yticks(y, names, fontsize=6)
    ax2.set_xlim(-0.3, 7.0)
    ax2.set_ylim(-0.6, len(rows) - 0.4)
    ax2.set_xlabel("calls per genome")
    ax2.set_title("No Tier-1 host call", loc="left", fontsize=8, pad=10)
    ax2.legend(handles=[mpl.lines.Line2D([], [], marker=m, ls="", color=TIER_COLOR[k], mec="k", mew=0.5,
                                         label=TIER_LABEL[k]) for k, m in (("Tier1_high_priority", "o"), ("Tier2", "s"))],
               loc="lower right", fontsize=5.5, frameon=False, handletextpad=0.2)


    # c: genome span inside oversized clusters (dropped whole by the classifier)
    ax3 = fig.add_subplot(gs[2])
    pg = pd.read_csv(IN / "survey/oversized_clusters_per_genome.tsv", sep="\t")
    panel = pg[pg.cohort != "extra"].oversized_span_bp / 1e6
    ax3.hist(panel, bins=np.arange(0, 1100, 40), color="0.75", edgecolor="0.4", lw=0.4)
    mark = {"pos_control_whole_genome": ("host 3559", "#b2182b"), "HG01868_standin": ("HG01868", "#2166ac"),
            "pos_control_10slice_regression": ("10-slice\nregression", "0.2")}
    ymax = ax3.get_ylim()[1] * 1.35
    ax3.set_ylim(0, ymax)
    side = {"pos_control_whole_genome": ("left", 12, 0.97), "HG01868_standin": ("right", -12, 0.97),
            "pos_control_10slice_regression": ("left", 12, 0.74)}
    for _, r in pg[pg.cohort == "extra"].iterrows():
        name, col = mark[r["sample"]]
        v = r.oversized_span_bp / 1e6
        ha, dx, yf = side[r["sample"]]
        ax3.axvline(v, color=col, lw=1.0)
        ax3.text(v + dx, ymax * yf, f"{name}\n{v:.0f} Mb", fontsize=5.5, color=col, va="top", ha=ha)
    ax3.set_xlabel("genome span in oversized clusters (Mb)")
    ax3.set_ylabel("panel genomes")
    ax3.set_title(f"Oversized clusters: median {panel.median():.0f} Mb", loc="left", fontsize=8, pad=10)
    ax3.set_xlim(-20, 1050)
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    for a_, letter in ((ax, "a"), (ax2, "b"), (ax3, "c")):    # letters left of each panel's labels
        bb = a_.get_tightbbox(rend).transformed(fig.transFigure.inverted())
        fig.text(max(bb.x0, a_.get_position().x0 - 0.06), bb.y1 + 0.01, letter, fontweight="bold",
                 fontsize=10, ha="left", va="bottom")
    fig.savefig(OUT / "pc_full_panel_figure.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / "pc_full_panel_figure.pdf", bbox_inches="tight")
    return fig


def main(IN, OUT):
    IN, OUT = Path(IN), Path(OUT)
    OUT.mkdir(parents=True, exist_ok=True)
    t = event_table(IN)
    t.to_csv(OUT / "pc_full_panel_events.tsv", sep="\t", index=False)
    host, summary, top = host_tables(IN)
    summary.to_csv(OUT / "pc_host_calls_summary.tsv", sep="\t", index=False)
    top.to_csv(OUT / "pc_host_tier12_calls.tsv", sep="\t", index=False)
    fig = figure(t, summary, IN, OUT)
    return t, summary, top, fig


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
