#!/usr/bin/env python3
"""Per-event score breakdown of the positive control (roadmap item 5).

One horizontal bar per planted event: the feature contributions (points) of the
event's best-scoring CxSV locus, positive to the right of 0 and negative to the
left, the total as a black diamond, and the tier thresholds of
config/scoring.yaml as dashed lines. Events that are not CxSVs (E11, E12) are
listed as "not a CxSV (not scored)".

    plot_score_breakdown.py SCORE_PREFIX OUT_PREFIX [--scoring config/scoring.yaml]

SCORE_PREFIX is the score_pos_control.py prefix (<workdir>/score), so
<prefix>.events.tsv and <workdir>/score.features.tsv are read. Writes
OUT_PREFIX.png and OUT_PREFIX.pdf. Needs matplotlib and pyyaml.
"""
import argparse
import csv
import os

import yaml

GROUP_COLOURS = {          # outer level = hue family, inner = shades (feature order)
    "privacy": ["#1f5f8b", "#7fb2d6"],
    "genes": ["#b5412a", "#e08a6d", "#f2c2ae", "#7a2a1b"],
    "event": ["#6a6a6a", "#b0b0b0"],
    "evidence": ["#4d7d3a", "#94bf7f", "#c9e0bd"],
    "chromothripsis": ["#6b4c9a", "#9b7fc4", "#c5b3e0", "#e3d9f1", "#3f2a63"],
    "phenotype": ["#a07a1c", "#d9b65c"],
    "inheritance": ["#2f8f8f"],
}
LABELS = {
    "panel_privacy": "panel privacy (x callability)", "catalogue_frequency": "external catalogues",
    "gene_dosage": "gene dosage (ClinGen HI/TS, LOEUF)", "region_dosage": "ClinGen region dosage (above gene)", "exons": "exons disrupted", "gene_count": "gene count",
    "no_genes": "no gene affected", "complexity_class": "complexity class", "size": "size",
    "caller_agreement": "caller agreement", "read_support": "read support", "artefact_context": "repeat / segdup",
    "ct_interleaving": "1 interleaving", "ct_breakpoint_tier": "2 breakpoint tier",
    "ct_cn_oscillation": "5 copy-number change", "ct_clustering": "7 clustering",
    "ct_bfb_modifier": "8 fold-back (BFB) modifier", "hpo_match": "HPO match", "omim": "OMIM",
    "inheritance": "inheritance",
}


def read(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def plot(events_tsv, features_tsv, scoring_yaml, out_prefix):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D

    events = read(events_tsv)
    feats = read(features_tsv)
    tiers = yaml.safe_load(open(scoring_yaml))["tiers"]
    by_rank = {}
    for f in feats:
        by_rank.setdefault(f["Rank"], []).append(f)
    order = [f["feature"] for f in feats[: len({x["feature"] for x in feats})]]
    groups = {f["feature"]: f["group"] for f in feats}
    colour, used = {}, {}
    for feat in order:
        g = groups[feat]
        k = used.get(g, 0)
        colour[feat] = GROUP_COLOURS[g][k % len(GROUP_COLOURS[g])]
        used[g] = k + 1

    rows = [e for e in events]
    fig, ax = plt.subplots(figsize=(7.2, 0.36 * len(rows) + 1.6))
    nonzero = set()
    lo, hi = 0.0, 0.0
    for i, e in enumerate(rows):
        y = len(rows) - 1 - i
        if e["rank"] in (".", ""):
            ax.text(0.02, y, "not a CxSV (not scored)", va="center", ha="left", fontsize=6, color="#555555")
            continue
        pos, neg = 0.0, 0.0
        for f in by_rank[e["rank"]]:
            p = float(f["points"])
            if p == 0:
                continue
            nonzero.add(f["feature"])
            if p > 0:
                ax.barh(y, p, left=pos, color=colour[f["feature"]], height=0.66, edgecolor="white", linewidth=0.4)
                pos += p
            else:
                ax.barh(y, p, left=neg, color=colour[f["feature"]], height=0.66, edgecolor="white", linewidth=0.4)
                neg += p
        total = float(e["score"])
        lo, hi = min(lo, neg), max(hi, pos)
        ax.plot(total, y, marker="D", color="black", markersize=4, zorder=5)
        ax.text(max(pos, total) + 0.1, y, f"{total:.2f}  {e['tier'].replace('_high_priority', '')}",
                va="center", ha="left", fontsize=6)
    labels = [f"{e['event']}  {e['gene'] if e['gene'] != 'intergenic' else 'intergenic'}" for e in rows]
    ax.set_yticks(range(len(rows))[::-1])
    ax.set_yticklabels(labels)
    for lab, e in zip(ax.get_yticklabels(), rows):
        if e["gene"] != "intergenic":
            lab.set_fontstyle("italic")
    ax.axvline(0, color="black", linewidth=0.6)
    for t in tiers:
        if -10 < float(t["min"]) < 10 and float(t["min"]) > 0:
            ax.axvline(float(t["min"]), color="#888888", linestyle="--", linewidth=0.6)
            ax.text(float(t["min"]), len(rows) - 0.45, t["name"].split("_")[0].replace("Tier", "Tier "),
                    ha="center", va="bottom", fontsize=6, color="#555555")
    ax.set_xlim(lo - 0.2, hi + 1.0)
    ax.set_ylim(-0.7, len(rows) + 0.2)
    ax.set_xlabel("points (positive = towards pathogenic)")
    ax.set_title("Planted gene events score in Tier 1; intergenic and public events score below them",
                 loc="left", fontsize=8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    handles = [Patch(facecolor=colour[f], label=LABELS.get(f, f)) for f in order if f in nonzero]
    handles.append(Line2D([], [], marker="D", color="black", linestyle="none", markersize=4, label="total"))
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False, fontsize=6,
              title="feature (chromothripsis criteria numbered)", title_fontsize=6, labelspacing=0.6)
    fig.tight_layout()
    fig.savefig(out_prefix + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(out_prefix + ".pdf", bbox_inches="tight")
    return fig


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("score_prefix")
    ap.add_argument("out_prefix")
    ap.add_argument("--scoring", default=os.path.join(os.path.dirname(__file__), "..", "..", "config", "scoring.yaml"))
    a = ap.parse_args()
    import matplotlib
    matplotlib.use("Agg")
    matplotlib.rcParams.update({"font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7,
                                "xtick.labelsize": 6, "ytick.labelsize": 6, "font.family": "sans-serif"})
    plot(a.score_prefix + ".events.tsv", os.path.join(os.path.dirname(a.score_prefix), "score.features.tsv"),
         a.scoring, a.out_prefix)


if __name__ == "__main__":
    main()
