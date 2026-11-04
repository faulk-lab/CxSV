#!/usr/bin/env python3
"""Sensitivity (literature events called, test A) against healthy Tier-1 calls per genome
(97-genome LOFO subset) for each long_span_rule, one point per rule
(docs/PAIRED_LINKAGE.md). Input: docs/item13b/rule_comparison.tsv.

  python calibration/oversized_clusters/fig_rule_tradeoff.py [TABLE] [OUT.png]
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
tab = sys.argv[1] if len(sys.argv) > 1 else os.path.join(REPO, "docs", "item13b", "rule_comparison.tsv")
out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(REPO, "docs", "figures", "item13", "fig_rule_tradeoff.png")
plt.rcParams.update({"font.size": 7, "axes.titlesize": 8, "axes.labelsize": 8, "xtick.labelsize": 6,
                     "ytick.labelsize": 6, "legend.fontsize": 6, "axes.spines.top": False, "axes.spines.right": False})
STYLE = {"v1.0 (reference)": ("v1.0 (panel-v1.0)", "#8d99ae", (0.0, -4.8, "center")),
         "exclude": ("exclude (main, default)", "#264653", (0.004, -3.4, "left")),
         "paired": ("paired (this branch)", "#1b6ca8", (0.006, 2.0, "left")),
         "breakends": ("breakends (de23a0e)", "#e76f51", (-0.006, 1.6, "right"))}
t = pd.read_csv(tab, sep="\t")
fig, ax = plt.subplots(figsize=(4.2, 3.4))
for _, r in t.iterrows():
    lab, col, (dx, dy, ha) = STYLE[r["rule"]]
    x = r["healthy_t1_per_genome"]
    ax.scatter(x, r["lit_called"], s=46, color=col, edgecolor="k", lw=0.4, zorder=3)
    ax.scatter(x, r["lit_tier1"], s=30, marker="D", facecolor="white", edgecolor=col, lw=1.2, zorder=3)
    ax.text(x + dx, r["lit_called"] + dy, lab, fontsize=6, ha=ha, color=col)
ax.set_xlabel("healthy panel genomes: Tier 1 calls per genome (97, LOFO)")
ax.set_ylabel("literature events called (of 63)")
ax.set_xlim(0.17, 0.40)
ax.set_ylim(0, 28)
ax.legend(handles=[Line2D([], [], marker="o", ls="", color="#555", label="events called"),
                   Line2D([], [], marker="D", ls="", mfc="white", mec="#555", label="events at Tier 1")],
          frameon=False, loc="center right")
ex, br = (t.set_index("rule").loc[k, "healthy_t1_per_genome"] for k in ("exclude", "breakends"))
ax.set_title(f"Only the breakends rule buys sensitivity, at +{100 * (br / ex - 1):.0f}% healthy Tier 1", loc="left")
fig.savefig(out, dpi=300, bbox_inches="tight")
print(out)
