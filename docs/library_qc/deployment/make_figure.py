#!/usr/bin/env python3
"""Figure for docs/library_qc/DEPLOYMENT.md: Tier 1 before/after the conservative gate in the
deployment runs (a) and the gate on HG01879's Tier-1 candidates (b).
Usage: make_figure.py HG01879_runB_cxsv_scores.tsv OUT.png
Panel-a counts are the run results recorded in DEPLOYMENT.md (HG01868, positive control:
bundled panel-v1.1; HG04176, HG01879: run B, exclusion test release)."""
import re, sys
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

scores, out = sys.argv[1:3]
d = pd.read_csv(scores, sep="\t", low_memory=False)
gate = d["Library_QC_Gate"].astype(str)
c = d[gate.str.startswith(("kept", "demoted"))].copy()
cat = lambda s: ("demoted to Tier 2" if s.startswith("demoted") else
                 "kept: copy-number change" if "copy-number" in s else "kept: >half multi-caller junctions")
c["cat"] = c["Library_QC_Gate"].map(cat)
frac = lambda s: (lambda m: int(m[1]) / int(m[2]) if m else np.nan)(re.search(r"(\d+)/(\d+) junctions", s))
c["mf"] = c["Library_QC_Gate"].map(frac)
before = [0, 15, 0, int((c["cat"] != "").sum())]
after = [0, 15, 0, int((~c["cat"].str.startswith("demoted")).sum())]
labels = ["HG01868\n(PASS)", "Positive\ncontrol (PASS)", "HG04176\n(FLAG)", "HG01879\n(FLAG)"]

plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.1), gridspec_kw=dict(width_ratios=[1, 1.35], wspace=0.38))
ax = axes[0]; x = np.arange(4); w = 0.38
ax.bar(x - w / 2, before, w, color="#bdbdbd", label="Tier 1 before the gate")
ax.bar(x + w / 2, after, w, color="#3b6fb6", label="Tier 1 reported")
for xi, (b, a) in enumerate(zip(before, after)):
    for off, v in ((-w / 2, b), (w / 2, a)):
        ax.text(xi + off, v + 1.2, str(v), ha="center", va="bottom", fontsize=6)
ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=7); ax.set_ylabel("Tier-1 loci"); ax.set_ylim(0, 84)
ax.legend(frameon=False, loc="upper left"); ax.set_title("Gate acts only on flagged libraries", loc="left")
ax.text(-0.12, 1.06, "a", transform=ax.transAxes, fontweight="bold", fontsize=10)
ax = axes[1]
col = {"demoted to Tier 2": "#d98c2b", "kept: >half multi-caller junctions": "#3b6fb6", "kept: copy-number change": "#6a3d9a"}
mk = {"demoted to Tier 2": "o", "kept: >half multi-caller junctions": "s", "kept: copy-number change": "D"}
rng = np.random.default_rng(1)
for k in col:
    s = c[c["cat"] == k]
    ax.scatter(s["Total_Score"].astype(float), s["mf"].fillna(-0.15) + rng.uniform(-0.025, 0.025, len(s)), s=16,
               color=col[k], marker=mk[k], label=f"{k} (n={len(s)})", edgecolor="white", linewidth=0.4, zorder=3)
ax.axhline(0.5, color="grey", lw=0.7, ls="--", zorder=1)
ax.text(c["Total_Score"].astype(float).max(), 0.53, "gate: > 0.5", ha="right", va="bottom", fontsize=6, color="grey")
ax.set_xlabel("Total score (points unchanged by the gate)")
ax.set_ylabel("Junctions called by ≥ 2 callers (fraction)")
ax.set_yticks([-0.15, 0, 0.25, 0.5, 0.75, 1.0]); ax.set_yticklabels(["CN kept,\nnot computed", "0", "0.25", "0.5", "0.75", "1"])
ax.legend(frameon=False, loc="center left", bbox_to_anchor=(0.12, 0.72)); ax.set_ylim(-0.24, 1.08); ax.margins(x=0.05)
n_dem = int(c["cat"].str.startswith("demoted").sum())
ax.set_title(f"HG01879: {n_dem} of {len(c)} Tier-1 candidates demoted", loc="left")
ax.text(-0.12, 1.06, "b", transform=ax.transAxes, fontweight="bold", fontsize=10)
fig.savefig(out, dpi=300, bbox_inches="tight")
