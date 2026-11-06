#!/usr/bin/env python3
"""Figure for docs/library_qc/DEPLOYMENT.md (rebuild 2026-11-06): z of the two flag metrics under the
auto-detected group (qc_groups.tsv) and under the previous bundle's ONT_1kGP_1019 default
(../deployment/deployment_library_qc_rows.tsv). Usage: make_figure_groups.py QC_GROUPS PREVIOUS OUT.png"""
import sys
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

qg, prev, out = sys.argv[1:4]
q = pd.read_csv(qg, sep="\t").set_index("sample"); pv = pd.read_csv(prev, sep="\t").set_index("sample")
order = ["HG01868", "pos_control", "PC_host_3559", "HG04176", "HG01879_PacBio"]
lab = {"HG01868": "HG01868 (stand-in)", "pos_control": "Positive control", "PC_host_3559": "PC host 3559",
       "HG04176": "HG04176", "HG01879_PacBio": "HG01879 PacBio"}
col = {"ONT_LRSC_R9": "#3b6fb6", "ONT_LRSC_R10": "#d98c2b", "PacBio_HiFi": "#6a3d9a", "ONT_1kGP_1019": "#7f7f7f"}
plt.rcParams.update({"font.size": 7, "axes.spines.top": False, "axes.spines.right": False})
fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.9), sharey=True, gridspec_kw=dict(wspace=0.08))
y = np.arange(len(order))[::-1]
for ax, m, t in zip(axes, ["z_svim_only_bnd_per_gb", "z_svim_only_any_per_gb"],
                    ["SVIM-only BND/TRA per Gb", "SVIM-only records per Gb"]):
    ax.axvline(6, color="#b2182b", lw=0.8, ls="--", zorder=1); ax.axvline(0, color="0.85", lw=0.6, zorder=0)
    for yi, s in zip(y, order):
        r = q.loc[s]
        ax.scatter(r[m], yi, s=34, color=col[r.qc_group], zorder=3, edgecolor="white", linewidth=0.5)
        if s in pv.index:
            ax.scatter(pv.loc[s, m], yi, s=30, facecolor="white", edgecolor=col["ONT_1kGP_1019"], linewidth=1, zorder=2)
            ax.plot([pv.loc[s, m], r[m]], [yi, yi], color="0.75", lw=0.7, zorder=1)
    ax.set_xlabel(f"robust z, {t}"); ax.set_xlim(-2.5, 20.5)
    ax.text(6.2, y[0] + 0.45, "flag z > 6", color="#b2182b", fontsize=6, va="center")
axes[0].set_yticks(y); axes[0].set_yticklabels([lab[s] for s in order]); axes[0].set_ylim(-0.6, len(order) - 0.3)
h = [Line2D([], [], marker="o", ls="", color=col[g], label=f"detected group: {g}") for g in
     ("ONT_LRSC_R9", "ONT_LRSC_R10", "PacBio_HiFi")]
h.append(Line2D([], [], marker="o", ls="", markerfacecolor="white", markeredgecolor="0.5",
                label="previous bundle: ONT_1kGP_1019 default"))
axes[1].legend(handles=h, frameon=False, loc="lower right", fontsize=6)
axes[0].set_title("Each library is scored against its own platform group", loc="left")
fig.savefig(out, dpi=300, bbox_inches="tight")
