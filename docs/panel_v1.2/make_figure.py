"""Figure for docs/panel_v1.2/README.md (fig_panel_v12_r2.png).

Inputs (this directory): pon_members.tsv (release cn/pon/pon_members.tsv),
release_chrX_segments_per_genome.tsv (chrX gain/loss segments per library in the
panel-v1.2 candidate and panel-v1.2-r2 releases, cn/panel_cn_segments.tsv),
lofo_per_genome.tsv (full-panel LOFO Tier 1 per genome with sex call).
Run: python3 make_figure.py   (matplotlib, pandas, numpy)
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

pm = pd.read_csv("pon_members.tsv", sep="\t")
seg = pd.read_csv("release_chrX_segments_per_genome.tsv", sep="\t")
A = pm.merge(seg, on="key", how="left").fillna({c: 0 for c in seg.columns if c != "key"})
G = pd.read_csv("lofo_per_genome.tsv", sep="\t")
cols = {"male": "#4C72B0", "female": "#C44E52", "ambiguous": "#7F7F7F"}
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})
fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.9), gridspec_kw=dict(width_ratios=[1.3, 0.9, 1.1], wspace=0.5))

# (a) chrX ratio of female and ambiguous libraries with the three candidate female upper bounds
ax = axs[0]
fx = pm[(pm.x_ratio > 0.75) & (pm.x_ratio < 1.32)]
bins = np.arange(0.78, 1.30, 0.01)
ax.hist(fx[fx.sex_call == "female"].x_ratio, bins=bins, color=cols["female"], alpha=.85, label="female call (n=806)")
ax.hist(fx[fx.sex_call == "ambiguous"].x_ratio, bins=bins, color=cols["ambiguous"], label="ambiguous call")
ax.set_ylim(0, 150)
for b, ls in ((1.12, "-"), (1.18, "--"), (1.25, ":")):
    ax.axvline(b, color="k", lw=.8, ls=ls)
    ax.text(b + 0.004, 148, f"{b:.2f}", fontsize=6, rotation=90, va="top")
labels = [("1kGP_ONT_1019_CHM13/NA18946", "NA18946", 1.135, 40, 1.5),
          ("1kGP_LRSC_ONT_CHM13/NA19238", "NA19238 LRSC", 1.20, 75, 1.5),
          ("1kGP_ONT_1019_CHM13/HG02805", "HG02805", 1.27, 40, 1.5),
          ("1kGP_ONT_1019_CHM13/NA19238", "NA19238\nONT-1019", 1.065, 105, 3)]
for key, text, xt, yt, y0 in labels:
    x = pm.loc[pm.key == key, "x_ratio"].iloc[0]
    ax.annotate(text, (x, y0), xytext=(xt, yt), fontsize=6, ha="center",
                arrowprops=dict(arrowstyle="-", lw=.5, color="0.3"))
ax.set_xlabel("chrX depth ratio"); ax.set_ylabel("libraries")
ax.set_title("Three libraries lie\nbetween the candidate bounds", loc="left")
ax.legend(frameon=False, loc="upper left", fontsize=6); ax.set_xlim(0.78, 1.30)

# (b) chrX gain segments per genome, pooled vs sex-matched chrX reference
ax = axs[1]; w = 0.38
for i, sx in enumerate(["male", "female"]):
    s = A[A.sex_call == sx]; m0, m1 = s.v12_chrX_gain.mean(), s.v12r2_chrX_gain.mean()
    ax.bar(i - w / 2, m0, w, color=cols[sx], alpha=.3, edgecolor=cols[sx], hatch="///", lw=.6)
    ax.bar(i + w / 2, m1, w, color=cols[sx], edgecolor=cols[sx])
    ax.text(i - w / 2, m0 + 0.08, f"{m0:.2g}", ha="center", fontsize=6)
    ax.text(i + w / 2, m1 + 0.08, f"{m1:.2g}", ha="center", fontsize=6)
ax.set_xticks([0, 1]); ax.set_xticklabels(["male\n(n=795)", "female\n(n=806)"])
ax.set_ylabel("chrX gain segments per genome"); ax.set_ylim(0, 5.2)
ax.legend(handles=[Patch(facecolor="0.85", hatch="///", edgecolor="0.4", label="pooled chrX (v1.2)"),
                   Patch(facecolor="0.4", label="sex-matched (r2)")], frameon=False, fontsize=6, loc="upper right")
ax.set_title("Male chrX gains\nfall 8-fold", loc="left")

# (c) healthy Tier 1 per genome (full-panel LOFO) by platform and sex call
ax = axs[2]
plats = ["ONT 1kGP-1019", "ONT LRSC R10", "ONT LRSC R9", "PacBio HiFi"]
for j, sx in enumerate(["male", "female"]):
    vals = [G[(G.group == p) & (G.sex_call == sx)].T1_v12r2.mean() for p in plats]
    ax.scatter(np.arange(4) + (j - 0.5) * 0.22, vals, color=cols[sx], s=18, label=sx, zorder=3)
allv = G.T1_v12r2.mean()
ax.axhline(allv, color="0.3", lw=.8); ax.text(3.5, allv + 0.006, f"all {allv:.3f}", fontsize=6, ha="right")
ax.axhline(0.301, color="k", lw=.8, ls="--"); ax.text(3.5, 0.307, "stop 0.301", fontsize=6, ha="right")
ax.set_xticks(range(4)); ax.set_xticklabels(["ONT\n1019", "LRSC\nR10", "LRSC\nR9", "PacBio\nHiFi"])
ax.set_ylabel("healthy Tier 1 per genome"); ax.set_ylim(0, 0.34); ax.set_xlim(-0.5, 3.6)
ax.legend(frameon=False, fontsize=6, loc="lower left")
ax.set_title("Every platform and sex\nis below the stop", loc="left")
for a, l in zip(axs, "abc"):
    a.text(-0.18, 1.08, l, transform=a.transAxes, fontweight="bold", fontsize=10)
fig.savefig("fig_panel_v12_r2.png", dpi=300, bbox_inches="tight")
