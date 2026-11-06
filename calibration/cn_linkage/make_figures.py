#!/usr/bin/env python3
"""Figures for docs/CN_LINKAGE_PHASE1.md (roadmap item 15, Phase 1).

  make_figures.py --test-a literature_cn_summary.tsv --test-d testD_per_genome.tsv \\
      --genomes genomes.tsv --grid grid_results.tsv --out docs/cn_linkage

  fig_testA_literature.png   literature events: best tier, junctions only vs with CN linkage
  fig_testD_healthy.png      healthy Tier 1 per genome (LOFO) per stratum vs the stop rule
  fig_testB_sensitivity.png  simulation grid: detection and Tier 1 by coverage, platform and
                             triplicated-segment size
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "xtick.labelsize": 6,
                     "ytick.labelsize": 6, "legend.fontsize": 7, "axes.spines.top": False,
                     "axes.spines.right": False, "savefig.bbox": "tight"})
ORDER = ["Tier1_high_priority", "Tier2", "Tier3", "Benign_likely", "not_called"]
LAB = {"Tier1_high_priority": "Tier 1", "Tier2": "Tier 2", "Tier3": "Tier 3", "Benign_likely": "Likely benign",
       "not_called": "Not called"}
COL = {"Tier1_high_priority": "#1f4e79", "Tier2": "#4f81bd", "Tier3": "#9dc3e6", "Benign_likely": "#d9d9d9",
       "not_called": "white"}


def fig_a(path, out):
    A = pd.read_csv(path, sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.6))
    for ax, grp in zip(axes, (True, False)):
        g = A[A.dup_trp == grp]
        n_cn = int(g.cn_called.sum()); n_main = int(g.main_called.sum())
        title = (f"DUP-TRP/INV-DUP events (n={len(g)}):\n{n_cn} called with CN linkage, {n_main} without" if grp else
                 f"Other junction-level events (n={len(g)}):\n{n_cn} called with CN linkage, {n_main} without")
        handles = {}
        for yi, col in enumerate(("main_best_tier", "cn_best_tier")):
            left = 0
            vc = g[col].fillna("not_called").value_counts()
            for t in ORDER:
                n = int(vc.get(t, 0))
                if not n:
                    continue
                h = ax.barh(yi, n, left=left, color=COL[t], edgecolor="#404040", linewidth=0.5, height=0.6)
                handles[t] = h
                if n >= 2:
                    ax.text(left + n / 2, yi, str(n), ha="center", va="center", fontsize=6,
                            color="white" if t in ORDER[:2] else "black")
                left += n
        ax.set_yticks([0, 1])
        ax.set_yticklabels(["Junctions only\n(main)", "With CN linkage"] if grp else ["", ""])
        ax.set_xlim(0, len(g)); ax.set_xlabel("Events"); ax.set_title(title, loc="left")
    fig.legend([plt.Rectangle((0, 0), 1, 1, fc=COL[t], ec="#404040", lw=0.5) for t in ORDER],
               [LAB[t] for t in ORDER], loc="lower center", ncol=5, frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.1, 1, 1))
    fig.savefig(f"{out}/fig_testA_literature.png", dpi=300)


def fig_d(path, genomes, out):
    D = pd.read_csv(path, sep="\t")
    G = pd.read_csv(genomes, sep="\t")
    D["stratum"] = (D.cohort + "/" + D["sample"]).map(dict(zip(G.key, G.stratum)))
    S = D.groupby("stratum").agg(n=("sample", "size"), old=("tier1_old", "mean"), new=("tier1_new", "mean")).reset_index()
    S = pd.concat([S, pd.DataFrame([dict(stratum=f"All {len(D)}", n=len(D), old=D.tier1_old.mean(), new=D.tier1_new.mean())])])
    x = np.arange(len(S)); w = 0.36
    fig, ax = plt.subplots(figsize=(4.6, 2.6))
    ax.bar(x - w / 2, S.old, w, color="#bfbfbf", edgecolor="#404040", lw=0.5, label="Panel-v1.1 (junctions only)")
    ax.bar(x + w / 2, S.new, w, color="#1f4e79", edgecolor="#404040", lw=0.5, label="With CN-linked loci")
    lim = S.old.iloc[-1] + 0.05
    ax.plot([x[-1] - 0.45, x[-1] + 0.45], [lim] * 2, color="#c00000", lw=1, ls="--")
    ax.text(x[-1] + 0.5, lim, "stop rule\n(+0.05)", va="center", ha="left", fontsize=6, color="#c00000")
    for xi, (o, nw) in enumerate(zip(S.old, S.new)):
        ax.text(xi + w / 2, (nw + 0.02) if xi < len(S) - 1 else lim + 0.03, f"+{nw - o:.2f}", ha="center", va="bottom", fontsize=6)
    ax.set_xticks(x); ax.set_xticklabels([f"{s}\n(n={n})" for s, n in zip(S.stratum, S.n)])
    ax.set_ylabel("Tier 1 per healthy genome (LOFO)")
    ax.set_xlim(-0.6, len(S) + 0.2); ax.set_ylim(0, max(1.0, S.new.max() * 1.2))
    d = S.new.iloc[-1] - S.old.iloc[-1]
    per = S.iloc[:-1].assign(delta=lambda t: t.new - t.old)
    top = per.loc[per.delta.idxmax()]
    note = (f"\n{top.stratum} {top.delta:+.2f} = {int(round(top.delta * top.n))} locus in {int(top.n)} genomes"
            if top.delta > d else "")
    ax.set_title(f"CN linkage adds {d:+.2f} Tier 1 per healthy genome overall{note}", loc="left")
    ax.legend(frameon=False, loc="upper left", fontsize=6)
    fig.tight_layout(); fig.savefig(f"{out}/fig_testD_healthy.png", dpi=300)


def fig_b(path, out):
    R = pd.read_csv(path, sep="\t")
    for c in ("any_called", "cn_called"):
        R[c] = R[c].astype(str).str.lower().eq("true")
    R["tier1"] = R.best_tier.astype(str).str.startswith("Tier1")
    plats = [p for p in ("R10", "R9", "PacBio") if p in set(R.platform)]
    trips = sorted(R.trip.unique())
    cmap = plt.get_cmap("Blues")
    tc = {t: cmap(0.35 + 0.6 * i / max(1, len(trips) - 1)) for i, t in enumerate(trips)}
    fig, axes = plt.subplots(2, len(plats), figsize=(7.0, 4.0), sharex=True, sharey=True)
    names = {"R10": "ONT R10", "R9": "ONT R9", "PacBio": "PacBio HiFi"}
    for j, p in enumerate(plats):
        for i, (col, ylab) in enumerate((("any_called", "Detected as a CxSV locus"), ("tier1", "Tier 1"))):
            ax = axes[i, j]
            for t in trips:
                g = R[(R.platform == p) & (R.trip == t)].groupby("cov")[col].agg(["mean", "size"]).reset_index()
                ax.plot(g["cov"], g["mean"], "-o", color=tc[t], ms=3, lw=1.2, label=f"{int(t) // 1000} kb")
            ax.axhline(0.9, color="#c00000", lw=0.8, ls="--")
            ax.set_ylim(-0.04, 1.04); ax.set_xticks(sorted(R["cov"].unique()))
            if i == 0:
                ax.set_title(names.get(p, p), loc="left")
            if j == 0:
                ax.set_ylabel(f"{ylab}\n(fraction of cells)")
            if i == 1:
                ax.set_xlabel("Coverage (x)")
    axes[0, 0].text(10, 0.92, "criterion 0.9", fontsize=6, color="#c00000", va="bottom")
    axes[0, -1].legend(title="Triplicated\nsegment", frameon=False, loc="lower right", fontsize=6, title_fontsize=6)
    fig.tight_layout(); fig.savefig(f"{out}/fig_testB_sensitivity.png", dpi=300)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--test-a"); ap.add_argument("--test-d"); ap.add_argument("--genomes"); ap.add_argument("--grid")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    if a.test_a:
        fig_a(a.test_a, a.out)
    if a.test_d:
        fig_d(a.test_d, a.genomes, a.out)
    if a.grid:
        fig_b(a.grid, a.out)


if __name__ == "__main__":
    main()
