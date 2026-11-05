#!/usr/bin/env python3
"""Figures of docs/R10_EXCESS.md (R10-excess diagnosis of the panel-v1.1 candidate).

    figures_r10.py --diag-genomes diag_genomes.tsv --diag-loci diag_added_loci.tsv \
        --arms ARM=per_genome.tsv ... --noisy noisy.tsv --out-dir DIR

fig_r10_diagnosis.png   (a) SVIM-only translocation records with QUAL < 10 per genome
                        against coverage; (b) properties of the Tier-1 loci added by
                        panel-v1.1, 11 noisy genomes vs the other 1,608
fig_r10_candidates.png  healthy LOFO Tier-1 calls per genome by stratum for panel-v1.0,
                        panel-v1.1 and the candidate discovery rules
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

plt.rcParams.update({"font.size": 7, "axes.labelsize": 7, "axes.titlesize": 7.5, "xtick.labelsize": 6,
                     "ytick.labelsize": 6, "legend.fontsize": 6, "axes.spines.top": False,
                     "axes.spines.right": False, "font.family": "DejaVu Sans"})
ORDER = ["ONT 1kGP-1019", "ONT LRSC R9", "ONT LRSC R10", "PacBio HiFi"]
COLS = dict(zip(ORDER, ["#4C72B0", "#55A868", "#C44E52", "#8172B2"]))


def diagnosis(G, L, noisy, out):
    G = G.copy()
    G["noisy"] = [(a, b) in noisy for a, b in zip(G.cohort, G.sample_id)]
    G["svim_TRA_lowq"] = G.get("single_SVIM_TRA_q<10", 0)
    L = L.copy()
    L["noisy"] = [(a, b) in noisy for a, b in zip(L.cohort, L.sample_id)]
    num = lambda s: pd.to_numeric(s, errors="coerce")  # noqa: E731
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.2))
    ax = axs[0]
    for g in ORDER:
        d = G[(G.group == g) & ~G.noisy]
        ax.scatter(d.coverage_x, d.svim_TRA_lowq + 1, s=4, color=COLS[g], alpha=0.45, linewidths=0,
                   label=f"{g} (n={(G.group == g).sum()})")
        d = G[(G.group == g) & G.noisy]
        ax.scatter(d.coverage_x, d.svim_TRA_lowq + 1, s=22, facecolors="none", edgecolors=COLS[g], linewidths=1.0)
    ax.scatter([], [], s=22, facecolors="none", edgecolors="k", label="11 genomes with ≥10 Tier-1 calls in v1.0")
    ax.set_yscale("log")
    ax.set_yticks([1, 10, 100, 1000, 10000])
    ax.set_yticklabels(["0", "10", "100", "1k", "10k"])
    ax.set_xlabel("Genome coverage (x)")
    ax.set_ylabel("SVIM-only translocations, QUAL < 10")
    lg = ax.legend(frameon=False, fontsize=6, loc="lower right", markerscale=2.5)
    for h in lg.legend_handles:
        h.set_alpha(1)
    ax.set_title("Floor-quality SVIM translocations pile up in a few genomes", loc="left", fontsize=7.5)
    ax.margins(0.04)
    ax.text(-0.12, 1.04, "a", transform=ax.transAxes, fontweight="bold", fontsize=10)
    ax = axs[1]
    allsv = L.junctions_svim_only == L.N_Junctions
    feats = [("all junctions SVIM-only", lambda m: allsv[m].mean()),
             ("BND_cluster class", lambda m: (L.Cmplx_Class[m] == "BND_cluster").mean()),
             ("SVIM QUAL at floor (5)", lambda m: (num(L.median_svim_qual[m]) <= 5).mean()),
             ("alt reads < 0.15 x coverage", lambda m: (num(L.median_alt_over_cov[m]) < 0.15).mean()),
             ("no junction seen in panel", lambda m: (L.junctions_panel_carriers_0[m] == L.N_Junctions[m]).mean())]
    y = np.arange(len(feats))
    for m, col, off, lab in [(L.noisy, "#C44E52", -0.18, "11 noisy genomes"), (~L.noisy, "0.55", 0.18, "1,608 others")]:
        ax.barh(y + off, [f(m) for _, f in feats], height=0.34, color=col, label=lab)
    ax.set_yticks(y)
    ax.set_yticklabels([f for f, _ in feats])
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel(f"Share of added Tier-1 loci ({int(L.noisy.sum())} noisy, {int((~L.noisy).sum())} others)")
    ax.legend(frameon=False, fontsize=6, loc="lower right")
    ax.set_title("Their added Tier-1 loci are SVIM-only artefacts", loc="left", fontsize=7.5)
    ax.text(-0.12, 1.04, "b", transform=ax.transAxes, fontweight="bold", fontsize=10)
    fig.tight_layout()
    fig.savefig(out, dpi=300)


def candidates(arms, noisy, out):
    strata = [("all 1,619", lambda d: d.index == d.index), ("ONT LRSC R10 (247)", lambda d: d.group == "ONT LRSC R10"),
              ("1,608 others", lambda d: ~d.noisy), ("ONT 1kGP-1019 (1,017)", lambda d: d.group == "ONT 1kGP-1019")]
    fig, ax = plt.subplots(figsize=(7.2, 2.8))
    names = list(arms)
    w = 0.8 / len(names)
    cols = ["0.65", "#C44E52", "#4C72B0", "#55A868", "#8172B2"][:len(names)]
    for i, n in enumerate(names):
        d = arms[n]
        v = [d.loc[f(d), "T1"].mean() for _, f in strata]
        x = np.arange(len(strata)) + (i - (len(names) - 1) / 2) * w
        ax.bar(x, v, width=w * 0.95, color=cols[i], label=n)
        for xi, vi in zip(x, v):
            ax.text(xi, vi + 0.03, f"{vi:.2f}", ha="center", fontsize=5.5)
    ax.set_xticks(np.arange(len(strata)))
    ax.set_xticklabels([s for s, _ in strata])
    ax.set_ylabel("Tier-1 calls per healthy genome")
    ax.legend(frameon=False, fontsize=6, ncol=len(names), loc="lower center", bbox_to_anchor=(0.5, 1.0))
    ax.set_title("Leaving floor-quality SVIM breakends out of discovery removes the R10 excess",
                 loc="left", fontsize=7.5, pad=16)
    ax.margins(y=0.12)
    fig.tight_layout()
    fig.savefig(out, dpi=300)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--diag-genomes", required=True)
    ap.add_argument("--diag-loci", required=True)
    ap.add_argument("--arms", nargs="+", required=True, help="NAME=per_genome.tsv:T1column")
    ap.add_argument("--noisy", required=True, help="TSV cohort, sample_id of the 11 noisy genomes")
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    nz = pd.read_csv(a.noisy, sep="\t")
    noisy = set(zip(nz.cohort, nz.sample_id))
    diagnosis(pd.read_csv(a.diag_genomes, sep="\t"), pd.read_csv(a.diag_loci, sep="\t", low_memory=False), noisy,
              f"{a.out_dir}/fig_r10_diagnosis.png")
    arms = {}
    for spec in a.arms:
        name, rest = spec.split("=", 1)
        path, col = rest.rsplit(":", 1)
        d = pd.read_csv(path, sep="\t")
        d["T1"] = d[col]
        d["noisy"] = [(x, y) in noisy for x, y in zip(d.cohort, d.sample_id)]
        if "group" not in d:
            g = arms[next(iter(arms))][["cohort", "sample_id", "group"]]
            d = d.merge(g, on=["cohort", "sample_id"])
        arms[name] = d
    candidates(arms, noisy, f"{a.out_dir}/fig_r10_candidates.png")


if __name__ == "__main__":
    main()
