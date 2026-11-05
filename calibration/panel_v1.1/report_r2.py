#!/usr/bin/env python3
"""panel-v1.1 (r2) release record: stratum table, residual noisy genomes and the
summary figure (docs/panel_v1.1/README.md).

    report_r2.py --per-genome lofo_r2_per_genome.tsv --compare compare_v10_r2.tsv \
        --candidate-rates candidate_rates.tsv --noisy noisy_genomes.tsv --out-dir DIR

--per-genome   panel_rates_collect.py per-genome table of the r2 LOFO run
               (columns *_v1.0 = panel-v1.0 item-10 calibration, *_v1.1 = r2)
--compare      compare_builds.py panel-v1.0 vs r2 (blind span per genome)
--candidate-rates  docs/r10_excess/candidate_rates.tsv (v1.0 and v1.1-candidate rows)
--noisy        docs/r10_excess/noisy_genomes.tsv (the 11 genomes with >= 10 Tier 1
               in panel-v1.0, keyed by cohort + sample)
Writes DIR/lofo_r2_strata_summary.tsv, DIR/lofo_r2_residual_noisy.tsv,
DIR/lofo_r2_noisy11.tsv and DIR/fig_panel_v11_r2.png.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

GROUPS = [("all", None), ("ONT LRSC R10", "ONT LRSC R10"), ("ONT LRSC R9", "ONT LRSC R9"),
          ("ONT 1kGP-1019", "ONT 1kGP-1019"), ("PacBio HiFi", "PacBio HiFi")]
GREY, GOLD, BLUE = "#9a9a9a", "#d9a441", "#1f5fa8"
STOP = 0.457


def strata(g, noisy):
    rows = []
    masks = [(n, (g.group == v) if v else pd.Series(True, index=g.index)) for n, v in GROUPS]
    masks += [("1,608 others", ~noisy), ("noisy 11", noisy)]
    for name, m in masks:
        s = g[m]
        rows.append(dict(stratum=name, n=len(s), T1_v10=round(s["T1_v1.0"].mean(), 3),
                         T1_r2=round(s["T1_v1.1"].mean(), 3), T2_v10=round(s["T2_v1.0"].mean(), 2),
                         T2_r2=round(s["T2_v1.1"].mean(), 2), share_T1_r2=round((s["T1_v1.1"] > 0).mean(), 3),
                         loci_r2=round(s["loci_v1.1"].mean(), 1)))
    return pd.DataFrame(rows)


def figure(g, noisy, st, cr, cmp_, out):
    plt.rcParams.update({"font.size": 7, "axes.titlesize": 8, "axes.labelsize": 8, "xtick.labelsize": 6,
                         "ytick.labelsize": 6, "legend.fontsize": 6, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.titlelocation": "left"})
    cols = ["T1/genome all", "T1/genome R10", "T1/genome R9", "T1/genome ONT-1019", "T1/genome PacBio",
            "T1/genome 1,608 others"]
    v10 = cr.loc["panel-v1.0", cols].values
    v11 = cr.loc["panel-v1.1 as built", cols].values
    r2 = st.set_index("stratum").loc[["all", "ONT LRSC R10", "ONT LRSC R9", "ONT 1kGP-1019", "PacBio HiFi",
                                      "1,608 others"], "T1_r2"].values
    fig = plt.figure(figsize=(7.6, 3.0))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.45, 1, 1], wspace=0.5)
    ax = fig.add_subplot(gs[0])
    x, w = np.arange(len(cols)), 0.27
    ax.bar(x - w, v10, w, color=GREY, label="panel-v1.0")
    ax.bar(x, v11, w, color=GOLD, label="v1.1 candidate (64e10c1b1363)")
    ax.bar(x + w, r2, w, color=BLUE, label="panel-v1.1 r2 (c31ecc9f8fa2)")
    ax.axhline(STOP, ls="--", lw=0.8, color="k")
    ax.text(5.5, STOP + 0.06, f"stop limit {STOP}", ha="right", va="bottom", fontsize=6)
    ax.set_xticks(x)
    ax.set_xticklabels(["All", "R10", "R9", "ONT\n1019", "Pac\nBio", "Others"])
    ax.set_ylabel("Healthy Tier-1 calls per genome")
    ax.set_title("svim_lowq: Tier 1 below v1.0")
    ax.legend(frameon=False, loc="upper right")
    ax2 = fig.add_subplot(gs[1])
    res = g["T1_v1.1"] >= 10
    ax2.scatter(g.loc[~noisy & ~res, "T1_v1.0"], g.loc[~noisy & ~res, "T1_v1.1"], s=6, color=GREY, alpha=0.5,
                lw=0, label="other genomes")
    ax2.scatter(g.loc[noisy, "T1_v1.0"], g.loc[noisy, "T1_v1.1"], s=14, color=BLUE, label="11 noisy genomes")
    ax2.scatter(g.loc[res & ~noisy, "T1_v1.0"], g.loc[res & ~noisy, "T1_v1.1"], s=14, facecolor="white",
                edgecolor=BLUE, lw=0.9, label="newly noisy")
    for setter in (ax2.set_xscale, ax2.set_yscale):
        setter("symlog", linthresh=1)
    ax2.plot([0, 120], [0, 120], lw=0.6, color="k", ls=":")
    ax2.set_xlim(-0.3, 150)
    ax2.set_ylim(-0.3, 150)
    for ticks, labels in ((ax2.set_xticks, ax2.set_xticklabels), (ax2.set_yticks, ax2.set_yticklabels)):
        ticks([0, 1, 10, 100])
        labels(["0", "1", "10", "100"])
    offsets = [(-34, 4), (4, -12), (-38, -2)]
    for (_, rr), off in zip(g[res].sort_values("T1_v1.1", ascending=False).iterrows(), offsets):
        ax2.annotate(rr["sample_id"], (rr["T1_v1.0"], rr["T1_v1.1"]), xytext=off, textcoords="offset points",
                     fontsize=6)
    ax2.set_xlabel("Tier 1, panel-v1.0")
    ax2.set_ylabel("Tier 1, panel-v1.1 r2")
    ax2.set_title("Noisy genomes drop")
    ax2.legend(frameon=False, loc="lower right")
    ax3 = fig.add_subplot(gs[2])
    for col, lab, c in (("blind_span_v1.0", "panel-v1.0", GREY), ("blind_span_v1.1", "panel-v1.1 r2", BLUE)):
        v = np.sort(cmp_[col].values / 1e6)
        ax3.step(np.maximum(v, 0.01), np.arange(1, len(v) + 1) / len(v), color=c, lw=1.2, label=lab, where="post")
    ax3.set_xscale("log")
    ax3.set_xticks([0.01, 1, 100])
    ax3.set_xticklabels(["≤0.01", "1", "100"])
    ax3.set_xlabel("Oversized span per genome (Mb)")
    ax3.set_ylabel("Fraction of genomes")
    med0, med1 = cmp_["blind_span_v1.0"].median() / 1e6, cmp_["blind_span_v1.1"].median() / 1e6
    ax3.set_title(f"Blind span {med0:.0f} → {med1:.1f} Mb")
    ax3.legend(frameon=False, loc="center right")
    for a, l in zip((ax, ax2, ax3), "abc"):
        a.text(-0.18, 1.06, l, transform=a.transAxes, fontsize=9, fontweight="bold", va="bottom")
    fig.savefig(out, dpi=300, bbox_inches="tight")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("per-genome", "compare", "candidate-rates", "noisy", "out-dir"):
        ap.add_argument(f"--{k}", required=True)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    g = pd.read_csv(a.per_genome, sep="\t")
    nz = pd.read_csv(a.noisy, sep="\t")
    noisy = (g.cohort + ":" + g.sample_id).isin(set(nz.cohort + ":" + nz.sample_id))
    assert noisy.sum() == len(nz), "noisy genomes not all found"
    st = strata(g, noisy)
    st.to_csv(f"{a.out_dir}/lofo_r2_strata_summary.tsv", sep="\t", index=False)
    cols = ["cohort", "sample_id", "group", "coverage_x", "T1_v1.0", "T1_v1.1", "T2_v1.1"]
    g[g["T1_v1.1"] >= 10].sort_values("T1_v1.1", ascending=False)[cols].to_csv(
        f"{a.out_dir}/lofo_r2_residual_noisy.tsv", sep="\t", index=False)
    g[noisy].sort_values("T1_v1.1", ascending=False)[cols].to_csv(f"{a.out_dir}/lofo_r2_noisy11.tsv", sep="\t",
                                                                   index=False)
    cr = pd.read_csv(a.candidate_rates, sep="\t").set_index("arm")
    figure(g, noisy, st, cr, pd.read_csv(a.compare, sep="\t"), f"{a.out_dir}/fig_panel_v11_r2.png")
    print(st.to_string(index=False))


if __name__ == "__main__":
    main()
