#!/usr/bin/env python3
"""Per-platform calibration of the library QC flag (docs/LIBRARY_QC.md).

    calibrate.py --survey survey_all.tsv --groups lofo_r2_per_genome.tsv \
        --noisy docs/r10_excess/noisy_genomes.tsv --k 6 --out-dir OUT [--config-out config/library_qc.yaml]

--survey  one library_qc.py row per genome (sample = <cohort>__<sample_id>; extra__* rows
          are patients, scored but not used for calibration)
--groups  per-genome table with cohort, sample_id, group (platform stratum), coverage_x,
          T1_v1.0, T1_v1.1 (panel_rates_collect.py per-genome output)
For each group and flag metric m (svim_only_bnd_per_gb, svim_only_any_per_gb):
  y = log10(m + eps), x = log10(coverage_est)
  b = Theil-Sen slope, a = median(y - b x), r = y - a - b x, r0 = median(r), mad = 1.4826 MAD(r)
  z = (r - r0) / mad;  FLAG when any flag metric has z > k.
Writes OUT/library_qc_calibration.tsv (every genome: metrics, z, flag),
OUT/library_qc_flagged.tsv, OUT/library_qc_coverage_independence.tsv,
OUT/fig_library_qc.png and, with --config-out, the calibration block of config/library_qc.yaml.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr, theilslopes  # noqa: E402

FLAG_METRICS = ("svim_only_bnd_per_gb", "svim_only_any_per_gb")
READ_METRICS = ("frac_sa", "frac_sa_interchrom")
EPS = 1e-3
GROUP_KEYS = {"ONT 1kGP-1019": "ONT_1kGP_1019", "ONT LRSC R10": "ONT_LRSC_R10", "ONT LRSC R9": "ONT_LRSC_R9",
              "PacBio HiFi": "PacBio_HiFi"}


def fit(df, m):
    y, x = np.log10(df[m] + EPS), np.log10(df["coverage_est"])
    b = theilslopes(y, x)[0]
    a = float(np.median(y - b * x))
    r = y - a - b * x
    r0 = float(np.median(r))
    mad = float(1.4826 * np.median(np.abs(r - r0)))
    return {"a": round(a, 6), "b": round(float(b), 6), "r0": round(r0, 6), "mad": round(mad, 6)}


def z_of(df, m, c):
    return (np.log10(df[m] + EPS) - (c["a"] + c["b"] * np.log10(df["coverage_est"])) - c["r0"]) / c["mad"]


def write_config(path, cal, k, n):
    L = ["# Library QC (docs/LIBRARY_QC.md). NOT fingerprinted: QC never changes which loci",
         "# a genome has. Calibration written by calibration/library_qc/calibrate.py from the",
         f"# panel-v1.1 genomes ({n}); z = (log10(m + eps) - (a + b log10(coverage_est)) - r0) / mad.",
         "library_qc:",
         "  enabled: true",
         "  windows: config/library_qc_windows.hs1.bed",
         "  min_mapq: 20",
         f"  k: {k:g}",
         f"  eps: {EPS:g}",
         "  # calibration group of a patient (lib/cxsv_libgroup.py): sample-sheet qc_group when set,",
         "  # else detected from the alignment header, else from the sampled reads, else by platform",
         "  default_group: {ONT: ONT_1kGP_1019, PACBIO: PacBio_HiFi}",
         "  group_detection: {max_reads: 20000, min_reads: 200, r10_min_median_q: 30, r9_max_median_q: 28,",
         "                    uniform_q_min_fraction: 0.9, pacbio_min_fraction: 0.5}",
         "  # conservative mode for FLAGged libraries (on by default, PI decision 2026-11-05): a",
         "  # Tier-1 locus stays Tier 1 only with > min_multi_caller_fraction of its junctions",
         "  # called by >= 2 callers or >= min_cn_changed copy-number-changed segments;",
         "  # otherwise it is reported as Tier 2 (score and points unchanged)",
         "  conservative_mode: true",
         "  conservative: {min_multi_caller_fraction: 0.5, min_cn_changed: 1, demote_to: Tier2}",
         "  groups:"]
    for g, ms in cal.items():
        L.append(f"    {g}:")
        for m, c in ms.items():
            L.append(f"      {m}: {{a: {c['a']}, b: {c['b']}, r0: {c['r0']}, mad: {c['mad']}}}")
    with open(path, "w") as fh:
        fh.write("\n".join(L) + "\n")


def figure(d, cal, k, out):
    plt.rcParams.update({"font.size": 7, "axes.titlesize": 8, "axes.labelsize": 8, "xtick.labelsize": 6,
                         "ytick.labelsize": 6, "legend.fontsize": 6, "axes.spines.top": False,
                         "axes.spines.right": False, "axes.titlelocation": "left"})
    order = ["ONT 1kGP-1019", "ONT LRSC R10", "ONT LRSC R9", "PacBio HiFi"]
    labels = ["ONT\n1019", "LRSC\nR10", "LRSC\nR9", "PacBio\nHiFi"]
    GREY, BLUE, RED = "#9a9a9a", "#1f5fa8", "#c4522f"
    fig, axs = plt.subplots(1, 3, figsize=(7.6, 2.8), gridspec_kw={"wspace": 0.42})
    rng = np.random.default_rng(0)
    titles = ("SVIM-only BND/TRA per Gb", "SVIM-only records per Gb")
    for ax, m, t in zip(axs[:2], FLAG_METRICS, titles):
        for i, gname in enumerate(order):
            x = d[d.group == gname]
            jit = i + rng.uniform(-0.28, 0.28, len(x))
            zz = x[f"z_{m}"].clip(-6, 30)
            norm = ~x.flagged & ~x.interest
            ax.scatter(jit[norm.values], zz[norm], s=3, color=GREY, alpha=0.5, lw=0)
            fl = x.flagged & ~x.interest
            ax.scatter(jit[fl.values], zz[fl], s=12, facecolor="white", edgecolor=RED, lw=0.8)
            it = x.interest
            ax.scatter(jit[it.values], zz[it], s=12, color=BLUE, marker="D", lw=0,
                       edgecolor=RED)
            itf = x.interest & x.flagged
            ax.scatter(jit[itf.values], zz[itf], s=26, facecolor="none", edgecolor=RED, lw=0.8)
        ax.axhline(k, ls="--", lw=0.8, color="k")
        ax.text(3.45, k + 0.6, f"z = {k:g}", ha="right", fontsize=6)
        ax.set_xticks(range(4))
        ax.set_xticklabels(labels)
        ax.set_ylabel("Coverage-adjusted robust z")
        ax.set_title(t)
    h = [plt.Line2D([], [], ls="", marker="o", color=GREY, ms=3, label="panel genome"),
         plt.Line2D([], [], ls="", marker="o", mfc="white", mec=RED, ms=4, label="flagged"),
         plt.Line2D([], [], ls="", marker="D", color=BLUE, ms=4, label="v1.0-noisy or HG01867"),
         plt.Line2D([], [], ls="", marker="o", mfc="none", mec=RED, ms=6, label="… and flagged")]
    axs[2].legend(handles=h, frameon=False, loc="upper left")
    ax = axs[2]
    zmax = d[[f"z_{m}" for m in FLAG_METRICS]].max(axis=1).clip(-6, 30)
    ax.scatter(d.coverage_est, zmax, s=3, color=GREY, alpha=0.5, lw=0)
    ax.scatter(d.coverage_est[d.flagged], zmax[d.flagged], s=12, facecolor="white", edgecolor=RED, lw=0.8)
    ax.axhline(k, ls="--", lw=0.8, color="k")
    rho = spearmanr(d.coverage_est, zmax).correlation
    ax.set_xscale("log")
    ax.set_xticks([5, 10, 20, 50, 100])
    ax.set_xticklabels(["5", "10", "20", "50", "100"])
    ax.set_xlabel("Coverage (x, from the QC windows)")
    ax.set_ylabel("Larger of the two z")
    ax.set_title(f"Independent of coverage (ρ = {rho:.2f})")
    for a, l in zip(axs, "abc"):
        a.text(-0.2, 1.06, l, transform=a.transAxes, fontsize=9, fontweight="bold", va="bottom")
    fig.savefig(out, dpi=300, bbox_inches="tight")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("survey", "groups", "noisy", "out-dir"):
        ap.add_argument(f"--{k}", required=True)
    ap.add_argument("--k", type=float, default=6)
    ap.add_argument("--extra-interest", default="1kGP_LRSC_ONT_CHM13__HG01867")
    ap.add_argument("--config-out", default=None)
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    s = pd.read_csv(a.survey, sep="\t")
    g = pd.read_csv(a.groups, sep="\t")
    g["sample"] = g.cohort + "__" + g.sample_id
    nz = pd.read_csv(a.noisy, sep="\t")
    interest = set(nz.cohort + "__" + nz.sample_id) | set(a.extra_interest.split(","))
    d = s.merge(g[["sample", "group", "coverage_x", "T1_v1.0", "T1_v1.1"]], on="sample", how="left")
    panel = d[d.group.notna()].copy()
    cal = {}
    for gname, x in panel.groupby("group"):
        cal[GROUP_KEYS[gname]] = {m: fit(x, m) for m in FLAG_METRICS}
    for m in FLAG_METRICS:
        d[f"z_{m}"] = np.nan
        for gname, key in GROUP_KEYS.items():
            idx = d.group == gname
            d.loc[idx, f"z_{m}"] = z_of(d[idx], m, cal[key][m])
    d["flagged"] = (d[[f"z_{m}" for m in FLAG_METRICS]] > a.k).any(axis=1)
    d["interest"] = d["sample"].isin(interest)
    p = d[d.group.notna()]
    d.to_csv(f"{a.out_dir}/library_qc_calibration.tsv", sep="\t", index=False, float_format="%.6g")
    p[p.flagged | p.interest].sort_values(["group", "sample"]).to_csv(
        f"{a.out_dir}/library_qc_flagged.tsv", sep="\t", index=False, float_format="%.6g")
    rows = []
    for gname, x in p.groupby("group"):
        for m in FLAG_METRICS + READ_METRICS + tuple(f"z_{m}" for m in FLAG_METRICS):
            rows.append(dict(group=gname, metric=m, n=len(x),
                             spearman_vs_coverage=round(spearmanr(x[m], x.coverage_est).correlation, 3)))
    pd.DataFrame(rows).to_csv(f"{a.out_dir}/library_qc_coverage_independence.tsv", sep="\t", index=False)
    figure(p, cal, a.k, f"{a.out_dir}/fig_library_qc.png")
    if a.config_out:
        write_config(a.config_out, cal, a.k, len(p))
    print(f"panel {len(p)}; flagged {int(p.flagged.sum())}; of interest {int(p.interest.sum())}, flagged "
          f"{int((p.flagged & p.interest).sum())}; others flagged {int((p.flagged & ~p.interest).sum())}")
    print(d[~d.group.notna()][["sample", "coverage_est", "frac_sa"] + list(FLAG_METRICS)].to_string(index=False))


if __name__ == "__main__":
    main()
