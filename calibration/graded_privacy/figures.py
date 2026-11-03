#!/usr/bin/env python3
"""Figures of the graded-privacy calibration (item 10), from sweep.py outputs.

  fig_tier1_vs_parameter.png   Tier-1 calls per healthy genome (LOFO, panel-only
                               privacy) against the ramp start f_lo (one line per
                               f_hi) and the power exponent k, with the
                               chance-demotion risk of private Tier-1 loci on the
                               second row
  fig_matched_fraction.png     matched junction fraction of shared_junctions
                               loci (LOFO) by complexity class, and of the
                               shifted (null) loci
  fig_before_after.png         Tier-1 per genome before / after by platform and
                               coverage stratum (LOFO and LOIO)
Usage: figures.py SWEEP_DIR LOCI_TSV_GZ NULL_TSV_GZ OUTDIR
"""
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "xtick.labelsize": 6,
                     "ytick.labelsize": 6, "legend.fontsize": 7, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.titlelocation": "left", "savefig.dpi": 300,
                     "font.family": "DejaVu Sans"})
FOCAL = "#c0392b"
GREY = "#8a8a8a"


def panel_letter(ax, s):
    ax.text(-0.12, 1.06, s, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")


def fig_curve(sweep, chosen, out):
    c = sweep[(sweep["mode"] == "LOFO") & (sweep["fraction"] == chosen["fraction"])]
    before = c[c["shape"] == "none"]["panel_only|all|T1_mean"].iloc[0]
    ramp = c[c["shape"] == "ramp"]
    pw = c[c["shape"] == "power"]
    g = chosen["graded"]
    fig = plt.figure(figsize=(7.2, 4.8))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.25, 0.04, 1], height_ratios=[2, 1], wspace=0.35, hspace=0.12,
                          left=0.11, right=0.98, top=0.93, bottom=0.11)
    axes = np.empty((2, 2), dtype=object)
    axes[0, 0] = fig.add_subplot(gs[0, 0]); axes[1, 0] = fig.add_subplot(gs[1, 0], sharex=axes[0, 0])
    axes[0, 1] = fig.add_subplot(gs[0, 2]); axes[1, 1] = fig.add_subplot(gs[1, 2], sharex=axes[0, 1])
    cax = fig.add_subplot(gs[0, 1])
    plt.setp(axes[0, 0].get_xticklabels(), visible=False); plt.setp(axes[0, 1].get_xticklabels(), visible=False)
    his = sorted(ramp["f_hi"].unique())
    cmap = plt.get_cmap("viridis")
    for i, hi in enumerate(his):
        r = ramp[ramp["f_hi"] == hi].sort_values("f_lo")
        if len(r) < 2:
            continue
        col = cmap(i / max(1, len(his) - 1))
        focal = g["shape"] == "ramp" and abs(hi - g.get("f_hi", -1)) < 1e-9
        kw = dict(color=FOCAL if focal else col, lw=1.6 if focal else 0.7, alpha=1 if focal else 0.6, zorder=3 if focal else 2)
        axes[0, 0].plot(r["f_lo"], r["panel_only|all|T1_mean"], **kw)
        axes[1, 0].plot(r["f_lo"], 100 * r["demotion_risk"], **kw)
    risk0 = c[c["shape"] == "none"]["demotion_risk"].iloc[0]
    axes[0, 0].axhline(before, color=GREY, lw=0.8, ls="--")
    for ax in axes[1, :]:
        ax.axhline(100 * risk0, color=GREY, lw=0.8, ls="--")
    axes[0, 0].text(0.0, before, "current rule (no penalty)", color=GREY, fontsize=7, va="top", ha="left")
    pw = pw.sort_values("k")
    axes[0, 1].plot(pw["k"], pw["panel_only|all|T1_mean"], "o-", color="#34495e", ms=3, lw=1)
    axes[1, 1].plot(pw["k"], 100 * pw["demotion_risk"], "o-", color="#34495e", ms=3, lw=1)
    axes[0, 1].axhline(before, color=GREY, lw=0.8, ls="--")
    axes[0, 1].set_xscale("log", base=2)
    axes[1, 1].set_xscale("log", base=2)
    lim = 100 * (risk0 + chosen["selection"]["max_risk_increase"])
    for ax in axes[1, :]:
        ax.axhline(lim, color=FOCAL, lw=0.8, ls=":")
    axes[1, 0].text(0.0, lim, "limit: current + 1 point", color=FOCAL, fontsize=6, ha="left", va="bottom")
    axes[1, 0].text(0.0, 100 * risk0, "current rule", color=GREY, fontsize=6, ha="left", va="bottom")
    if g["shape"] == "ramp":
        sel = ramp[(abs(ramp["f_lo"] - g["f_lo"]) < 1e-9) & (abs(ramp["f_hi"] - g["f_hi"]) < 1e-9)].iloc[0]
        axes[0, 0].plot(sel["f_lo"], sel["panel_only|all|T1_mean"], "o", color=FOCAL, ms=5, zorder=4)
        axes[1, 0].plot(sel["f_lo"], 100 * sel["demotion_risk"], "o", color=FOCAL, ms=5, zorder=4)
        axes[0, 0].annotate(f"chosen: ramp {g['f_lo']:.2f}-{g['f_hi']:.2f}\n{sel['panel_only|all|T1_mean']:.2f} per genome",
                            (sel["f_lo"], sel["panel_only|all|T1_mean"]), xytext=(-110, 22), textcoords="offset points",
                            fontsize=7, color=FOCAL, arrowprops=dict(arrowstyle="-", color=FOCAL, lw=0.6))
    axes[0, 0].set_title("Linear ramp: penalty starts at f_lo, full at f_hi (colour)")
    axes[0, 1].set_title("Power rule: penalty share f^k")
    axes[0, 0].set_ylabel("Tier-1 calls per healthy genome\n(leave-one-family-out)")
    axes[1, 0].set_ylabel("Chance demotion of private\nTier-1 loci (%)")
    axes[1, 0].set_xlabel("Ramp start f_lo (matched junction fraction)")
    axes[1, 1].set_xlabel("Exponent k")
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(min(his), max(his)))
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label("f_hi", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    for ax, s in zip(axes.flat, "abcd"):
        ax.margins(x=0.03)
    panel_letter(axes[0, 0], "a")
    panel_letter(axes[0, 1], "b")
    fig.savefig(out)
    plt.close(fig)


def fig_fraction(loci, null, out):
    d = loci[(loci["mode"] == "LOFO") & (loci["Panel_Status"] == "shared_junctions")]
    order = d["Cmplx_Class"].value_counts()
    order = order[order >= 200].index.tolist()
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.0), gridspec_kw={"width_ratios": [2.2, 1]})
    data = [d.loc[d["Cmplx_Class"] == c, "Panel_Matched_Fraction"].astype(float).to_numpy() for c in order]
    parts = axes[0].violinplot(data, positions=range(len(order)), showmedians=True, widths=0.8)
    for b in parts["bodies"]:
        b.set_facecolor("#5d8aa8"); b.set_alpha(0.6)
    axes[0].set_xticks(range(len(order)))
    axes[0].set_xticklabels([f"{c}\n(n={len(x):,})" for c, x in zip(order, data)], rotation=45, ha="right")
    axes[0].set_ylabel("Matched junction fraction")
    axes[0].set_title("Partly matched (shared_junctions) loci in healthy genomes")
    nf = null.loc[null["null_status"] != "unshiftable", "null_matched_fraction"].astype(float)
    bins = np.linspace(0, 1, 21)
    real = d["Panel_Matched_Fraction"].astype(float)
    axes[1].hist(real, bins=bins, color="#5d8aa8", alpha=0.7, density=True, label=f"real loci (n={len(real):,})")
    axes[1].hist(nf[nf > 0], bins=bins, color=GREY, alpha=0.7, density=True,
                 label=f"shifted loci with f>0 (n={int((nf > 0).sum()):,})")
    axes[1].set_xlabel("Matched junction fraction")
    axes[1].set_ylabel("Density")
    axes[1].set_title(f"{100 * (nf > 0).mean():.0f}% of shifted loci match by chance")
    axes[1].legend(frameon=False, loc="upper left", fontsize=6)
    panel_letter(axes[0], "a")
    panel_letter(axes[1], "b")
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def fig_before_after(per_ind, out, lowcov):
    p = per_ind.copy()
    p["cov_bin"] = pd.cut(p["coverage_x"].astype(float), [0, lowcov, 10, 20, 30, 1000],
                          labels=[f"<{lowcov:g}x", f"{lowcov:g}-10x", "10-20x", "20-30x", ">=30x"], right=False)
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.6), sharey=True, gridspec_kw={"wspace": 0.06})
    for ax, key, title in ((axes[0], "platform_group", "By platform"), (axes[1], "cov_bin", "By coverage")):
        groups = [g for g in (p[key].cat.categories if hasattr(p[key], "cat") else sorted(p[key].unique()))
                  if (p[key] == g).any()]
        x = np.arange(len(groups))
        for j, (mode, col) in enumerate((("LOFO", "#2c3e50"), ("LOIO", "#7f8c8d"))):
            b = [p.loc[p[key] == g, f"{mode}_panel_only_T1_before"].mean() for g in groups]
            a = [p.loc[p[key] == g, f"{mode}_panel_only_T1_after"].mean() for g in groups]
            off = -0.18 + 0.36 * j
            ax.scatter(x + off, b, marker="o", facecolor="white", edgecolor=col, s=18, zorder=3,
                       label=f"{mode} before" if ax is axes[0] else None)
            ax.scatter(x + off, a, marker="o", color=col, s=18, zorder=3, label=f"{mode} after" if ax is axes[0] else None)
            for xi, bi, ai in zip(x + off, b, a):
                ax.plot([xi, xi], [bi, ai], color=col, lw=0.8, zorder=2)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{g}\n(n={int((p[key] == g).sum())})" for g in groups], rotation=30, ha="right")
        ax.set_title(title)
    axes[0].set_ylabel("Tier-1 calls per healthy genome")
    axes[0].legend(frameon=False, fontsize=6, loc="upper right")
    panel_letter(axes[0], "a")
    panel_letter(axes[1], "b")
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.30)
    fig.savefig(out)
    plt.close(fig)


def main(sweep_dir, loci_path, null_path, outdir):
    chosen = json.load(open(f"{sweep_dir}/chosen_parameters.json"))
    sweep = pd.read_csv(f"{sweep_dir}/sweep_curve.tsv", sep="\t")
    per_ind = pd.read_csv(f"{sweep_dir}/per_individual_leaveout.tsv", sep="\t")
    loci = pd.read_csv(loci_path, sep="\t", usecols=["mode", "Panel_Status", "Cmplx_Class", "Panel_Matched_Fraction"])
    null = pd.read_csv(null_path, sep="\t", usecols=["null_status", "null_matched_fraction"])
    fig_curve(sweep, chosen, f"{outdir}/fig_tier1_vs_parameter.png")
    fig_fraction(loci, null, f"{outdir}/fig_matched_fraction.png")
    fig_before_after(per_ind, f"{outdir}/fig_before_after.png", 5.0)


if __name__ == "__main__":
    main(*sys.argv[1:5])
