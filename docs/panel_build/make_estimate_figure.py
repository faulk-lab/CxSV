#!/usr/bin/env python3
"""Figure for docs/PANEL_BUILD.md: (a) per-rule CPU-hours of the measured test
genomes, (b) projected hours of the per-sample build step per dispatch option.

    python3 docs/panel_build/make_estimate_figure.py [--site SITE] [--registry REG] [--out PNG]

Defaults read the lab registry on the NAS; the estimates come from
scripts/panel_build/make_shards.py with the measured cost model in
docs/panel_build/measurements_20261029/. Needs matplotlib and pandas.
"""
import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "scripts", "panel_build"))
import make_shards as ms  # noqa: E402
import shardlib  # noqa: E402

MEAS = os.path.join(REPO, "docs", "panel_build", "measurements_20261029")
GENOMES = [("NA19734", "ONT 1kGP-1019 (CRAM)"), ("HG01435", "ONT LRSC R10"), ("HG00240", "ONT LRSC R9"),
           ("NA20819", "PacBio HiFi, 49 GB BAM"), ("HG01795", "PacBio HiFi, 373 GB BAM")]
RULES = [("run_svim", "SVIM (1 thread)", "#4C72B0"), ("run_cutesv", "cuteSV", "#DD8452"),
         ("run_sniffles", "Sniffles2", "#55A868"), ("callable_mask", "mosdepth callability", "#8172B3")]
FREYA_SPEEDS = {"freya": 0.56, "shadowfax": 0.69, "athena": 0.862, "betsy": 0.872, "keurig": 0.658}


def scenario(site, fleet_path, reg, cm, exclude=(), sets=(), nas=None, cohorts=None):
    fleet = shardlib.load_yaml(fleet_path)
    ms.apply_sets(fleet, list(sets))
    if nas:
        fleet["nas_aggregate_mb_s"] = nas
    return ms.plan(site, fleet, reg, cm, exclude, cohorts=cohorts)[4]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--site", default=os.path.join(REPO, "config/sites/faulk_lab_panel.yaml"))
    ap.add_argument("--fleet", default=os.path.join(REPO, "config/sites/faulk_lab_fleet.yaml"))
    ap.add_argument("--registry", default="/mnt/sargasso1/CxSV/datasets/registry/sample_registry.tsv")
    ap.add_argument("--out", default=os.path.join(REPO, "docs/figures/panel_build_cost_estimate.png"))
    a = ap.parse_args()
    site = shardlib.load_yaml(a.site)
    reg = [r for r in shardlib.read_sheet(a.registry)[1] if r.get("status") == "verified"]
    cm = ms.read_cost_model(os.path.join(MEAS, "cost_model.tsv"))
    pr = pd.read_csv(os.path.join(MEAS, "per_rule_cost.tsv"), sep="\t")

    A = scenario(site, a.fleet, reg, cm)
    opts = [("All 5 nodes now (recommended)", A["makespan_h"]),
            ("freya capped at 16 threads", scenario(site, a.fleet, reg, cm, sets=["freya.reserve_threads=16"])["makespan_h"]),
            ("freya nights only (18–08 h)", scenario(site, a.fleet, reg, cm, sets=[
                f"freya.speed={0.560 * 24 / 14:.3f}", f"freya.nic_mb_s={496 * 14 / 24:.0f}"])["makespan_h"]),
            ("without freya", scenario(site, a.fleet, reg, cm, exclude=["freya"])["makespan_h"]),
            ("NAS throttled to 70%", scenario(site, a.fleet, reg, cm, nas=440)["makespan_h"]),
            ("staged: ONT-1019, then LRSC + PacBio",
             scenario(site, a.fleet, reg, cm, cohorts=["1kGP_ONT_1019_CHM13"])["makespan_h"]
             + scenario(site, a.fleet, reg, cm, cohorts=["1kGP_LRSC_ONT_CHM13", "1kGP_PacBio_CHM13"])["makespan_h"]),
            ("all 5 nodes, 1.5x slower per sample", scenario(site, a.fleet, reg, cm, sets=[
                f"{n}.speed={v * 1.5:.3f}" for n, v in FREYA_SPEEDS.items()])["makespan_h"])]

    plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.3), gridspec_kw={"width_ratios": [1.15, 1]})
    labels = []
    for i, (sid, lab) in enumerate(GENOMES):
        d = pr[pr.sample_id == sid].set_index("rule")
        left = 0
        for r, nm, col in RULES:
            v = d.loc[r, "cpu_s"] / 3600
            ax1.barh(i, v, left=left, color=col, height=0.62, label=nm if i == 0 else None)
            left += v
        other = d.drop(index=[r for r, _, _ in RULES] + ["ALL (instance)"], errors="ignore").cpu_s.sum() / 3600
        ax1.barh(i, other, left=left, color="0.75", height=0.62, label="other per-sample rules" if i == 0 else None)
        ax1.text(left + other + 0.03, i, f"{d.loc['ALL (instance)', 'wall_s'] / 60:.0f} min", va="center", fontsize=6)
        labels.append(f"{lab}\n{sid}, {d.coverage_x.iloc[0]:.1f}x")
    ax1.set_yticks(range(len(GENOMES)))
    ax1.set_yticklabels(labels, fontsize=6)
    ax1.invert_yaxis()
    ax1.set_xlabel("CPU-hours per genome")
    ax1.set_xlim(0, 2.35)
    ax1.set_title("A full genome costs 0.7–1.9 CPU-h", loc="left")
    ax1.legend(frameon=False, fontsize=6, loc="upper center", bbox_to_anchor=(0.45, -0.2), ncol=3,
               handlelength=1.2, columnspacing=0.8)
    for i, (nm, h) in enumerate(opts):
        col = "#C44E52" if i == 0 else "0.45"
        ax2.plot([0, h], [i, i], color=col, lw=1)
        ax2.plot(h, i, "o", color=col, ms=4)
        ax2.text(h + 1.2, i, f"{h:.0f} h", va="center", fontsize=6)
    ax2.axvline(A["nas_floor_h"], color="0.3", ls="--", lw=0.7)
    ax2.text(A["nas_floor_h"] + 1, -0.6, f"NAS read floor\n{A['nas_tb']:.1f} TB at 630 MB/s",
             fontsize=6, ha="left", va="center", color="0.3")
    ax2.set_yticks(np.arange(len(opts)))
    ax2.set_yticklabels([o[0] for o in opts], fontsize=6)
    ax2.set_ylim(len(opts) - 0.4, -1.1)
    ax2.set_xlabel(f"projected hours, per-sample step ({A['samples']:,} files)")
    ax2.set_xlim(0, 72)
    ax2.set_title("Every option sits above the NAS read floor", loc="left")
    for ax, l in ((ax1, "a"), (ax2, "b")):
        ax.text(-0.02, 1.04, l, transform=ax.transAxes, fontweight="bold", fontsize=10, ha="right")
    fig.tight_layout()
    fig.savefig(a.out, dpi=300)
    for nm, h in opts:
        print(f"{nm}\t{h:.1f}")


if __name__ == "__main__":
    main()
