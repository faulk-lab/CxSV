#!/usr/bin/env python3
"""Callability of every panel junction site: how many panel samples could have shown it.

A site is callable in a sample when both of its breakpoints are callable in
that sample's CALLABLE mask (depth >= panel_callability.min_depth from reads
of MAPQ >= sv_params.min_mapq, on at least one `flank`-bp side of each
breakpoint; workflow/lib/cxsv_callability.py). An individual is callable when
any of its samples is.

Inputs
------
--sites              panel_junction_sites.tsv (shared schema; representative coordinates)
--callability-dir    directory with samples.tsv + <cohort>__<sample>.callable.bed.gz
--panel-individuals  panel_individuals.tsv (denominator, platforms)
--flank              bp (panel_callability.flank)

Outputs
-------
--out      panel_junction_callability.tsv, one row per site in the order of
           --sites: chrom1 pos1 chrom2 pos2 svtype svlen ac n_callable_samples
           n_callable n_individuals callable_fraction callable_<PLATFORM>...
           (n_callable = individuals; ac copied from the site for side-by-side reading)
--summary  panel_callability_summary.tsv: per-sample callable bp, and site-level
           totals (all / >= half / none of the individuals callable)
"""
import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import cxsv_callability as cc  # noqa: E402


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--sites", required=True)
    ap.add_argument("--callability-dir", required=True)
    ap.add_argument("--panel-individuals", required=True)
    ap.add_argument("--flank", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--summary", required=True)
    a = ap.parse_args()

    sites = read_tsv(a.sites)
    inds = {r["individual_id"]: r for r in read_tsv(a.panel_individuals)}
    samples = cc.read_samples(a.callability_dir)
    platforms = sorted({s["platform"] for s in samples})
    n_inds = len(inds)
    n_samples, callable_inds, by_platform = cc.callable_counts(a.callability_dir, sites, a.flank, inds)

    cols = ["chrom1", "pos1", "chrom2", "pos2", "svtype", "svlen", "ac", "n_callable_samples", "n_callable",
            "n_individuals", "callable_fraction"] + [f"callable_{p}" for p in platforms]
    n_all = n_half = n_none = 0
    with open(a.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for k, s in enumerate(sites):
            n = len(callable_inds[k])
            n_all += n == n_inds
            n_half += 2 * n >= n_inds
            n_none += n == 0
            row = [s["chrom1"], s["pos1"], s["chrom2"], s["pos2"], s["svtype"], s["svlen"], s["ac"],
                   int(n_samples[k]), n, n_inds, f"{n / n_inds:.6f}" if n_inds else "NA"]
            row += [len(by_platform[k].get(p, ())) for p in platforms]
            fh.write("\t".join(map(str, row)) + "\n")

    with open(a.summary, "w") as fh:
        fh.write("level\tname\tvalue\n")
        fh.write(f"definition\tflank_bp\t{a.flank}\n")
        fh.write(f"sites\tn_sites\t{len(sites)}\n")
        fh.write(f"sites\tcallable_in_all_individuals\t{n_all}\n")
        fh.write(f"sites\tcallable_in_at_least_half\t{n_half}\n")
        fh.write(f"sites\tcallable_in_none\t{n_none}\n")
        for s in samples:
            fh.write(f"sample_callable_bp\t{s['cohort']}/{s['sample_id']}\t{s['callable_bp']}\n")
    print(f"site callability: {len(sites)} sites x {len(samples)} samples / {n_inds} individuals "
          f"(flank {a.flank} bp): callable in all {n_all}, in >= half {n_half}, in none {n_none}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
