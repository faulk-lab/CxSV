#!/usr/bin/env python3
"""Segments and flanks of every patient CxSV locus, for read-depth copy number (criterion 5).

For each locus (cxsv_panel_status.tsv) the distinct junctions
(junctions/junctions.tsv, duplicates collapsed and SVIM composites dropped as in
workflow/lib/cxsv_classes.py) give the locus's breakends on each chromosome.
Between consecutive breakends lies a *segment* (trimmed by --margin bp at both
ends, so imprecise breakpoints do not bias the depth); outside the outermost
breakends a left and a right *flank* of --flank bp (after a --gap) give the
local baseline. Output BED (0-based, half-open), sorted:

    chrom start end locus_id kind index      kind: flank_left | segment | flank_right

The rule then runs `samtools bedcov` on the patient alignment over these
intervals; scripts/score_cxsv.py turns mean depth into copy number.
"""
import argparse
import csv
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import cxsv_classes as cc     # noqa: E402
import cxsv_junctions as cj   # noqa: E402


def rows_to_junctions(rows):
    return [cj.Junction(r["chrom1"], r["pos1"], r["chrom2"], r["pos2"], r["svtype"], r["svlen"],
                        r["supp_vec"], r["dosage"], r["record_ids"].split(",")) for r in rows]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--junctions", required=True)
    ap.add_argument("--loci", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--collapse-distance", type=int, default=300)
    ap.add_argument("--composite-fraction", type=float, default=0.25)
    ap.add_argument("--flank", type=int, default=20000)
    ap.add_argument("--gap", type=int, default=1000)
    ap.add_argument("--margin", type=int, default=200)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    order = cj.contig_order(a.fai)
    lengths = {}
    with open(a.fai) as fh:
        for line in fh:
            f = line.split("\t")
            lengths[f[0]] = int(f[1])
    by_locus = defaultdict(list)
    for r in cj.read_table(a.junctions):
        for L in r["loci"]:
            by_locus[L].append(r)
    with open(a.loci, newline="") as fh:
        loci = [f"{r['Chr']}:{r['Start']}-{r['End']}" for r in csv.DictReader(fh, delimiter="\t")]
    p = {**cc.DEFAULTS, "collapse_distance": a.collapse_distance, "composite_fraction": a.composite_fraction}
    out = []
    for L in loci:
        kept, _ = cc.distinct_junctions(rows_to_junctions(by_locus.get(L, [])), order, p)
        per_chrom = defaultdict(list)
        for chrom, pos in cc.breakends(kept, a.collapse_distance):
            per_chrom[chrom].append(pos)
        for chrom, ps in per_chrom.items():
            ps.sort()
            n = lengths.get(chrom, max(ps) + a.flank + a.gap + 1)
            lo, hi = max(0, ps[0] - a.gap - a.flank), max(0, ps[0] - a.gap)
            if hi > lo:
                out.append((chrom, lo, hi, L, "flank_left", 0))
            for i, (s, e) in enumerate(zip(ps, ps[1:]), 1):
                s2, e2 = s + a.margin, e - a.margin
                if e2 > s2:
                    out.append((chrom, s2, e2, L, "segment", i))
            lo, hi = min(n, ps[-1] + a.gap), min(n, ps[-1] + a.gap + a.flank)
            if hi > lo:
                out.append((chrom, lo, hi, L, "flank_right", len(ps)))
    out.sort(key=lambda r: (order.get(r[0], 1 << 30), r[1], r[2], r[3], r[4]))
    with open(a.out, "w") as fh:
        for r in out:
            fh.write("\t".join(map(str, r)) + "\n")
    print(f"{len(loci)} loci -> {len(out)} depth intervals", file=sys.stderr)


if __name__ == "__main__":
    main()
