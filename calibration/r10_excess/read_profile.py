#!/usr/bin/env python3
"""Read-level profile of panel alignments (R10-excess diagnosis).

For each genome: alignment header provenance (@RG, @PG ids) and, over one
window (default chr20:30-32 Mb), primary reads, fraction carrying a
supplementary alignment (SA tag), fraction with >= 2 SA segments, median read
length, fraction MAPQ < 20. Run in the pipeline image (pysam).

    read_profile.py --sheets DIR --out OUT.tsv cohort:sample ...
"""
import argparse
import csv
import statistics

import pysam

ap = argparse.ArgumentParser()
ap.add_argument("--sheets", required=True)
ap.add_argument("--region", default="chr20:30000000-32000000")
ap.add_argument("--ref", default="/mnt/sargasso1/CxSV/refs/hs1.fa")
ap.add_argument("--out", required=True)
ap.add_argument("genomes", nargs="+")
a = ap.parse_args()
chrom, rng = a.region.split(":")
s0, e0 = map(int, rng.split("-"))
with open(a.out, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["cohort", "sample_id", "reads", "frac_SA", "frac_SA_2plus", "median_len", "frac_mapq_lt20", "rg", "pg_ids"])
    for g in a.genomes:
        c, s = g.split(":")
        aln = next(r["alignment"] for r in csv.DictReader(open(f"{a.sheets}/{c}.tsv"), delimiter="\t") if r["sample_id"] == s)
        f = pysam.AlignmentFile(aln, reference_filename=a.ref)
        hd = f.header.to_dict()
        rg = ";".join(f"{k}={v}" for k, v in (hd.get("RG") or [{}])[0].items())[:200]
        pg = ",".join(sorted({p.get("ID", "") for p in hd.get("PG", [])}))[:200]
        n = sa = sa2 = lowq = 0
        lens = []
        for r in f.fetch(chrom, s0, e0):
            if r.is_secondary or r.is_supplementary or r.is_unmapped:
                continue
            n += 1
            lens.append(r.query_length or r.infer_read_length() or 0)
            lowq += r.mapping_quality < 20
            if r.has_tag("SA"):
                sa += 1
                sa2 += r.get_tag("SA").count(";") >= 2
        w.writerow([c, s, n, f"{sa / n:.4f}", f"{sa2 / n:.4f}", int(statistics.median(lens)), f"{lowq / n:.4f}", rg, pg])
        fh.flush()
