#!/usr/bin/env python3
"""Compare the item-13 downstream rerun of a sample subset with panel-v1.0 (copies kept
in <sample>/v1.0/ by downstream_rerun.sh): CxSV loci, oversized span, junction table
(every column but locus_id must be identical: the panel junction sites do not change),
cxsv_only.vcf records, and the per-sample cost of the re-run rules (benchmarks/).
Usage: compare_rerun.py SAMPLES_ROOT SUBSET.tsv OUT.tsv"""
import csv, glob, os, sys

root, sub, out = sys.argv[1:4]


def summ(p):
    n, span, iv = 0, 0, []
    for r in csv.DictReader(open(p), delimiter="\t"):
        if r["IsCxSV"] == "YES":
            n += 1
        if r["CxSV_Criteria"] == "FILTERED_oversized_cluster":
            iv.append((r["Chr"], int(r["Start"]), int(r["End"])))
    cur = None
    for c, s, e in sorted(iv):
        if cur and cur[0] == c and s <= cur[2]:
            cur[2] = max(cur[2], e)
        else:
            if cur:
                span += cur[2] - cur[1]
            cur = [c, s, e]
    if cur:
        span += cur[2] - cur[1]
    return n, span, len(iv)


def junc(p):
    rows = []
    with open(p) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        k = h.index("locus_id")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rows.append(tuple(f[:k] + f[k + 1:]))
    return rows


def nrec(p):
    return sum(1 for l in open(p) if l[0] != "#")


cols = ["cohort", "sample", "cxsv_v10", "cxsv_new", "oversized_bp_v10", "oversized_bp_new", "n_oversized_v10",
        "n_oversized_new", "junction_rows", "junctions_identical_but_locus", "cxsv_only_records_v10",
        "cxsv_only_records_new", "downstream_cpu_s", "downstream_wall_s", "downstream_max_rss_mb"]
with open(out, "w") as fo:
    fo.write("\t".join(cols) + "\n")
    for r in csv.DictReader(open(sub), delimiter="\t"):
        d = f"{root}/{r['cohort']}/{r['sample_id']}"
        a, b = summ(f"{d}/v1.0/cxsv_summary.tsv"), summ(f"{d}/cxsv/cxsv_summary.tsv")
        ja, jb = junc(f"{d}/v1.0/junctions.tsv"), junc(f"{d}/junctions/junctions.tsv")
        cpu = wall = rss = 0.0
        for bm in glob.glob(f"{d}/benchmarks/sample_*.tsv"):
            x = list(csv.DictReader(open(bm), delimiter="\t"))[-1]
            num = lambda v: float(v) if v not in ("", "NA", "-") else 0.0
            wall += num(x["s"]); cpu += num(x["cpu_time"]); rss = max(rss, num(x["max_rss"]))
        fo.write("\t".join(map(str, [r["cohort"], r["sample_id"], a[0], b[0], a[1], b[1], a[2], b[2], len(jb),
                                     int(ja == jb), nrec(f"{d}/v1.0/cxsv_only.vcf"), nrec(f"{d}/cxsv/cxsv_only.vcf"),
                                     round(cpu, 1), round(wall, 1), round(rss, 1)])) + "\n")
