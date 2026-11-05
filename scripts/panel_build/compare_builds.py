#!/usr/bin/env python3
"""Per-genome comparison of two panel builds' per-sample CxSV outputs
(panel-v1.0 vs the downstream-only rebuild panel-v1.1).

    compare_builds.py OLD/results/panel NEW/results/panel MANIFEST OUT.tsv

Per genome: CxSV loci (cxsv_master_table.tsv rows), oversized clusters
(cxsv_summary.tsv, CxSV_Criteria FILTERED_oversized_cluster) and their summed
span (the blind span), records listed in long_span_records.tsv (new build),
junction count, and whether junctions.tsv is identical apart from the locus_id
column (sha256 of the remaining columns).
"""
import csv
import hashlib
import os
import sys

OVERSIZED = "FILTERED_oversized_cluster"


def summary(root, c, s):
    d = f"{root}/samples/{c}/{s}"
    with open(f"{d}/cxsv/cxsv_master_table.tsv") as fh:
        loci = sum(1 for _ in fh) - 1
    n_over = span = 0
    with open(f"{d}/cxsv/cxsv_summary.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if OVERSIZED in r["CxSV_Criteria"]:
                n_over += 1
                span += int(r["Size"])
    h = hashlib.sha256()
    nj = 0
    with open(f"{d}/junctions/junctions.tsv") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        keep = [i for i, x in enumerate(hdr) if x != "locus_id"]
        h.update("\t".join(hdr[i] for i in keep).encode())
        for line in fh:
            f = line.rstrip("\n").split("\t")
            h.update(("\n" + "\t".join(f[i] for i in keep)).encode())
            nj += 1
    lsr = f"{d}/cxsv/long_span_records.tsv"
    n_long = (sum(1 for _ in open(lsr)) - 1) if os.path.exists(lsr) else None
    return loci, n_over, span, nj, h.hexdigest()[:16], n_long


def main():
    old, new, man, out = sys.argv[1:5]
    rows = list(csv.DictReader(open(man), delimiter="\t"))
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["cohort", "sample_id", "coverage_x", "loci_v1.0", "loci_v1.1", "oversized_v1.0", "oversized_v1.1",
                    "blind_span_v1.0", "blind_span_v1.1", "junctions_v1.0", "junctions_v1.1",
                    "junctions_identical_ex_locus_id", "long_span_records_v1.1"])
        for i, r in enumerate(rows):
            a = summary(old, r["cohort"], r["sample_id"])
            b = summary(new, r["cohort"], r["sample_id"])
            w.writerow([r["cohort"], r["sample_id"], r["coverage_x"], a[0], b[0], a[1], b[1], a[2], b[2], a[3], b[3],
                        "yes" if a[4] == b[4] else "NO", b[5]])
            if i % 200 == 0:
                print(i, file=sys.stderr)


if __name__ == "__main__":
    main()
