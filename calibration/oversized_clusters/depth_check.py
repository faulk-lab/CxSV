#!/usr/bin/env python3
"""Independent read-depth check of giant (> max_cluster_size) DEL/DUP calls.

For each test (depth_tests.tsv: sample, a same-cohort control genome without a giant
call over the interval, 4 interior and 2-4 flanking positions chosen outside segdups
and satellites), count MAPQ>=20 primary reads in 20 kb windows with samtools in both
genomes. Copy-number ratio = (sample interior / sample flank) / (control interior /
control flank); a heterozygous deletion gives ~0.5, a heterozygous duplication ~1.5,
an unchanged interval ~1.0.

Usage: depth_check.py TESTS.tsv DATASETS_ROOT REF OUT.tsv [THREADS]
"""
import csv, subprocess, sys
from concurrent.futures import ThreadPoolExecutor

tests, root, ref, out = sys.argv[1:5]
threads = int(sys.argv[5]) if len(sys.argv) > 5 else 8
W = 10000


def count(aln, chrom, p):
    r = subprocess.run(["samtools", "view", "-c", "-q", "20", "-F", "0xF04", "-T", ref,
                        f"{root}/{aln}", f"{chrom}:{max(1, p - W)}-{p + W}"],
                       capture_output=True, text=True)
    return int(r.stdout.strip()) if r.returncode == 0 and r.stdout.strip() else -1


def one(t):
    inner = list(map(int, t["inner"].split(",")))
    flank = list(map(int, t["flank"].split(",")))
    res = {}
    for who in ("s", "c"):
        a = t[f"aln_{who}"]
        res[f"{who}_inner"] = [count(a, t["chrom"], p) for p in inner]
        res[f"{who}_flank"] = [count(a, t["chrom"], p) for p in flank]
    return t, res


rows = list(csv.DictReader(open(tests), delimiter="\t"))
with open(out, "w") as fh, ThreadPoolExecutor(threads) as ex:
    fh.write("test_id\ts_inner\ts_flank\tc_inner\tc_flank\tcn_ratio\n")
    for t, r in ex.map(one, rows):
        si, sf, ci, cf = (sum(r[k]) for k in ("s_inner", "s_flank", "c_inner", "c_flank"))
        ok = min(r["s_inner"] + r["s_flank"] + r["c_inner"] + r["c_flank"]) >= 0
        ni, nf = len(r["s_inner"]), len(r["s_flank"])
        ratio = (si / ni / (sf / nf)) / (ci / ni / (cf / nf)) if ok and sf and ci and cf else float("nan")
        fh.write("\t".join([t["test_id"]] + [",".join(map(str, r[k])) for k in ("s_inner", "s_flank", "c_inner", "c_flank")]
                           + [f"{ratio:.3f}"]) + "\n")
        fh.flush()
