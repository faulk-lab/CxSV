#!/usr/bin/env python3
"""Per-genome discovery summary of downstream_rerun.sh arms (one per long_span_rule) on
the item-13 evaluation subset, and, optionally, a check that an arm reproduces an
earlier run of the same rule (cxsv_summary.tsv and clusters.bed identical).

  compare_rules.py SUBSET.tsv OUT.tsv ARM=SAMPLES_ROOT ... [--same ARM=OTHER_SAMPLES_ROOT ...]

Per genome and arm: CxSV loci (IsCxSV YES), oversized clusters and their summed Size
(the blind span; for clusters joined by paired long records Size is the sum of their
parts, classify_cxsv.py), long records by the rule applied (long_span_records.tsv), and
clusters joined by paired records (clusters.bed rows with an 8th column).
"""
import csv
import hashlib
import os
import sys


def md5(p):
    return hashlib.md5(open(p, "rb").read()).hexdigest() if os.path.exists(p) else "missing"


def main():
    args = sys.argv[1:]
    sub, out = args[0], args[1]
    arms, same = [], []
    rest = args[2:]
    target = arms
    for a in rest:
        if a == "--same":
            target = same
            continue
        target.append(a.split("=", 1))
    samples = [(r["cohort"], r["sample_id"], r["group"], r["coverage_x"])
               for r in csv.DictReader(open(sub), delimiter="\t")]
    with open(out, "w") as fh:
        fh.write("arm\tcohort\tsample\tgroup\tcoverage_x\tcxsv_loci\toversized_clusters\tblind_span_bp\t"
                 "long_records\tlong_paired\tlong_excluded\tlong_breakends\tjoined_clusters\n")
        for arm, root in arms:
            for c, s, g, cov in samples:
                d = f"{root}/{c}/{s}/cxsv"
                rows = list(csv.DictReader(open(f"{d}/cxsv_summary.tsv"), delimiter="\t"))
                loci = sum(r["IsCxSV"] == "YES" for r in rows)
                over = [r for r in rows if r["CxSV_Criteria"] == "FILTERED_oversized_cluster"]
                rule = {}
                if os.path.exists(f"{d}/long_span_records.tsv"):
                    for r in csv.DictReader(open(f"{d}/long_span_records.tsv"), delimiter="\t"):
                        rule[r["RULE"]] = rule.get(r["RULE"], 0) + 1
                joined = set()
                with open(f"{d}/clusters.bed") as cf:
                    for line in cf:
                        f = line.rstrip("\n").split("\t")
                        if len(f) > 7:
                            joined.add(f[6])
                fh.write("\t".join(map(str, [arm, c, s, g, cov, loci, len(over), sum(int(r["Size"]) for r in over),
                                             sum(rule.values()), rule.get("paired", 0), rule.get("exclude", 0),
                                             rule.get("breakends", 0), len(joined)])) + "\n")
    roots = dict(arms)
    for arm, other in same:
        diff = [(c, s, f) for c, s, _, _ in samples for f in ("cxsv_summary.tsv", "clusters.bed")
                if md5(f"{roots[arm]}/{c}/{s}/cxsv/{f}") != md5(f"{other}/{c}/{s}/cxsv/{f}")]
        print(f"{arm} vs {other}: {len(samples) - len({(c, s) for c, s, _ in diff})}/{len(samples)} genomes identical"
              + (f"; differ: {diff[:6]}" if diff else ""))


if __name__ == "__main__":
    main()
