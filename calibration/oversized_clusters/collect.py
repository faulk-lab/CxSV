#!/usr/bin/env python3
"""Collect characterise.py per-sample outputs into panel-wide tables.

Writes under --out:
  genome_variants.tsv        one row per sample x rule (blind span, clusters, CxSV loci, changes)
  oversized_clusters.tsv.gz  every oversized cluster of the current rule, with residual
                             oversized span and CxSV loci inside under each rule
  giant_records.tsv.gz       every record whose interval exceeds max_cluster_size
  new_loci_by_class.tsv      CxSV loci new (or lost) under each rule, by class and
                             whether a member record exceeds max_cluster_size
"""
import argparse
import glob
import gzip
import json
import os
from collections import Counter

ap = argparse.ArgumentParser()
ap.add_argument("--root", required=True)
ap.add_argument("--out", required=True)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)

dirs = sorted(glob.glob(f"{a.root}/*/*/done.ok"))
cols = None
n_ident = 0
with open(f"{a.out}/genome_variants.tsv", "w") as gv:
    for d in dirs:
        s = json.load(open(d.replace("done.ok", "summary.json")))
        n_ident += bool(s["v0_partition_identical"])
        for v, x in s["variants"].items():
            row = {"cohort": s["cohort"], "sample": s["sample"], "individual_id": s["individual_id"],
                   "platform": s["platform"], "coverage_x": s["coverage_x"],
                   "v0_partition_identical": int(s["v0_partition_identical"]), "variant": v, **x}
            if cols is None:
                cols = list(row)
                gv.write("\t".join(cols) + "\n")
            gv.write("\t".join(str(row[c]) for c in cols) + "\n")


def cat(name, out):
    hdr = None
    with gzip.open(out, "wt") as fo:
        for d in dirs:
            p = d.replace("done.ok", name)
            with open(p) as fh:
                h = fh.readline()
                if hdr is None and h.strip() != "cohort":
                    hdr = h
                    fo.write(h)
                for line in fh:
                    fo.write(line)


cat("oversized_clusters.tsv", f"{a.out}/oversized_clusters.tsv.gz")
cat("giant_records.tsv", f"{a.out}/giant_records.tsv.gz")
cnt = Counter()
for d in dirs:
    with open(d.replace("done.ok", "loci_by_variant.tsv")) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            v, cls, maxlen, incur = f[0], f[4], int(f[6]), f[8]
            if v == "V0_current":
                continue
            cnt[(v, cls, int(maxlen > 500000))] += 1
with open(f"{a.out}/new_loci_by_class.tsv", "w") as fh:
    fh.write("variant\tclass\thas_giant_member\tn_loci\n")
    for (v, cls, g), n in sorted(cnt.items()):
        fh.write(f"{v}\t{cls}\t{g}\t{n}\n")
print(len(dirs), "samples;", n_ident, "with the released partition reproduced exactly")
