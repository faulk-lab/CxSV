#!/usr/bin/env python3
"""
link_paired_clusters.py
=======================
Join the two end clusters of every paired long record into one cluster
(cxsv_params.long_span_rule paired; roadmap item 13 follow-up).

vcf_to_bed.py writes a paired record (a DEL/DUP/INV longer than span_link_max
whose two ends each link to other records) as two breakend rows, ID and ID_end.
`bedtools cluster -d cluster_distance` puts them in the clusters at the two
ends; the record is an edge between those clusters, so they are one CxSV
candidate. Clusters connected by such edges (transitively) get one cluster id
(the first one in file order), and every row of a joined cluster gets an 8th
column with its original bedtools cluster id, so classify_cxsv.py can size the
cluster by its parts (the sum of the parts' spans) and not by the distance
between them: the long record's interval is never part of the size.

Rows of clusters that are not joined are written unchanged (7 columns), so with
long_span_rule exclude or breakends (no record with RULE paired) the output is
byte-identical to the input.

Usage
-----
python link_paired_clusters.py --clusters clusters.raw.bed --long long_span_records.tsv \\
    --out clusters.bed
"""
import argparse
import csv
import sys


def read_paired(path):
    """IDs of long records with RULE paired in long_span_records.tsv."""
    with open(path, newline="") as fh:
        return {r["ID"] for r in csv.DictReader(fh, delimiter="\t") if r.get("RULE") == "paired"}


def link(lines, paired):
    """Rows (list of str) with joined cluster ids; see the module docstring."""
    rows = [l.rstrip("\n").split("\t") for l in lines if l.strip()]
    ends = {}                                     # paired id -> [cluster ids of its rows]
    for f in rows:
        vid = f[4]
        rid = vid[:-4] if vid.endswith("_end") else vid
        if rid in paired:
            ends.setdefault(rid, []).append(f[6])
    parent, first = {}, {}
    for i, f in enumerate(rows):
        first.setdefault(f[6], i)

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for cids in ends.values():
        for c in cids[1:]:
            a, b = find(cids[0]), find(c)
            if a != b:                            # keep the cluster seen first as the root
                if first[a] <= first[b]:
                    parent[b] = a
                else:
                    parent[a] = b
    size = {}
    for c in first:
        size[find(c)] = size.get(find(c), 0) + 1
    out = []
    for f in rows:
        root = find(f[6])
        if size[root] > 1:
            out.append("\t".join(f[:6] + [root, f[6]]))
        else:
            out.append("\t".join(f))
    return out, sum(1 for r, n in size.items() if n > 1), sum(n for n in size.values() if n > 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--clusters", required=True, help="bedtools cluster output (7 columns)")
    ap.add_argument("--long", required=True, help="long_span_records.tsv (vcf_to_bed.py --long-out)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    paired = read_paired(a.long)
    with open(a.clusters) as fh:
        out, n_joined, n_parts = link(fh.readlines(), paired)
    with open(a.out, "w") as fh:
        fh.write("".join(l + "\n" for l in out))
    print(f"link_paired_clusters: {len(paired)} paired records; {n_joined} joined clusters "
          f"from {n_parts} parts", file=sys.stderr)


if __name__ == "__main__":
    main()
