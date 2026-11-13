#!/usr/bin/env python3
"""Check that every locus id in CxSV tables is well formed (bundle validation check 10; ROADMAP item 18).

    check_locus_ids.py --fai hs1.fa.fai TABLE.tsv [TABLE.tsv ...]

A locus id is well formed when it is CONTIG:START-END with CONTIG a sequence of the reference
.fai and 0 <= START < END <= length(CONTIG). The id column is Locus_ID or locus_id (comment
lines '#' skipped). Prints one line per malformed id and a summary; exit 1 if any.
"""
import argparse
import csv
import re
import sys

ID = re.compile(r"^(\S+):(\d+)-(\d+)$")


def problem(lid, lengths):
    m = ID.match(lid or "")
    if not m:
        return "not CONTIG:START-END"
    c, s, e = m.group(1), int(m.group(2)), int(m.group(3))
    if c not in lengths:
        return f"contig {c} not in the reference"
    if not s < e:
        return "start >= end"
    if e > lengths[c]:
        return f"end {e} past the end of {c} ({lengths[c]})"
    return None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--fai", required=True)
    ap.add_argument("tables", nargs="+")
    a = ap.parse_args(argv)
    lengths = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(a.fai) if l.strip()}
    n = bad = 0
    for t in a.tables:
        with open(t, newline="") as fh:
            rows = csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t")
            col = next((c for c in ("Locus_ID", "locus_id") if c in (rows.fieldnames or [])), None)
            if col is None:
                continue
            seen = set()
            for r in rows:
                lid = r[col]
                if lid in seen:
                    continue
                seen.add(lid)
                n += 1
                p = problem(lid, lengths)
                if p:
                    bad += 1
                    print(f"MALFORMED\t{t}\t{lid}\t{p}")
    print(f"{n} locus ids in {len(a.tables)} table(s): {bad} malformed")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
