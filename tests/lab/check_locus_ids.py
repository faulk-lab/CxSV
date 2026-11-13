#!/usr/bin/env python3
"""Check that every locus id in CxSV tables is well formed (bundle validation check 10; ROADMAP item 18).

    check_locus_ids.py --fai hs1.fa.fai TABLE.tsv [TABLE.tsv ...]

A locus id is well formed when it is CONTIG:START-END with CONTIG a sequence of the reference
.fai, 0 <= START < length(CONTIG) and START < END <= length(CONTIG) + 1. The id column is
Locus_ID or locus_id (comment lines '#' skipped). Prints one line per malformed id and a
summary; exit 1 if any.

END = length + 1 is accepted (PI decision 22, 2026-11-13). vcf_to_bed.py writes a BND breakend as
[POS - bnd_window, POS + bnd_window) without clamping to the contig, so a breakend on the last
base of a contig gives END = length + 1 (24 panel loci in panel-v1.2-r2 and in panel-v1.3; fix in
calling code for v1.4, ROADMAP item 19). Anything further past the end, a start on or past the end
of the contig, zero length and contigs absent from the reference (chrN) are still malformed.
"""
import argparse
import csv
import re
import sys

ID = re.compile(r"^(\S+):(\d+)-(\d+)$")
END_SLACK = 1   # END = length + 1 from the unclamped BND window (ROADMAP item 19)


def problem(lid, lengths):
    m = ID.match(lid or "")
    if not m:
        return "not CONTIG:START-END"
    c, s, e = m.group(1), int(m.group(2)), int(m.group(3))
    if c not in lengths:
        return f"contig {c} not in the reference"
    if not s < e:
        return "start >= end"
    if s >= lengths[c]:
        return f"start {s} off the end of {c} ({lengths[c]})"
    if e > lengths[c] + END_SLACK:
        return f"end {e} past the end of {c} ({lengths[c]}; allowed {lengths[c] + END_SLACK})"
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
