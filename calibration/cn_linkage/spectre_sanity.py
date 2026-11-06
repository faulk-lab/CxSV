#!/usr/bin/env python3
"""Spectre with and without a PoN blacklist against the copy-number track's segments
(roadmap item 15 Phase 2, step 4; a sanity check, not a benchmark).

    spectre_sanity.py --segments segments_all.tsv --truth L01.bed --out summary.tsv \
        --calls KEY=RUN=cnv.bed [...] --blacklist KEY=blacklist.bed [...]

cnv.bed: chrom start end type copies (spectre_compare.sh). segments_all.tsv: the Phase 1
segments (key, chrom, start, end, state, ploidy, ...). Per key and run: calls by type,
called Mb, the fraction of called bp inside the blacklist, calls matching one of our
segments (same direction, 50 % reciprocal overlap) and, for keys with --truth, the
fraction of the truth intervals covered by Spectre DUP calls and by our gain segments.
"""
import argparse
import csv
from collections import defaultdict


def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def covered(a, ivs):
    """bp of interval a=(s,e) covered by merged ivs."""
    return sum(max(0, min(a[1], e) - max(a[0], s)) for s, e in ivs)


def recip(a, b, f=0.5):
    ov = min(a[1], b[1]) - max(a[0], b[0])
    return ov > 0 and ov >= f * (a[1] - a[0]) and ov >= f * (b[1] - b[0])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--segments", required=True)
    ap.add_argument("--calls", action="append", required=True, help="KEY=RUN=cnv.bed")
    ap.add_argument("--blacklist", action="append", default=[], help="KEY=blacklist.bed")
    ap.add_argument("--truth", action="append", default=[], help="KEY=truth.bed (chrom start end)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    ours = defaultdict(list)
    with open(a.segments, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            d = "gain" if int(r["state"]) > int(r["ploidy"]) else "loss"
            ours[r["key"]].append((r["chrom"], int(r["start"]), int(r["end"]), d))
    bl = {}
    for spec in a.blacklist:
        k, p = spec.split("=", 1)
        iv = defaultdict(list)
        for l in open(p):
            f = l.split()
            iv[f[0]].append((int(f[1]), int(f[2])))
        bl[k] = {c: merge(v) for c, v in iv.items()}
    truth = {}
    for spec in a.truth:
        k, p = spec.split("=", 1)
        truth[k] = [(f[0], int(f[1]), int(f[2])) for f in (l.split() for l in open(p)) if f]
    cols = ["key", "run", "n_calls", "n_DEL", "n_DUP", "called_Mb", "frac_called_bp_in_blacklist",
            "n_our_segments", "calls_matching_ours", "ours_matched", "truth_bp_covered_spectre_DUP",
            "truth_bp_covered_ours_gain"]
    rows = []
    for spec in a.calls:
        k, run, p = spec.split("=", 2)
        calls = []
        for l in open(p):
            f = l.rstrip("\n").split("\t")
            if len(f) >= 4 and f[2]:
                calls.append((f[0], int(f[1]), int(f[2]), f[3]))
        tot = sum(e - s for _, s, e, _ in calls)
        inbl = sum(covered((s, e), bl.get(k, {}).get(c, [])) for c, s, e, _ in calls)
        o = ours.get(k, [])
        def same(c, x):
            return c[0] == x[0] and ((c[3] == "DUP") == (x[3] == "gain")) and recip(c[1:3], x[1:3])
        m_calls = sum(any(same(c, x) for x in o) for c in calls)
        m_ours = sum(any(same(c, x) for c in calls) for x in o)
        t_sp = t_us = "NA"
        if k in truth:
            tb = sum(e - s for _, s, e in truth[k])
            dup = defaultdict(list); gain = defaultdict(list)
            for c, s, e, t in calls:
                if t == "DUP":
                    dup[c].append((s, e))
            for c, s, e, d in o:
                if d == "gain":
                    gain[c].append((s, e))
            t_sp = round(sum(covered((s, e), merge(dup[c])) for c, s, e in truth[k]) / tb, 3)
            t_us = round(sum(covered((s, e), merge(gain[c])) for c, s, e in truth[k]) / tb, 3)
        rows.append([k, run, len(calls), sum(c[3] == "DEL" for c in calls), sum(c[3] == "DUP" for c in calls),
                     round(tot / 1e6, 2), round(inbl / tot, 3) if tot else "NA", len(o), m_calls, m_ours, t_sp, t_us])
    with open(a.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    for r in rows:
        print("\t".join(map(str, r)))


if __name__ == "__main__":
    main()
