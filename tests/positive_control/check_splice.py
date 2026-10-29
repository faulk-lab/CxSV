#!/usr/bin/env python3
"""Depth-continuity check for the spliced positive control.

For every replace window (+/- 25 kb) compare host and spliced depth in 2 kb bins.
Writes splice_depth_bins.tsv and prints one summary line per event:
  inside  = spliced/host median depth ratio inside the window
  edge    = worst spliced/host ratio in the four bins straddling the two window edges
Copy-number changes are expected inside a window (deletions ~0.5x, duplications
~1.5x, triplications ~2x over the affected segment), so 'inside' is a whole-window
median; edges sit on unrearranged flank and should be ~1.0.

Usage: check_splice.py replace_all.bed host.cram spliced.cram hs1.fa out.tsv
"""
import statistics
import subprocess
import sys
from collections import defaultdict

PAD, BIN = 25000, 2000


def depth_bins(cram, ref, c, s, e):
    out = subprocess.run(f"samtools depth -a -Q 0 --reference {ref} -r {c}:{s + 1}-{e} {cram}",
                         shell=True, check=True, capture_output=True, text=True).stdout
    bins = defaultdict(list)
    for line in out.splitlines():
        _, p, d = line.split("\t")
        bins[(int(p) - 1 - s) // BIN].append(int(d))
    return {k: statistics.mean(v) for k, v in bins.items()}


def main(bed, host, spliced, ref, out):
    rows = [l.split("\t")[:4] for l in open(bed) if l.strip()]
    with open(out, "w") as o:
        o.write("event\tchrom\tbin_start\tbin_end\tposition\thost_depth\tspliced_depth\n")
        for c, s, e, eid in rows:
            s, e, eid = int(s), int(e), eid.strip()
            lo, hi = max(0, s - PAD), e + PAD
            h, p = depth_bins(host, ref, c, lo, hi), depth_bins(spliced, ref, c, lo, hi)
            inside_r, edge_r = [], []
            for k in sorted(h):
                b0, b1 = lo + k * BIN, lo + (k + 1) * BIN
                pos = "inside" if (b0 >= s and b1 <= e) else ("outside" if (b1 <= s or b0 >= e) else "edge")
                near_edge = abs(b0 - s) <= BIN or abs(b0 - e) <= BIN or abs(b1 - s) <= BIN or abs(b1 - e) <= BIN
                hd, pd = h[k], p.get(k, 0.0)
                o.write(f"{eid}\t{c}\t{b0}\t{b1}\t{pos}\t{hd:.2f}\t{pd:.2f}\n")
                if hd > 0:
                    if pos == "inside":
                        inside_r.append(pd / hd)
                    if near_edge and pos != "inside":
                        edge_r.append(pd / hd)
            ins = statistics.median(inside_r) if inside_r else float("nan")
            worst = max(edge_r, key=lambda r: abs(r - 1)) if edge_r else float("nan")
            print(f"{eid}\t{c}:{s}-{e}\tinside={ins:.2f}\tedge_worst={worst:.2f}", flush=True)


if __name__ == "__main__":
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    main(*sys.argv[1:])
