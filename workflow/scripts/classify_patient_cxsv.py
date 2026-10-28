#!/usr/bin/env python3
"""Label each patient CxSV locus public or private against the panel.

For every locus in the patient's cxsv_master_table.tsv, counts the distinct
panel individuals with a CxSV call that overlaps it after padding both sides by
--window bp (the same tolerance as the previous `bedtools window` filter). A
locus with at least --min-carriers carriers is public; otherwise private.

The carrier count and frequency are kept on every row, so a rare-but-seen
locus is distinguishable from one absent from the panel. Breakpoint-level
matching and callability weighting replace interval overlap in roadmap item 3.

Outputs
-------
--out-all      master table + Panel_Carriers, Panel_Carrier_Freq,
               Panel_Carriers_By_Platform, Panel_Status
--out-private  the private rows only (same columns)
--out-bed      private loci as BED (Chr, Start, End) for VCF extraction
"""
import argparse
import bisect
import csv
import sys
from collections import defaultdict


def interval(row):
    start, end = max(0, int(row["Start"])), int(row["End"])
    if end < start:
        start, end = end, start
    if end <= start:
        end = start + 1
    return start, end


def load_panel(calls_tsv):
    by_chr = defaultdict(list)
    with open(calls_tsv, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            by_chr[r["Chr"]].append((int(r["Start"]), int(r["End"]), r["Individual_ID"], r["Platform"]))
    index = {}
    for chrom, rows in by_chr.items():
        rows.sort()
        index[chrom] = ([r[0] for r in rows], rows, max(r[1] - r[0] for r in rows))
    return index


def panel_individuals(summary_tsv):
    with open(summary_tsv, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["level"] == "total":
                return int(r["n_individuals"])
    raise ValueError(f"{summary_tsv}: no total row")


def carriers(index, chrom, start, end, window):
    if chrom not in index:
        return set(), {}
    starts, rows, maxlen = index[chrom]
    lo, hi = start - window, end + window
    i = bisect.bisect_left(starts, lo - maxlen)
    found, by_plat = set(), defaultdict(set)
    while i < len(rows) and rows[i][0] <= hi:
        s, e, ind, plat = rows[i]
        if e >= lo:
            found.add(ind)
            by_plat[plat].add(ind)
        i += 1
    return found, by_plat


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--master-table", required=True)
    ap.add_argument("--panel-calls", required=True)
    ap.add_argument("--panel-summary", required=True)
    ap.add_argument("--window", type=int, default=1000)
    ap.add_argument("--min-carriers", type=int, default=1)
    ap.add_argument("--out-all", required=True)
    ap.add_argument("--out-private", required=True)
    ap.add_argument("--out-bed", required=True)
    a = ap.parse_args()

    index = load_panel(a.panel_calls)
    n_panel = panel_individuals(a.panel_summary)
    with open(a.master_table, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols = reader.fieldnames + ["Panel_Carriers", "Panel_Carrier_Freq",
                                    "Panel_Carriers_By_Platform", "Panel_Status"]
        rows = list(reader)

    n_priv = 0
    with open(a.out_all, "w", newline="") as fa, open(a.out_private, "w", newline="") as fp, \
            open(a.out_bed, "w") as fb:
        wa = csv.DictWriter(fa, cols, delimiter="\t", lineterminator="\n")
        wp = csv.DictWriter(fp, cols, delimiter="\t", lineterminator="\n")
        wa.writeheader()
        wp.writeheader()
        for r in rows:
            start, end = interval(r)
            found, by_plat = carriers(index, r["Chr"], start, end, a.window)
            r["Panel_Carriers"] = len(found)
            r["Panel_Carrier_Freq"] = f"{len(found) / n_panel:.6f}" if n_panel else "NA"
            r["Panel_Carriers_By_Platform"] = ";".join(f"{p}:{len(v)}" for p, v in sorted(by_plat.items())) or "none"
            r["Panel_Status"] = "public" if len(found) >= a.min_carriers else "private"
            wa.writerow(r)
            if r["Panel_Status"] == "private":
                wp.writerow(r)
                fb.write(f"{r['Chr']}\t{start}\t{end}\n")
                n_priv += 1
    print(f"{a.master_table}: {len(rows)} CxSV loci, {n_priv} private, {len(rows) - n_priv} public "
          f"(panel of {n_panel} individuals, window {a.window} bp, min carriers {a.min_carriers})",
          file=sys.stderr)


if __name__ == "__main__":
    main()
