#!/usr/bin/env python3
"""
filter_table_by_bed.py
========================
Keep or exclude rows of a TSV table that overlap regions in a BED file, with
an optional bp window padding each region before the overlap check. Chr/
Start/End columns are detected from the header (falls back to columns 0/1/2
if the expected names aren't found).

Used by workflow/rules/private_filter.smk to derive each query cohort's
private_cxsv_master_table.tsv: rows of cxsv_master_table.tsv matching a
locus in private_cxsv_population.bed (already computed with the configured
slop tolerance) are kept.

Usage
-----
python filter_table_by_bed.py \
    --table  results/patient_batchA/final/cxsv_master_table.tsv \
    --bed    results/patient_batchA/final/private_cxsv_population.bed \
    --mode   keep \
    --window 0 \
    --output results/patient_batchA/final/private_cxsv_master_table.tsv
"""

import argparse
import sys


def load_regions(bed_path: str, window: int = 0) -> list[tuple[str, int, int]]:
    """Parse a BED file into (chrom, start, end) tuples, padded by `window` bp."""
    regions = []
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            try:
                start = max(0, int(parts[1]) - window)
                end   = int(parts[2]) + window
            except ValueError:
                continue
            regions.append((parts[0], start, end))
    return regions


def overlaps(chrom: str, start: int, end: int,
             regions: list[tuple[str, int, int]]) -> bool:
    """Simple linear scan — fast enough for typical CxSV locus counts (<1000)."""
    for r_chrom, r_start, r_end in regions:
        if chrom == r_chrom and start <= r_end and end >= r_start:
            return True
    return False


def main():
    parser = argparse.ArgumentParser(
        description="Keep or exclude table rows overlapping a BED file."
    )
    parser.add_argument("--table",  required=True, help="Input TSV table (header required)")
    parser.add_argument("--bed",    required=True, help="BED file of regions to compare against")
    parser.add_argument("--mode",   required=True, choices=("keep", "exclude"),
                        help="keep: write rows overlapping the BED; "
                             "exclude: write rows NOT overlapping it")
    parser.add_argument("--window", type=int, default=0,
                        help="bp padding added to each side of every BED region (default 0)")
    parser.add_argument("--output", required=True, help="Output TSV table")
    args = parser.parse_args()

    regions = load_regions(args.bed, window=args.window)
    if not regions:
        sys.stderr.write(f"WARNING: No regions loaded from {args.bed}.\n")

    with open(args.table) as fh:
        header = fh.readline()
        cols = header.rstrip("\n").split("\t")
        try:
            chr_idx   = cols.index("Chr")
            start_idx = cols.index("Start")
            end_idx   = cols.index("End")
        except ValueError:
            chr_idx, start_idx, end_idx = 0, 1, 2

        n_written = 0
        n_total   = 0
        with open(args.output, "w") as out:
            out.write(header)
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= max(chr_idx, start_idx, end_idx):
                    continue
                n_total += 1
                try:
                    chrom = parts[chr_idx]
                    start = int(parts[start_idx])
                    end   = int(parts[end_idx])
                except ValueError:
                    continue

                hit = overlaps(chrom, start, end, regions)
                keep = hit if args.mode == "keep" else not hit

                if keep:
                    out.write(line)
                    n_written += 1

    print(f"[{args.mode}] wrote {n_written}/{n_total} rows to {args.output} "
          f"(window={args.window} bp, {len(regions)} regions from {args.bed}).")


if __name__ == "__main__":
    main()
