#!/usr/bin/env python3
"""
extract_cxsv_vcf.py
===================
Extract records from a merged VCF that overlap CxSV loci defined in
cxsv_summary.tsv (output of the classify_cxsv Snakemake rule).

Usage
-----
python extract_cxsv_vcf.py \
    --vcf results/final/long_read_filtered.vcf \
    --summary results/final/cxsv_summary.tsv \
    --output results/final/cxsv_only.vcf
"""

import argparse
import sys
from collections import defaultdict


def load_cxsv_regions(summary_path: str) -> list[tuple[str, int, int]]:
    """
    Parse cxsv_summary.tsv and return a list of (chrom, start, end)
    tuples for all rows where IsCxSV == YES.

    Column layout (0-based after split):
      0  ClusterID
      1  Chr
      2  Start
      3  End
      4  Size
      5  N_Breakpoints
      6  N_SVtypes
      7  SV_Types
      8  Has_Inter      ← NEW: YES if any BND/TRA mate is inter-chromosomal
      9  IsCxSV         ← shifted from col 8
    """
    regions = []
    with open(summary_path) as fh:
        header = fh.readline()
        # Detect column indices dynamically from header
        cols = header.rstrip("\n").split("\t")
        try:
            is_cxsv_idx = cols.index("IsCxSV")
            chr_idx     = cols.index("Chr")
            start_idx   = cols.index("Start")
            end_idx     = cols.index("End")
        except ValueError:
            # Fall back to positional defaults if header differs
            is_cxsv_idx, chr_idx, start_idx, end_idx = 9, 1, 2, 3

        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= is_cxsv_idx:
                continue
            if parts[is_cxsv_idx] == "YES":
                try:
                    regions.append((parts[chr_idx], int(parts[start_idx]), int(parts[end_idx])))
                except ValueError:
                    continue
    return regions


def overlaps(chrom: str, pos: int, end: int,
             regions: list[tuple[str, int, int]]) -> bool:
    """
    Return True if (chrom, pos, end) overlaps any region.
    Simple linear scan — fast enough for typical CxSV locus counts (<1000).
    For large region sets, replace with an interval tree.
    """
    for r_chrom, r_start, r_end in regions:
        if chrom == r_chrom and pos <= r_end and end >= r_start:
            return True
    return False


def main():
    parser = argparse.ArgumentParser(description="Extract CxSV records from a VCF.")
    parser.add_argument("--vcf",     required=True, help="Input VCF file")
    parser.add_argument("--summary", required=True, help="cxsv_summary.tsv")
    parser.add_argument("--output",  required=True, help="Output CxSV-only VCF")
    args = parser.parse_args()

    regions = load_cxsv_regions(args.summary)
    if not regions:
        sys.stderr.write("WARNING: No CxSV regions found in summary. Output will contain only the header.\n")

    n_written = 0
    with open(args.vcf) as vcf_in, open(args.output, "w") as vcf_out:
        for line in vcf_in:
            # Always write header lines
            if line.startswith("#"):
                vcf_out.write(line)
                continue

            parts  = line.split("\t")
            chrom  = parts[0]
            pos    = int(parts[1])
            info   = parts[7]

            # Parse END from INFO; fall back to POS + SVLEN
            end = pos
            for field in info.split(";"):
                if field.startswith("END="):
                    try:
                        end = int(field[4:])
                    except ValueError:
                        pass
                    break
            if end == pos:
                for field in info.split(";"):
                    if field.startswith("SVLEN="):
                        try:
                            svlen = abs(int(field[6:]))
                            end   = pos + svlen
                        except ValueError:
                            pass
                        break

            if overlaps(chrom, pos, end, regions):
                vcf_out.write(line)
                n_written += 1

    print(f"Wrote {n_written} CxSV records to {args.output}.")


if __name__ == "__main__":
    main()
