#!/usr/bin/env python3
"""
filter_vcf_by_bed.py
=====================
Keep or exclude VCF records that overlap regions in a BED file, with an
optional bp window padding each region before the overlap check.

Used by workflow/rules/private_filter.smk to derive each query cohort's
per-sample private CxSV VCF: records in {sample}_cxsv_only.vcf that do NOT
overlap (within the configured slop window) any locus in the pooled public
reference master BED are kept as candidate pathogenic/private calls.

Usage
-----
python filter_vcf_by_bed.py \
    --vcf    results/patient_batchA/per_sample/SAMPLE_cxsv_only.vcf \
    --bed    results/_public_reference/master_cxsv.bed \
    --mode   exclude \
    --window 1000 \
    --output results/patient_batchA/per_sample/SAMPLE_private_cxsv.vcf
"""

import argparse
import sys


def load_regions(bed_path: str, window: int = 0) -> list[tuple[str, int, int]]:
    """
    Parse a BED file into a list of (chrom, start, end) tuples, each padded
    by `window` bp on both sides (clamped at 0). Only the first three
    columns are used; anything else in the BED is ignored.
    """
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


def vcf_record_span(parts: list[str]) -> tuple[str, int, int]:
    """Extract (chrom, pos, end) from a split VCF data line, using END or SVLEN."""
    chrom = parts[0]
    pos   = int(parts[1])
    info  = parts[7]

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
                    end = pos + abs(int(field[6:]))
                except ValueError:
                    pass
                break
    return chrom, pos, end


def main():
    parser = argparse.ArgumentParser(
        description="Keep or exclude VCF records overlapping a BED file."
    )
    parser.add_argument("--vcf",    required=True, help="Input VCF file")
    parser.add_argument("--bed",    required=True, help="BED file of regions to compare against")
    parser.add_argument("--mode",   required=True, choices=("keep", "exclude"),
                        help="keep: write records overlapping the BED; "
                             "exclude: write records NOT overlapping it")
    parser.add_argument("--window", type=int, default=0,
                        help="bp padding added to each side of every BED region (default 0)")
    parser.add_argument("--output", required=True, help="Output VCF file")
    args = parser.parse_args()

    regions = load_regions(args.bed, window=args.window)
    if not regions:
        sys.stderr.write(f"WARNING: No regions loaded from {args.bed}.\n")

    n_written = 0
    n_total   = 0
    with open(args.vcf) as vcf_in, open(args.output, "w") as vcf_out:
        for line in vcf_in:
            if line.startswith("#"):
                vcf_out.write(line)
                continue

            parts = line.split("\t")
            if len(parts) < 8:
                continue
            n_total += 1

            chrom, pos, end = vcf_record_span(parts)
            hit = overlaps(chrom, pos, end, regions)
            keep = hit if args.mode == "keep" else not hit

            if keep:
                vcf_out.write(line)
                n_written += 1

    print(f"[{args.mode}] wrote {n_written}/{n_total} records to {args.output} "
          f"(window={args.window} bp, {len(regions)} regions from {args.bed}).")


if __name__ == "__main__":
    main()
