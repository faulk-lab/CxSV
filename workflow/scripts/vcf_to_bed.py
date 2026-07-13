#!/usr/bin/env python3
"""
vcf_to_bed.py
=============
Convert a merged SV VCF to a 6-column BED file suitable for bedtools
clustering, with full support for BND and TRA records.

Standard SV types (DEL, DUP, INV, INS)
---------------------------------------
  Coordinates come from POS and INFO/END.
  If END < POS (a known bcftools quirk for some callers), coordinates
  are swapped. If END is missing, INFO/SVLEN is used as a fallback.

BND records
-----------
  VCF BND ALT fields encode the mate position in one of four bracket
  notations:
      ]chr:pos]N   N[chr:pos[   [chr:pos[N   N]chr:pos]
  We parse the mate chromosome and position from the ALT field and emit
  TWO BED records per BND pair:
    1. The breakend position itself (POS ± a small window)
    2. The mate position (mate_chr, mate_pos ± window)
  This ensures both ends of the translocation/rearrangement are
  represented in the clustering step so inter-chromosomal CxSV loci
  are captured.

TRA records
-----------
  Some callers (e.g. SVIM) emit TRA with CHR2 and END2 in INFO rather
  than the bracket notation. We parse CHR2/END (or END2) and emit two
  records as for BND.

Output columns
--------------
  1  CHROM      chromosome of this interval
  2  START      0-based start
  3  END        0-based end
  4  SVTYPE     canonical type (DEL, DUP, INV, INS, BND, TRA)
  5  ID         variant ID
  6  IS_INTER   "inter" if the two ends are on different chromosomes, else "intra"

Usage
-----
  python vcf_to_bed.py --vcf input.vcf --out output.bed [--window 1]
"""

import argparse
import re
import sys

# bp added on each side of a BND/TRA breakpoint for interval representation
DEFAULT_WINDOW = 1


def parse_info(info_str: str) -> dict:
    """Parse INFO field into a dict. Flags get value True."""
    d = {}
    for field in info_str.split(";"):
        if "=" in field:
            k, v = field.split("=", 1)
            d[k] = v
        else:
            d[field] = True
    return d


def parse_bnd_alt(alt: str):
    """
    Parse a BND ALT string and return (mate_chr, mate_pos) or (None, None).
    Handles all four VCF BND notations:
      ]chr:pos]N   [chr:pos[N   N]chr:pos]   N[chr:pos[
    """
    m = re.search(r'[\[\]]([^[\]]+):(\d+)[\[\]]', alt)
    if m:
        return m.group(1), int(m.group(2))
    return None, None


def resolve_standard(chrom, pos, info, svtype):
    """
    Return (start, end) for DEL/DUP/INV/INS.
    Applies END swap fix and SVLEN fallback.
    Returns None if coordinates cannot be resolved.
    """
    pos = int(pos)
    end_raw = info.get("END", ".")
    svlen_raw = info.get("SVLEN", ".")

    if end_raw not in (".", "", "0"):
        try:
            end = int(end_raw)
        except ValueError:
            end = None
    else:
        end = None

    if end is None:
        # Fallback: use absolute SVLEN
        if svlen_raw not in (".", ""):
            try:
                end = pos + abs(int(svlen_raw.split(",")[0]))
            except ValueError:
                return None
        else:
            return None

    # Swap if inverted (bcftools quirk)
    if end < pos:
        pos, end = end, pos

    if end <= pos:
        return None
    if end - pos > 500_000_000:
        return None

    return pos, end


def emit_records(chrom, start, end, svtype, vid, is_inter="intra"):
    """Yield a single BED record as a tab-separated string."""
    yield f"{chrom}\t{start}\t{end}\t{svtype}\t{vid}\t{is_inter}"


def process_vcf(vcf_path: str, window: int):
    """
    Generator: yield BED lines for each VCF record.
    BND/TRA yield 2 lines each (both breakend positions).
    """
    with open(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue

            chrom   = parts[0]
            pos_str = parts[1]
            vid     = parts[2]
            alt     = parts[4]
            info_str = parts[7]

            pos  = int(pos_str)
            info = parse_info(info_str)
            svtype = info.get("SVTYPE", "UNK").upper()

            # ── Standard SV types ──────────────────────────────────────────
            if svtype in ("DEL", "DUP", "INV", "INS"):
                coords = resolve_standard(chrom, pos, info, svtype)
                if coords is None:
                    sys.stderr.write(
                        f"WARN: Cannot resolve END for {vid} ({chrom}:{pos}), skipping.\n"
                    )
                    continue
                start, end = coords
                yield from emit_records(chrom, start, end, svtype, vid, "intra")

            # ── BND records ────────────────────────────────────────────────
            elif svtype == "BND":
                mate_chr, mate_pos = parse_bnd_alt(alt)

                # Record 1: the breakend itself (window around POS)
                b1_start = max(0, pos - window)
                b1_end   = pos + window
                is_inter = "inter" if (mate_chr and mate_chr != chrom) else "intra"
                yield from emit_records(chrom, b1_start, b1_end, "BND", vid, is_inter)

                # Record 2: mate position
                if mate_chr and mate_pos:
                    b2_start = max(0, mate_pos - window)
                    b2_end   = mate_pos + window
                    yield from emit_records(
                        mate_chr, b2_start, b2_end, "BND", f"{vid}_mate", is_inter
                    )
                else:
                    sys.stderr.write(
                        f"WARN: Could not parse mate from BND ALT for {vid}: {alt}\n"
                    )

            # ── TRA records ────────────────────────────────────────────────
            elif svtype == "TRA":
                # SVIM uses CHR2 + END; some other callers use CHR2 + END2
                chr2     = info.get("CHR2", info.get("CHROM2", None))
                end2_raw = info.get("END2", info.get("END", None))

                # Fall back to BND-style ALT parsing if INFO fields missing
                if chr2 is None or end2_raw is None:
                    mate_chr, mate_pos = parse_bnd_alt(alt)
                    chr2  = mate_chr
                    end2_raw = str(mate_pos) if mate_pos else None

                # Record 1: the TRA origin
                b1_start = max(0, pos - window)
                b1_end   = pos + window
                is_inter = "inter" if (chr2 and chr2 != chrom) else "intra"
                yield from emit_records(chrom, b1_start, b1_end, "TRA", vid, is_inter)

                # Record 2: the TRA destination
                if chr2 and end2_raw:
                    try:
                        dest_pos = int(end2_raw)
                        b2_start = max(0, dest_pos - window)
                        b2_end   = dest_pos + window
                        yield from emit_records(
                            chr2, b2_start, b2_end, "TRA", f"{vid}_dest", is_inter
                        )
                    except ValueError:
                        sys.stderr.write(
                            f"WARN: Cannot parse TRA destination for {vid}, skipping dest record.\n"
                        )
                else:
                    sys.stderr.write(
                        f"WARN: No CHR2/END2 found for TRA record {vid} ({chrom}:{pos}).\n"
                    )

            else:
                # Unknown type — emit a window around POS as best effort
                b_start = max(0, pos - window)
                b_end   = pos + window
                yield from emit_records(chrom, b_start, b_end, svtype, vid, "intra")


def main():
    parser = argparse.ArgumentParser(description="Convert SV VCF to BED with BND/TRA support.")
    parser.add_argument("--vcf",    required=True, help="Input VCF")
    parser.add_argument("--out",    required=True, help="Output BED")
    parser.add_argument("--window", type=int, default=DEFAULT_WINDOW,
                        help=f"bp window around BND/TRA breakpoints (default: {DEFAULT_WINDOW})")
    args = parser.parse_args()

    n_written = 0
    with open(args.out, "w") as out_fh:
        for bed_line in process_vcf(args.vcf, args.window):
            out_fh.write(bed_line + "\n")
            n_written += 1

    print(f"vcf_to_bed: wrote {n_written} BED records to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
