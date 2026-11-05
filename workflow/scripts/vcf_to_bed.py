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

Long span records (roadmap item 13)
-----------------------------------
  Written as one interval, a DEL/DUP/INV longer than cxsv_params.max_cluster_size
  joins every SV it spans into one cluster that classify_cxsv.py then drops whole
  as oversized: a 25.7 Mb record across the chr17 pericentromere put 224 records
  into one 25.8 Mb cluster. For records longer than --span-link-max bp
  (= max_cluster_size, so no cluster that is not oversized changes):
    --long-span-rule exclude    the record is not
        written: it takes no part in CxSV discovery, as before, when any cluster
        containing it was dropped; it stays a junction of the sample
        (junctions.tsv) and a merged call (merged_sv.vcf)
    --long-span-rule breakends  the record is written as its two breakends (POS
        and END +- --window bp; the END row's ID gets the suffix "_end"), like a
        BND, so it links only to SVs near its own breakpoints and is a junction
        of the loci there (classify_cxsv.py strips "_end" like "_mate")
  Either way the record is listed in --long-out (ID CHROM POS END SVTYPE LENGTH
  SUPP_VEC RULE). --span-link-max 0 restores the interval rule for every record.

Single-caller breakend filter (cxsv_params.single_caller_bnd_rule, docs/R10_EXCESS.md)
--------------------------------------------------------------------------------------
  A few panel genomes carry thousands of SVIM-only translocation records at
  SVIM's quality floor (QUAL 5): 3,753-14,258 in 7 LRSC R10 genomes against a
  median of 64 in the other R10 genomes. Clustered with real SVs they make
  BND_cluster / CN_plus_orientation_change loci out of simple deletions.
  --single-caller-bnd-rule decides which merged records take no part in CxSV
  discovery (they stay junctions and merged calls, like long-span records):
    keep            none (command-line default; discovery as in panel-v1.0)
    svim_lowq       SVIM-only (SUPP_VEC 010) BND/TRA records with QUAL < --svim-min-qual
                    (the pipeline default since 2026-11-05, config/calling.yaml)
    svim_lowq_any   SVIM-only records of any type with QUAL < --svim-min-qual
    single_bnd      every single-caller BND/TRA record, whatever its quality
  SUPP_VEC order is the SURVIVOR file order: cuteSV, SVIM, Sniffles2.
  Excluded records are listed in --filtered-out.

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
  python vcf_to_bed.py --vcf input.vcf --out output.bed [--window 1] \
      [--span-link-max 500000 --long-span-rule exclude --long-out long_span_records.tsv]
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


LONG_RULES = ("exclude", "breakends")
SINGLE_BND_RULES = ("keep", "svim_lowq", "svim_lowq_any", "single_bnd")
SVIM_ONLY = "010"          # SURVIVOR file order: cuteSV, SVIM, Sniffles2


def single_caller_excluded(rule, svtype, supp_vec, qual, svim_min_qual):
    """True when --single-caller-bnd-rule leaves this merged record out of discovery."""
    if rule == "keep":
        return False
    bndlike = svtype in ("BND", "TRA")
    if rule == "single_bnd":
        return bndlike and supp_vec.count("1") == 1
    lowq = qual is not None and qual < svim_min_qual
    if rule == "svim_lowq":
        return bndlike and supp_vec == SVIM_ONLY and lowq
    if rule == "svim_lowq_any":
        return supp_vec == SVIM_ONLY and lowq
    raise ValueError(f"single-caller rule must be one of {SINGLE_BND_RULES}, not {rule!r}")


def process_vcf(vcf_path: str, window: int, span_link_max: int = 0, counts: dict = None,
                long_rule: str = "breakends", long_rows: list = None,
                single_rule: str = "keep", svim_min_qual: float = 10, filtered_rows: list = None):
    """
    Generator: yield BED lines for each VCF record.
    BND/TRA yield 2 lines each (both breakend positions). A DEL/DUP/INV longer
    than span_link_max (> 0) yields 2 breakend lines (long_rule "breakends") or
    none ("exclude"); counts[long_rule] counts them and long_rows, if given,
    receives (ID, CHROM, POS, END, SVTYPE, LENGTH, SUPP_VEC, RULE) for each.
    """
    if long_rule not in LONG_RULES:
        raise ValueError(f"long_rule must be one of {LONG_RULES}, not {long_rule!r}")
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
            if single_rule != "keep":
                try:
                    qual = float(parts[5])
                except ValueError:
                    qual = None
                supp_vec = str(info.get("SUPP_VEC", ""))
                if single_caller_excluded(single_rule, svtype, supp_vec, qual, svim_min_qual):
                    if counts is not None:
                        counts["single_caller_filtered"] = counts.get("single_caller_filtered", 0) + 1
                    if filtered_rows is not None:
                        filtered_rows.append((vid, chrom, pos, svtype, supp_vec, parts[5], single_rule))
                    continue

            # ── Standard SV types ──────────────────────────────────────────
            if svtype in ("DEL", "DUP", "INV", "INS"):
                coords = resolve_standard(chrom, pos, info, svtype)
                if coords is None:
                    sys.stderr.write(
                        f"WARN: Cannot resolve END for {vid} ({chrom}:{pos}), skipping.\n"
                    )
                    continue
                start, end = coords
                if span_link_max > 0 and svtype != "INS" and end - start > span_link_max:
                    if counts is not None:
                        counts[long_rule] = counts.get(long_rule, 0) + 1
                    if long_rows is not None:
                        long_rows.append((vid, chrom, start, end, svtype, end - start,
                                          info.get("SUPP_VEC", "."), long_rule))
                    if long_rule == "exclude":
                        continue
                    yield from emit_records(chrom, max(0, start - window), start + window, svtype, vid, "intra")
                    yield from emit_records(chrom, max(0, end - window), end + window, svtype,
                                            f"{vid}_end", "intra")
                    continue
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
    parser.add_argument("--span-link-max", type=int, default=0,
                        help="DEL/DUP/INV longer than this many bp follow --long-span-rule "
                             "(cxsv_params.span_link_max; 0 = every record is one interval)")
    parser.add_argument("--long-span-rule", choices=LONG_RULES, default="breakends",
                        help="exclude: leave such records out of clustering; breakends: write "
                             "their two breakends (cxsv_params.long_span_rule)")
    parser.add_argument("--long-out", default=None,
                        help="TSV listing the records longer than --span-link-max and the rule applied")
    parser.add_argument("--single-caller-bnd-rule", choices=SINGLE_BND_RULES, default="keep",
                        help="single-caller breakend records left out of discovery (see the module docstring; the pipeline passes cxsv_params.single_caller_bnd_rule, default svim_lowq)")
    parser.add_argument("--svim-min-qual", type=float, default=10,
                        help="QUAL below which an SVIM-only record is low quality (svim_lowq rules)")
    parser.add_argument("--filtered-out", default=None,
                        help="TSV listing the records left out by --single-caller-bnd-rule")
    args = parser.parse_args()

    n_written = 0
    counts, long_rows = {}, []
    with open(args.out, "w") as out_fh:
        filtered_rows = []
        for bed_line in process_vcf(args.vcf, args.window, args.span_link_max, counts,
                                    args.long_span_rule, long_rows, args.single_caller_bnd_rule,
                                    args.svim_min_qual, filtered_rows):
            out_fh.write(bed_line + "\n")
            n_written += 1
    if args.long_out:
        with open(args.long_out, "w") as fh:
            fh.write("ID\tCHROM\tPOS\tEND\tSVTYPE\tLENGTH\tSUPP_VEC\tRULE\n")
            for r in long_rows:
                fh.write("\t".join(map(str, r)) + "\n")

    if args.filtered_out:
        with open(args.filtered_out, "w") as fh:
            fh.write("ID\tCHROM\tPOS\tSVTYPE\tSUPP_VEC\tQUAL\tRULE\n")
            for r in filtered_rows:
                fh.write("\t".join(map(str, r)) + "\n")
    print(f"vcf_to_bed: wrote {n_written} BED records to {args.out}", file=sys.stderr)
    print(f"vcf_to_bed: {len(filtered_rows)} records left out of discovery by the single-caller rule "
          f"{args.single_caller_bnd_rule}", file=sys.stderr)
    print(f"vcf_to_bed: {len(long_rows)} DEL/DUP/INV records longer than {args.span_link_max} bp "
          f"({args.long_span_rule})", file=sys.stderr)


if __name__ == "__main__":
    main()
