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
    --long-span-rule paired     paired-breakend linkage: the record is written as
        its two breakends (as for "breakends") only when BOTH ends independently
        link to another record that takes part in discovery, i.e. lie within
        --cluster-distance bp of it (the bedtools cluster -d rule). Records that
        duplicate it (both ends within --collapse-distance of its two ends: the
        same junction called twice, e.g. a Sniffles2 BND pair rewritten as INV
        beside an SVIM INV) do not count. Long records are tested together until
        no more fail (a long partner must itself be paired). A record that fails
        is not written, as under "exclude". The two end clusters of a paired
        record are joined into one cluster afterwards (link_paired_clusters.py),
        so the record is an edge between them; its interval is never a clustering
        interval, and classify_cxsv.py sizes such a cluster by its parts.
  Either way the record is listed in --long-out (ID CHROM POS END SVTYPE LENGTH
  SUPP_VEC RULE; RULE is the rule applied: exclude, breakends, or for "paired"
  either paired or exclude). --span-link-max 0 restores the interval rule for
  every record.

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
      [--span-link-max 500000 --long-span-rule exclude|breakends|paired --long-out long_span_records.tsv]
      [--cluster-distance 50000 --collapse-distance 300]   (used by "paired" only)
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


LONG_RULES = ("exclude", "breakends", "paired")


def _near(row_start, row_end, pos, distance):
    """A BED row [row_start, row_end] lies within distance bp of point pos."""
    return max(0, row_start - pos, pos - row_end) <= distance


def paired_records(rows, long_rows, cluster_distance, collapse_distance):
    """IDs of the long records whose two ends each lie within cluster_distance of
    another record taking part in discovery (paired-breakend linkage).

    rows       BED lines ("CHROM\tSTART\tEND\tSVTYPE\tID\tIS_INTER") of every record,
               long records written as their two breakend rows (ID, ID_end)
    long_rows  (ID, CHROM, POS, END, ...) of the long records

    A record whose rows touch both ends of the long record within collapse_distance
    is the same junction (a duplicate) and is no partner. A long partner counts only
    while it is itself paired: long records that fail are removed and the rest are
    tested again until none fails.
    """
    by_chrom = {}
    points = {}                                   # canonical id -> [(chrom, pos)]
    for line in rows:
        c, s, e, _t, vid, _i = line.split("\t")[:6]
        rid = _canonical(vid)
        by_chrom.setdefault(c, []).append((int(s), int(e), rid))
        points.setdefault(rid, []).extend([(c, int(s)), (c, int(e))])
    for c in by_chrom:
        by_chrom[c].sort()
    longs = {r[0]: (r[1], int(r[2]), int(r[3])) for r in long_rows}

    def duplicate(rid, chrom, p1, p2):
        pts = points.get(rid, ())
        return (any(c == chrom and abs(p - p1) <= collapse_distance for c, p in pts)
                and any(c == chrom and abs(p - p2) <= collapse_distance for c, p in pts))

    def linked(lid, chrom, p, other_end, alive):
        for s, e, rid in by_chrom.get(chrom, ()):
            if s > p + cluster_distance:
                break
            if rid == lid or not _near(s, e, p, cluster_distance):
                continue
            if rid in longs and rid not in alive:
                continue
            if duplicate(rid, chrom, p, other_end):
                continue
            return True
        return False

    alive = set(longs)
    while True:
        failed = {lid for lid in alive
                  if not (linked(lid, longs[lid][0], longs[lid][1], longs[lid][2], alive)
                          and linked(lid, longs[lid][0], longs[lid][2], longs[lid][1], alive))}
        if not failed:
            return alive
        alive -= failed


def _canonical(vid):
    for suf in ("_mate", "_dest", "_end"):
        if vid.endswith(suf):
            return vid[: -len(suf)]
    return vid


def process_vcf(vcf_path: str, window: int, span_link_max: int = 0, counts: dict = None,
                long_rule: str = "breakends", long_rows: list = None,
                cluster_distance: int = 50000, collapse_distance: int = 300):
    """
    Generator: yield BED lines for each VCF record (see the module docstring).
    long_rule "paired" reads the whole VCF first (paired_records), then yields
    every row except those of the long records that are not paired.
    """
    if long_rule not in LONG_RULES:
        raise ValueError(f"long_rule must be one of {LONG_RULES}, not {long_rule!r}")
    if long_rule != "paired":
        yield from _process_vcf(vcf_path, window, span_link_max, counts, long_rule, long_rows)
        return
    lrows = []
    rows = list(_process_vcf(vcf_path, window, span_link_max, None, "breakends", lrows))
    keep = paired_records(rows, lrows, cluster_distance, collapse_distance)
    for r in lrows:
        rule = "paired" if r[0] in keep else "exclude"
        if counts is not None:
            counts[rule] = counts.get(rule, 0) + 1
        if long_rows is not None:
            long_rows.append(tuple(r[:7]) + (rule,))
    drop = {r[0] for r in lrows} - keep
    for line in rows:
        f = line.split("\t")
        if f[3] in ("DEL", "DUP", "INV") and _canonical(f[4]) in drop:
            continue
        yield line


def _process_vcf(vcf_path: str, window: int, span_link_max: int = 0, counts: dict = None,
                 long_rule: str = "breakends", long_rows: list = None):
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
                             "their two breakends; paired: write their two breakends only when "
                             "both ends link to other records (cxsv_params.long_span_rule)")
    parser.add_argument("--cluster-distance", type=int, default=50000,
                        help="paired: an end links to a record within this many bp (cxsv_params.cluster_distance)")
    parser.add_argument("--collapse-distance", type=int, default=300,
                        help="paired: a record with rows within this many bp of both ends is the same "
                             "junction, not a partner (junction_params.collapse_distance)")
    parser.add_argument("--long-out", default=None,
                        help="TSV listing the records longer than --span-link-max and the rule applied")
    args = parser.parse_args()

    n_written = 0
    counts, long_rows = {}, []
    with open(args.out, "w") as out_fh:
        for bed_line in process_vcf(args.vcf, args.window, args.span_link_max, counts,
                                    args.long_span_rule, long_rows, args.cluster_distance,
                                    args.collapse_distance):
            out_fh.write(bed_line + "\n")
            n_written += 1
    if args.long_out:
        with open(args.long_out, "w") as fh:
            fh.write("ID\tCHROM\tPOS\tEND\tSVTYPE\tLENGTH\tSUPP_VEC\tRULE\n")
            for r in long_rows:
                fh.write("\t".join(map(str, r)) + "\n")

    print(f"vcf_to_bed: wrote {n_written} BED records to {args.out}", file=sys.stderr)
    print(f"vcf_to_bed: {len(long_rows)} DEL/DUP/INV records longer than {args.span_link_max} bp "
          f"({args.long_span_rule}: " + ", ".join(f"{k} {v}" for k, v in sorted(counts.items())) + ")",
          file=sys.stderr)


if __name__ == "__main__":
    main()
