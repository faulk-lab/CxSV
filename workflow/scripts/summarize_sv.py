#!/usr/bin/env python3
"""
summarize_sv.py
===============
Per-sample SV summary report. Produces a human-readable .txt file with:

  SECTION 1 — SV counts by type per caller
              (DEL / DUP / INV / INS / BND / TRA / UNKNOWN)
              Columns: Sniffles2 | CuteSV | SVIM | Merged(>=1 caller)

  SECTION 2 — CxSV summary for this sample
              - Total CxSV loci
              - Breakdown by which criteria fired (C1-C5)
              - Breakdown by complexity class (DEL+INV, BND_cluster, etc.)

  SECTION 3 — Caller concordance from SURVIVOR SUPP_VEC
              How many merged SVs were called by 1, 2, or all 3 callers.
              SUPP_VEC bit order: CuteSV=bit0, SVIM=bit1, Sniffles2=bit2.

Important note on genomic context:
  CxSVs are NOT filtered by genomic location. Intergenic, repeat-rich, and
  segmental duplication loci are INCLUDED -- these grey-zone regions are where
  chromoanagenesis mechanisms (FoSTeS/MMBIR, chromothripsis) preferentially
  act. The Genomic_Context column in the master table is annotation-only.

Usage
-----
python summarize_sv.py
    --merged   results/per_sample/SAMPLE_merged_sv.vcf
    --cutesv   results/callers/SAMPLE_cutesv_filtered.vcf
    --sniffles results/callers/SAMPLE_sniffles_filtered.vcf
    --svim     results/callers/SAMPLE_svim_filtered.vcf
    --cxsv     results/per_sample/SAMPLE_cxsv_only.vcf
    --summary  results/final/cxsv_summary.tsv
    --sample   SAMPLE_NAME
    --output   results/per_sample/SAMPLE_sv_summary.txt
"""

import argparse
import sys
from collections import defaultdict


def count_sv_types(vcf_path):
    counts = defaultdict(int)
    try:
        with open(vcf_path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 8:
                    continue
                svtype = "UNKNOWN"
                for field in parts[7].split(";"):
                    if field.startswith("SVTYPE="):
                        svtype = field[7:].strip()
                        break
                counts[svtype] += 1
    except FileNotFoundError:
        pass
    return dict(counts)


def parse_supp_vec(vcf_path):
    """Parse SUPP_VEC from SURVIVOR merged VCF. Bit order: CuteSV=0 SVIM=1 Sniffles2=2."""
    counts = defaultdict(int)
    try:
        with open(vcf_path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 8:
                    continue
                vec = None
                for field in parts[7].split(";"):
                    if field.startswith("SUPP_VEC="):
                        vec = field[9:].strip()
                        break
                if vec:
                    n = vec.count("1")
                    label = f"{n}_caller{'s' if n != 1 else ''}"
                    counts[label] += 1
                else:
                    counts["no_SUPP_VEC"] += 1
    except FileNotFoundError:
        pass
    return dict(counts)


def classify_complexity(sv_types_str, n_bp, criteria_str):
    """Complexity class — must stay in sync with annotate_cxsv.py."""
    types    = set(sv_types_str.upper().split(","))
    criteria = set(criteria_str.split(","))
    if n_bp >= 10 or "C5_chromoanagenesis_cluster" in criteria:
        return "Chromothripsis_like"
    if {"DEL","INV","DUP"}.issubset(types): return "DEL+DUP+INV"
    if {"DEL","INV"}.issubset(types):       return "DEL+INV"
    if {"DUP","INV"}.issubset(types):       return "DUP+INV"
    if {"INS","DEL"}.issubset(types):       return "INS+DEL"
    if "BND" in types and n_bp >= 3:        return "BND_cluster"
    if n_bp >= 5:                           return "Multi_breakpoint_cluster"
    if "C2_nested_SV"         in criteria:  return "Nested_SV"
    if "C3_overlapping_types" in criteria:  return "Overlapping_types"
    if "C4_CN_plus_orient"    in criteria:  return "CN_plus_orientation_change"
    return "Complex_other"


def load_cxsv_for_sample(summary_path, sample_cxsv_vcf):
    """Cross-reference population cxsv_summary.tsv with per-sample cxsv VCF."""
    cxsv_regions = []
    try:
        with open(summary_path) as fh:
            cols = fh.readline().rstrip("\n").split("\t")
            try:
                ix = {c: i for i, c in enumerate(cols)}
                ci_chr      = ix["Chr"]
                ci_start    = ix["Start"]
                ci_end      = ix["End"]
                ci_is_cxsv  = ix["IsCxSV"]
                ci_criteria = ix["CxSV_Criteria"]
                ci_types    = ix["SV_Types"]
                ci_nbp      = ix["N_Breakpoints"]
            except KeyError as e:
                sys.stderr.write(f"WARN: missing column in cxsv_summary.tsv: {e}\n")
                return 0, {}, {}
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= ci_is_cxsv or parts[ci_is_cxsv] != "YES":
                    continue
                cxsv_regions.append((
                    parts[ci_chr], int(parts[ci_start]), int(parts[ci_end]),
                    parts[ci_criteria], parts[ci_types], int(parts[ci_nbp])
                ))
    except FileNotFoundError:
        return 0, {}, {}

    if not cxsv_regions:
        return 0, {}, {}

    sample_hits = set()
    try:
        with open(sample_cxsv_vcf) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 8:
                    continue
                chrom = parts[0]
                pos   = int(parts[1])
                end   = pos
                for field in parts[7].split(";"):
                    if field.startswith("END="):
                        try:
                            end = int(field[4:])
                        except ValueError:
                            pass
                        break
                for i, region in enumerate(cxsv_regions):
                    rc, rs, re = region[0], region[1], region[2]
                    if rc == chrom and pos <= re and end >= rs:
                        sample_hits.add(i)
    except FileNotFoundError:
        pass

    criteria_counts = defaultdict(int)
    class_counts    = defaultdict(int)
    for i in sample_hits:
        rc, rs, re, criteria_str, sv_types, n_bp = cxsv_regions[i]
        for c in criteria_str.split(","):
            if c and c != "NONE":
                criteria_counts[c] += 1
        cmplx = classify_complexity(sv_types, n_bp, criteria_str)
        class_counts[cmplx] += 1

    return len(sample_hits), dict(criteria_counts), dict(class_counts)


CRITERIA_LABELS = {
    "C1_min_breakpoints":         "C1  >=3 breakpoints in 50 kb window",
    "C2_nested_SV":               "C2  Nested SV (smaller inside larger)",
    "C3_overlapping_types":       "C3  Overlapping SV types at same locus",
    "C4_CN_plus_orient":          "C4  Copy-number shift + orientation change",
    "C5_chromoanagenesis_cluster":"C5  Chromoanagenesis cluster (10-50 kb)",
}


def rule(char="-", width=74):
    return char * width


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--merged",   required=True)
    parser.add_argument("--cutesv",   required=True)
    parser.add_argument("--sniffles", required=True)
    parser.add_argument("--svim",     required=True)
    parser.add_argument("--cxsv",     required=True)
    parser.add_argument("--summary",  required=True)
    parser.add_argument("--sample",   required=True)
    parser.add_argument("--output",   required=True)
    args = parser.parse_args()

    caller_counts = {
        "Sniffles2":          count_sv_types(args.sniffles),
        "CuteSV":             count_sv_types(args.cutesv),
        "SVIM":               count_sv_types(args.svim),
        "Merged (>=1 caller)":count_sv_types(args.merged),
    }
    all_types   = sorted({t for c in caller_counts.values() for t in c})
    supp_counts = parse_supp_vec(args.merged)
    total_cxsv, criteria_counts, class_counts = load_cxsv_for_sample(
        args.summary, args.cxsv
    )
    cw = 22

    with open(args.output, "w") as out:

        out.write(rule("=") + "\n")
        out.write(f"  SV SUMMARY REPORT -- Sample: {args.sample}\n")
        out.write(rule("=") + "\n\n")

        # SECTION 1
        out.write("SECTION 1 -- SV Counts by Type and Caller\n")
        out.write(rule() + "\n")
        hdr = f"{'SV Type':<14}" + "".join(f"{k:>{cw}}" for k in caller_counts)
        out.write(hdr + "\n")
        out.write(rule("-") + "\n")
        for svtype in all_types:
            row = f"{svtype:<14}"
            row += "".join(
                f"{caller_counts[k].get(svtype, 0):>{cw},}" for k in caller_counts
            )
            out.write(row + "\n")
        out.write(rule("-") + "\n")
        tot_row = f"{'TOTAL':<14}"
        tot_row += "".join(
            f"{sum(caller_counts[k].values()):>{cw},}" for k in caller_counts
        )
        out.write(tot_row + "\n\n")

        # SECTION 2
        out.write("SECTION 2 -- CxSV Summary\n")
        out.write(rule() + "\n")
        out.write(f"  Total CxSV loci in this sample: {total_cxsv}\n")
        out.write(f"  (Cross-referenced against population cxsv_summary.tsv)\n\n")

        out.write("  Criteria breakdown (loci can satisfy multiple criteria):\n")
        out.write(f"  {'Criterion':<52} {'Count':>6}\n")
        out.write(f"  {rule('-', 58)}\n")
        for key, label in CRITERIA_LABELS.items():
            n = criteria_counts.get(key, 0)
            out.write(f"  {label:<52} {n:>6}\n")

        out.write(f"\n  Complexity class breakdown:\n")
        out.write(f"  {'Class':<45} {'Count':>6}\n")
        out.write(f"  {rule('-', 51)}\n")
        if class_counts:
            for cmplx, n in sorted(class_counts.items(), key=lambda x: -x[1]):
                out.write(f"  {cmplx:<45} {n:>6}\n")
        else:
            out.write("  No CxSV loci detected for this sample.\n")
        out.write("\n")

        # SECTION 3
        out.write("SECTION 3 -- Caller Concordance (SURVIVOR SUPP_VEC)\n")
        out.write(rule() + "\n")
        out.write("  Bit order in SUPP_VEC: bit0=CuteSV, bit1=SVIM, bit2=Sniffles2\n\n")
        for key in ["1_caller", "2_callers", "3_callers", "no_SUPP_VEC"]:
            n = supp_counts.get(key, 0)
            if n > 0:
                out.write(f"  {'Supported by ' + key:<35} {n:>7,}\n")
        if not supp_counts:
            out.write("  Could not parse SUPP_VEC (check merged VCF header).\n")
        out.write("\n")

        out.write(rule("-") + "\n")
        out.write("  CxSV criteria: >=3 bp (C1) | nested SV (C2) | overlapping\n")
        out.write("  types (C3) | CN+orient change (C4) | 10-50kb cluster (C5)\n")
        out.write("  Intergenic / repeat / segdup loci are INCLUDED -- most\n")
        out.write("  CxSVs are expected in non-coding grey-zone regions.\n")
        out.write(rule("-") + "\n")

    print(f"Summary written to {args.output}")


if __name__ == "__main__":
    main()
