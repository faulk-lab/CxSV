#!/usr/bin/env python3
"""
annotate_cxsv.py
================
Build the CxSV master table and population BED file.

IMPORTANT — no filtering by genomic context
--------------------------------------------
CxSVs are annotated with their genomic context (exonic, intronic, repeat,
segdup, intergenic) for informational purposes ONLY. No loci are removed
based on context. The majority of CxSVs are expected in intergenic, repeat-
rich, and segmental duplication regions — the 'grey zone' of the genome
where chromoanagenesis mechanisms (FoSTeS/MMBIR, chromothripsis) act.

Master table columns (tab-separated, with header)
--------------------------------------------------
 1  Chr
 2  Start
 3  End
 4  Size
 5  Breakpoint_1
 6  Breakpoint_2
 7  SV_Types
 8  Cmplx_Class
 9  CxSV_Criteria          — which of C1-C5 fired
10  Seq_Type               — always "long" for this workflow
11  Has_Inter_Chromosomal  — YES/NO from classify_cxsv
12  Genomic_Context        — annotation-only label (see above)
13  Gene_Overlap           — gene name(s) or NA
14  Nearest_Gene
15  Distance_To_Nearest_Gene
16  Exon_Overlap           — yes/no
17  Repeat_Annotation      — repeat class(es) or NA
18  SegDup_Overlap         — fraction 0.0-1.0

Population BED (4 columns)
---------------------------
 1  Chr
 2  Start
 3  End
 4  Exon_Disrupt  — Y or N
"""

import argparse
import subprocess
import sys
import os
import re
import tempfile
from collections import defaultdict


# ─── bedtools helpers ─────────────────────────────────────────────────────────

def bedtools_intersect(a_bed: str, b_bed: str, extra: str = "") -> list:
    cmd = f"bedtools intersect -a {a_bed} -b {b_bed} {extra}"
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return r.stdout.strip().split("\n") if r.stdout.strip() else []


def bedtools_closest(a_bed: str, b_bed: str) -> list:
    cmd = f"bedtools sort -i {a_bed} | bedtools closest -a stdin -b {b_bed} -D ref"
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return r.stdout.strip().split("\n") if r.stdout.strip() else []


def coverage_fraction(a_bed: str, b_bed: str) -> dict:
    """Return {col4_id -> overlap_fraction} using bedtools coverage."""
    cmd = f"bedtools coverage -a {a_bed} -b {b_bed}"
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    fractions = {}
    for line in r.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) >= 4:
            fractions[parts[3]] = float(parts[-1])
    return fractions


# ─── Complexity classification ────────────────────────────────────────────────

def classify_complexity(sv_types_str: str, n_breakpoints: int,
                        criteria_str: str = "") -> str:
    """
    Assign a complexity class label.
    Uses the CxSV_Criteria string from classify_cxsv.py when available so
    C5 (chromoanagenesis cluster) is correctly labelled even when bp count < 10.
    Priority order reflects biological severity/complexity.
    """
    types    = set(sv_types_str.upper().replace(" ", "").split(","))
    n        = n_breakpoints
    criteria = set(criteria_str.split(",")) if criteria_str else set()

    if n >= 10 or "C5_chromoanagenesis_cluster" in criteria:
        return "Chromothripsis_like"
    if {"DEL", "INV", "DUP"}.issubset(types):
        return "DEL+DUP+INV"
    if {"DEL", "INV"}.issubset(types):
        return "DEL+INV"
    if {"DUP", "INV"}.issubset(types):
        return "DUP+INV"
    if {"INS", "DEL"}.issubset(types):
        return "INS+DEL"
    if "BND" in types and n >= 3:
        return "BND_cluster"
    if n >= 5:
        return "Multi_breakpoint_cluster"
    # Catch loci that passed via structural criteria even with low bp count
    if "C2_nested_SV"         in criteria: return "Nested_SV"
    if "C3_overlapping_types" in criteria: return "Overlapping_types"
    if "C4_CN_plus_orient"    in criteria: return "CN_plus_orientation_change"
    return "Complex_other"


# ─── Genomic context ──────────────────────────────────────────────────────────

def assign_genomic_context(uid: str,
                           gene_overlap: dict,
                           exon_coverage: dict,
                           nearest_dist: dict,
                           repeat_annot: dict,
                           segdup_frac: dict,
                           promoter_window: int = 2000,
                           exon_frac_threshold: float = 0.10) -> str:
    """
    Assign a single genomic context label for display/stratification.

    Priority: Exonic > Promoter > Intronic > Repeat_intergenic >
              SegDup_intergenic > Intergenic

    IMPORTANT CHANGE: Exonic requires >=10% of the CxSV locus to be
    covered by exon sequence (not just any 1bp overlap). This prevents
    large CxSV clusters from being labelled Exonic simply because they
    span a region that contains an exon somewhere within it.

    This label is ANNOTATION ONLY. It is never used to filter CxSVs.
    Most CxSVs are expected to be Repeat_intergenic or Intergenic
    because chromoanagenesis acts in the grey zone of the genome.
    """
    # Exonic: require substantial exon coverage fraction
    if exon_coverage.get(uid, 0.0) >= exon_frac_threshold:
        return "Exonic"

    has_gene = bool(gene_overlap.get(uid))
    dist_raw = nearest_dist.get(uid, "NA")
    try:
        dist_val = abs(int(dist_raw))
    except (ValueError, TypeError):
        dist_val = None

    if has_gene:
        # dist_val from bedtools closest -D ref:
        #   0        = locus overlaps the gene body (intronic)
        #   negative = locus is upstream of gene (potential promoter)
        #   positive = locus is downstream of gene
        # Re-read the raw signed distance for directional logic
        dist_raw = nearest_dist.get(uid, "NA")
        try:
            dist_signed = int(dist_raw)
        except (ValueError, TypeError):
            dist_signed = None

        # Promoter: upstream within promoter_window AND outside gene body
        if (dist_signed is not None
                and dist_signed < 0
                and abs(dist_signed) <= promoter_window):
            return "Promoter"
        # Everything else overlapping a gene body = Intronic
        return "Intronic"

    if repeat_annot.get(uid):
        return "Repeat_intergenic"

    if segdup_frac.get(uid, 0.0) > 0.5:
        return "SegDup_intergenic"

    return "Intergenic"


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Annotate CxSV loci and build master table + population BED."
    )
    parser.add_argument("--summary",   required=True,
                        help="cxsv_summary.tsv from classify_cxsv.py")
    parser.add_argument("--vcf",       required=True,
                        help="cxsv_only.vcf")
    parser.add_argument("--genes",     required=True,
                        help="Genes BED (col4 = gene_name)")
    parser.add_argument("--exons",     required=True,
                        help="Exons BED (col4 = gene_name)")
    parser.add_argument("--repeats",   required=True,
                        help="RepeatMasker BED (col4 = repeat_class)")
    parser.add_argument("--segdups",   required=True,
                        help="SegDup BED")
    parser.add_argument("--out-table", required=True,
                        help="Output master table TSV")
    parser.add_argument("--out-bed",   required=True,
                        help="Output population BED")
    parser.add_argument("--seq-type",  default="long",
                        help="Sequencing type label (default: long)")
    parser.add_argument("--promoter-window", type=int, default=2000,
                        help="bp upstream of TSS considered a promoter (default: 2000)")
    args = parser.parse_args()

    # ── 1. Load CxSV summary ──────────────────────────────────────────────────
    cxsv_loci = []
    with open(args.summary) as fh:
        raw_header = fh.readline()
        cols = raw_header.rstrip("\n").split("\t")
        try:
            ix = {c: i for i, c in enumerate(cols)}
            ci_id       = ix["ClusterID"]
            ci_chr      = ix["Chr"]
            ci_start    = ix["Start"]
            ci_end      = ix["End"]
            ci_size     = ix["Size"]
            ci_nbp      = ix["N_Breakpoints"]
            ci_types    = ix["SV_Types"]
            ci_inter    = ix.get("Has_Inter", None)
            ci_is_cxsv  = ix["IsCxSV"]
            # CxSV_Criteria column — present in classify_cxsv.py output
            ci_criteria = ix.get("CxSV_Criteria", None)
        except KeyError as e:
            sys.exit(f"ERROR: Expected column not found in cxsv_summary.tsv: {e}")

        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= ci_is_cxsv:
                continue
            if parts[ci_is_cxsv] != "YES":
                continue
            cxsv_loci.append({
                "id":        parts[ci_id],
                "chr":       parts[ci_chr],
                "start":     int(parts[ci_start]),
                "end":       int(parts[ci_end]),
                "size":      int(parts[ci_size]),
                "n_bp":      int(parts[ci_nbp]),
                "sv_types":  parts[ci_types],
                "has_inter": parts[ci_inter] if ci_inter is not None else "NA",
                "criteria":  parts[ci_criteria] if ci_criteria is not None else "",
            })

    if not cxsv_loci:
        sys.stderr.write(
            "WARNING: No CxSV loci found in cxsv_summary.tsv. "
            "Outputs will be header-only.\n"
        )

    # ── 2. Temporary query BED for bedtools calls ─────────────────────────────
    tmp_dir   = tempfile.mkdtemp()
    query_bed = os.path.join(tmp_dir, "query.bed")
    with open(query_bed, "w") as fh:
        for locus in cxsv_loci:
            fh.write(
                f"{locus['chr']}\t{locus['start']}\t{locus['end']}\t{locus['id']}\n"
            )

    # ── 3. Gene body overlap ──────────────────────────────────────────────────
    gene_overlap: dict = defaultdict(list)
    for line in bedtools_intersect(query_bed, args.genes, "-wa -wb"):
        parts = line.split("\t")
        if len(parts) >= 8:
            uid, gene = parts[3], parts[7]
            if gene not in gene_overlap[uid]:
                gene_overlap[uid].append(gene)

    # ── 4. Nearest gene + distance ────────────────────────────────────────────
    nearest_gene: dict = {}
    nearest_dist: dict = {}
    for line in bedtools_closest(query_bed, args.genes):
        parts = line.split("\t")
        if len(parts) >= 9:
            uid               = parts[3]
            nearest_gene[uid] = parts[7]
            nearest_dist[uid] = parts[-1]

    # ── 5. Exon coverage fraction ─────────────────────────────────────────────
    # Use coverage fraction rather than any-overlap so that a CxSV locus
    # spanning 50 kb is only called "Exonic" if >=10% of its length is
    # covered by exon sequence — not just because it touches one exon.
    exon_coverage: dict = coverage_fraction(query_bed, args.exons)

    # ── 6. Repeat annotation ──────────────────────────────────────────────────
    repeat_annot: dict = defaultdict(list)
    for line in bedtools_intersect(query_bed, args.repeats, "-wa -wb"):
        parts = line.split("\t")
        if len(parts) >= 8:
            uid, rep_class = parts[3], parts[7]
            if rep_class not in repeat_annot[uid]:
                repeat_annot[uid].append(rep_class)

    # ── 7. SegDup overlap fraction ────────────────────────────────────────────
    segdup_frac = coverage_fraction(query_bed, args.segdups)

    # ── 8. Parse VCF for breakpoint positions (BND/TRA-aware) ────────────────
    bp1: dict = {}
    bp2: dict = {}
    with open(args.vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, vid, _, alt_allele, _, _, info_str = parts[:8]
            uid = re.sub(r'(_mate|_dest)$', '', vid).split(".")[0]

            if uid not in bp1:
                bp1[uid] = f"{chrom}:{pos}"
            else:
                # BND ALT bracket notation
                m = re.search(r'[\[\]]([^\[\]]+):(\d+)[\[\]]', alt_allele)
                if m:
                    bp2[uid] = f"{m.group(1)}:{m.group(2)}"
                else:
                    # TRA: CHR2/END2 in INFO
                    info_d = {}
                    for f in info_str.split(";"):
                        if "=" in f:
                            k, v = f.split("=", 1)
                            info_d[k] = v
                    chr2 = info_d.get("CHR2", chrom)
                    end2 = info_d.get("END2", info_d.get("END", pos))
                    bp2[uid] = f"{chr2}:{end2}"

    # ── 9. Write master table ─────────────────────────────────────────────────
    os.makedirs(os.path.dirname(args.out_table) or ".", exist_ok=True)

    header_cols = [
        "Chr", "Start", "End", "Size",
        "Breakpoint_1", "Breakpoint_2",
        "SV_Types", "Cmplx_Class", "CxSV_Criteria",
        "Seq_Type", "Has_Inter_Chromosomal",
        "Genomic_Context",              # annotation-only; CxSVs never filtered here
        "Gene_Overlap", "Nearest_Gene", "Distance_To_Nearest_Gene",
        "Exon_Overlap", "Repeat_Annotation", "SegDup_Overlap"
    ]

    n_context: dict = defaultdict(int)   # tally for logging

    with open(args.out_table, "w") as out:
        out.write("\t".join(header_cols) + "\n")

        for locus in cxsv_loci:
            uid      = locus["id"]
            genes    = ",".join(gene_overlap[uid]) if gene_overlap[uid] else "NA"
            ng       = nearest_gene.get(uid, "NA")
            nd       = nearest_dist.get(uid, "NA")
            exon     = "yes" if exon_coverage.get(uid, 0.0) >= 0.10 else "no"
            rep      = ",".join(repeat_annot[uid]) if repeat_annot[uid] else "NA"
            segdup   = f"{segdup_frac.get(uid, 0.0):.4f}"
            criteria = locus.get("criteria", "")
            cmplx    = classify_complexity(locus["sv_types"], locus["n_bp"], criteria)
            context  = assign_genomic_context(
                uid, gene_overlap, exon_coverage, nearest_dist,
                repeat_annot, segdup_frac, args.promoter_window
            )
            b1 = bp1.get(uid, f"{locus['chr']}:{locus['start']}")
            b2 = bp2.get(uid, f"{locus['chr']}:{locus['end']}")

            n_context[context] += 1

            row = [
                locus["chr"], str(locus["start"]), str(locus["end"]),
                str(locus["size"]),
                b1, b2,
                locus["sv_types"], cmplx,
                criteria if criteria else "NONE",
                args.seq_type,
                locus.get("has_inter", "NA"),
                context,
                genes, ng, nd,
                exon, rep, segdup
            ]
            out.write("\t".join(row) + "\n")

    # ── 10. Write population BED ──────────────────────────────────────────────
    with open(args.out_bed, "w") as bed:
        for locus in cxsv_loci:
            uid   = locus["id"]
            exon  = "Y" if exon_coverage.get(uid, 0.0) >= 0.10 else "N"
            # Clamp coordinates: no negative starts, no inverted intervals
            start = max(0, locus["start"])
            end   = locus["end"]
            if end < start:
                start, end = end, start   # swap inverted coordinates
            if end <= start:
                end = start + 1           # ensure non-zero interval length
            bed.write(f"{locus['chr']}\t{start}\t{end}\t{exon}\n")

    # ── Summary to stderr ─────────────────────────────────────────────────────
    print(f"Master table: {args.out_table}  ({len(cxsv_loci)} CxSV loci)")
    print(f"Population BED: {args.out_bed}")
    print("\nGenomic context breakdown (annotation-only — no loci filtered):")
    for ctx, n in sorted(n_context.items(), key=lambda x: -x[1]):
        pct = 100 * n / len(cxsv_loci) if cxsv_loci else 0
        print(f"  {ctx:<30} {n:>5}  ({pct:.1f}%)")
    print(
        "\nNote: Intergenic/repeat/segdup loci are expected to be the majority.\n"
        "These are included in all downstream analyses."
    )


if __name__ == "__main__":
    main()
