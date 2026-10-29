#!/usr/bin/env python3
"""
classify_cxsv.py
================
Classify one sample's SV clusters as CxSV loci (criteria C1-C5) and label them.

The structure of each cluster is computed on DISTINCT JUNCTIONS, not on VCF
records (workflow/lib/cxsv_classes.py, roadmap item 5): duplicate and
reciprocal records (same svtype and contigs, both ends within
--collapse-distance; e.g. Sniffles2's two records per inversion breakpoint)
are one junction, SVIM's composite inversions between two offset inversion
junctions are not counted, and breakends are the distinct break positions.
A simple balanced inversion is therefore one junction with two breakends and
not a CxSV (positive control E12).

  C1  >= --min-bp junctions, or >= 2 junctions with >= --min-bp breakends and
      an orientation change (INV/BND)
  C2  nested span junctions          C3  overlapping span junctions of 2 types
  C4  DEL/DUP plus INV/BND           C5  >= --chromoanag-min-sv junctions in 10-50 kb

Genomic context (intergenic, intronic, repeat, segdup) is never a criterion.

Inputs
------
--bed   clusters.bed (bedtools cluster on merged_sv.bed):
        CHROM START END SVTYPE ID IS_INTER CLUSTER_ID  (cluster membership; locus span)
--vcf   merged_sv.vcf (the records' junctions and SUPP_VEC, cxsv_junctions.parse_vcf)
--fai   reference .fai (contig order of junction ends)

Output TSV (tab-separated, with header)
---------------------------------------
  ClusterID Chr Start End Size
  N_Breakpoints     distinct breakends (break positions; criterion-2 count X)
  N_SVtypes SV_Types   called types of the distinct junctions (TRA written BND)
  Has_Inter IsCxSV CxSV_Criteria
  N_Records         merged-VCF records in the cluster
  N_Junctions       distinct junctions after collapsing duplicates/reciprocals
  N_Composite       single-caller composite records not counted
  Max_Interleaved   largest interleaved component (criterion 1)
  N_Inter_Junctions translocation junctions (ends on two chromosomes)
  Complexity_Tier   simple (<3 breakends) / intermediate (3-5) / extreme (>5)
  Cmplx_Class       class label (cxsv_classes.complexity_class); NA if not a CxSV
"""

import argparse
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import cxsv_classes as cc     # noqa: E402
import cxsv_junctions as cj   # noqa: E402

HEADER = ["ClusterID", "Chr", "Start", "End", "Size", "N_Breakpoints", "N_SVtypes", "SV_Types",
          "Has_Inter", "IsCxSV", "CxSV_Criteria", "N_Records", "N_Junctions", "N_Composite",
          "Max_Interleaved", "N_Inter_Junctions", "Complexity_Tier", "Cmplx_Class"]


def read_clusters(path):
    """{cluster: {"chrom", "start", "end", "ids": set of canonical record ids, "inter"}}."""
    clusters = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                if line.strip():
                    sys.stderr.write(f"WARN: short line skipped: {line!r}\n")
                continue
            chrom, start, end, _svtype, vid, is_inter, cid = parts[:7]
            c = clusters.setdefault(cid, {"chrom": chrom, "start": int(start), "end": int(end),
                                          "ids": set(), "inter": False, "order": len(clusters)})
            c["start"] = min(c["start"], int(start))
            c["end"] = max(c["end"], int(end))
            c["ids"].add(cj.canonical_id(vid))
            c["inter"] = c["inter"] or is_inter == "inter"
    return clusters


def main():
    p = argparse.ArgumentParser(description="Classify SV clusters as CxSV loci (C1-C5) on distinct junctions.")
    p.add_argument("--bed", required=True, help="clusters.bed from bedtools cluster")
    p.add_argument("--vcf", required=True, help="merged_sv.vcf (records' junctions, SUPP_VEC)")
    p.add_argument("--fai", required=True, help="reference .fai (contig order)")
    p.add_argument("--output", required=True)
    p.add_argument("--collapse-distance", type=int, default=cc.DEFAULTS["collapse_distance"])
    p.add_argument("--ins-len-ratio", type=float, default=cc.DEFAULTS["ins_len_ratio"])
    p.add_argument("--composite-fraction", type=float, default=cc.DEFAULTS["composite_fraction"])
    p.add_argument("--min-bp", type=int, default=3)
    p.add_argument("--max-window", type=int, default=50000, help="(clustering distance; informational)")
    p.add_argument("--chromoanag-min", type=int, default=10000)
    p.add_argument("--chromoanag-max", type=int, default=50000)
    p.add_argument("--chromoanag-min-sv", type=int, default=3)
    p.add_argument("--extreme-min-breakends", type=int, default=cc.DEFAULTS["extreme_min_breakends"])
    p.add_argument("--min-interleaved", type=int, default=cc.DEFAULTS["min_interleaved"])
    p.add_argument("--interleave-alt-min-intra", type=int, default=cc.DEFAULTS["interleave_alt_min_intra"])
    p.add_argument("--interleave-alt-min-inter", type=int, default=cc.DEFAULTS["interleave_alt_min_inter"])
    p.add_argument("--max-cluster-size", type=int, default=500000,
                   help=("Clusters spanning more than this many bp come from bedtools chaining "
                         "(SV A -> B -> C ... over Mb), not from one localised rearrangement; "
                         "they are written with IsCxSV=NO. 0 disables the filter."))
    a = p.parse_args()
    params = {"collapse_distance": a.collapse_distance, "ins_len_ratio": a.ins_len_ratio,
              "composite_fraction": a.composite_fraction, "min_breakpoints": a.min_bp,
              "chromoanag_min_window": a.chromoanag_min, "chromoanag_max_window": a.chromoanag_max,
              "chromoanag_min_svs": a.chromoanag_min_sv, "extreme_min_breakends": a.extreme_min_breakends,
              "min_interleaved": a.min_interleaved, "interleave_alt_min_intra": a.interleave_alt_min_intra,
              "interleave_alt_min_inter": a.interleave_alt_min_inter}

    order = cj.contig_order(a.fai)
    by_id = {}
    for j in cj.parse_vcf(a.vcf, order):
        by_id.setdefault(j.record_ids[0], j)
    clusters = read_clusters(a.bed)
    sys.stderr.write(f"classify_cxsv: {len(clusters)} clusters, {len(by_id)} resolvable records\n")

    rows, n_cxsv = [], 0
    crit_counts, class_counts = defaultdict(int), defaultdict(int)
    for cid, c in sorted(clusters.items(), key=lambda kv: kv[1]["order"]):
        size = c["end"] - c["start"]
        js = [by_id[i] for i in sorted(c["ids"]) if i in by_id]
        res = cc.analyse_locus(js, c["chrom"], order, params)
        has_inter = "YES" if (c["inter"] or res["n_inter"]) else "NO"
        types = res["types"]
        base = [cid, c["chrom"], str(c["start"]), str(c["end"]), str(size), str(res["n_breakends"]),
                str(len(types)), ",".join(types) if types else "NA", has_inter]
        extra = [str(len(c["ids"])), str(res["n_junctions"]), str(res["n_composite"]),
                 str(res["max_interleaved"]), str(res["n_inter"]), res["tier"]]
        if a.max_cluster_size > 0 and size > a.max_cluster_size:
            rows.append(base + ["NO", "FILTERED_oversized_cluster"] + extra + ["NA"])
            continue
        if res["is_cxsv"]:
            n_cxsv += 1
            class_counts[res["class"]] += 1
            for k in res["criteria"]:
                crit_counts[k] += 1
        rows.append(base + ["YES" if res["is_cxsv"] else "NO",
                            ",".join(res["criteria"]) if res["criteria"] else "NONE"]
                    + extra + [res["class"]])

    def chr_sort_key(row):
        s = row[1][3:] if row[1].startswith("chr") else row[1]
        return (0, int(s), "") if s.isdigit() else (1, 0, s)

    rows.sort(key=lambda r: (chr_sort_key(r), int(r[2]), int(r[3]), r[0]))
    with open(a.output, "w") as out:
        out.write("\t".join(HEADER) + "\n")
        for r in rows:
            out.write("\t".join(r) + "\n")

    sys.stderr.write(f"classify_cxsv: {n_cxsv} CxSV / {len(clusters)} clusters\n")
    for k in cc.CRITERIA:
        sys.stderr.write(f"  {k}: {crit_counts[k]}\n")
    for k in cc.CLASSES:
        if class_counts[k]:
            sys.stderr.write(f"  class {k}: {class_counts[k]}\n")


if __name__ == "__main__":
    main()
