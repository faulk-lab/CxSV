#!/usr/bin/env python3
"""Write one sample's junction table: every merged SV call as a breakpoint junction.

Reads the sample's merged_sv.vcf (SURVIVOR merge of the three callers) and its
CxSV clustering (cxsv/clusters.bed + cxsv/cxsv_summary.tsv), and writes
junctions/junctions.tsv with one row per junction (duplicate and reciprocal
breakend records collapsed, workflow/lib/cxsv_junctions.py):

  chrom1 pos1 chrom2 pos2 svtype svlen   the junction (shared-schema conventions)
  locus_id    CxSV locus/loci it belongs to ("Chr:Start-End" as in
              cxsv_master_table.tsv), or "." for SVs outside CxSV loci
  supp_vec    callers supporting it (cuteSV, SVIM, Sniffles2), OR-ed over merged records
  dosage      1 (heterozygous) or 2 (homozygous), from the caller genotypes
  n_records, record_ids   merged_sv.vcf records collapsed into this junction

The whole SV set is written, not only CxSV loci: the panel needs every
junction to count carriers of a patient locus's constituent SVs.

Usage: extract_junctions.py --vcf merged_sv.vcf --clusters clusters.bed
           --summary cxsv_summary.tsv --fai ref.fa.fai --collapse-distance 50 --out junctions.tsv
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import cxsv_junctions as cj  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--clusters", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--collapse-distance", type=int, required=True)
    ap.add_argument("--ins-len-ratio", type=float, default=0.5)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    order = cj.contig_order(a.fai)
    n_rec = sum(1 for l in open(a.vcf) if not l.startswith("#") and l.strip())
    js = cj.sample_junctions(a.vcf, order, a.collapse_distance, a.clusters, a.summary, a.ins_len_ratio)
    cj.write_table(a.out, js)
    in_loci = sum(1 for j in js if j.locus_ids)
    print(f"{a.vcf}: {n_rec} records -> {len(js)} junctions "
          f"({n_rec - sum(len(j.record_ids) for j in js)} unresolved records skipped; "
          f"{in_loci} junctions in CxSV loci)", file=sys.stderr)


if __name__ == "__main__":
    main()
