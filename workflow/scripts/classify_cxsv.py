#!/usr/bin/env python3
"""
classify_cxsv.py
================
Classify SV clusters as CxSV using all five project criteria.

CxSV Definition
---------------
A locus is a Complex Structural Variant (CxSV) if it satisfies
ANY of the following:

  C1 — ≥3 unique breakpoints within the 50 kb cluster window
       (BND/TRA _mate/_dest records are deduplicated back to the
       canonical variant ID before counting)

  C2 — Nested SV: at least one variant's [start, end] interval is
       wholly contained within a different variant's interval in the
       same cluster (e.g. a DEL sitting entirely inside an INV)

  C3 — Overlapping SV types: ≥2 distinct SV types whose coordinate
       intervals physically share at least 1 bp (not just co-clustered
       within 50 kb, but genuinely overlapping on the genome)

  C4 — Copy-number shift + orientation change: the cluster contains
       BOTH a CN-altering SV (DEL or DUP) AND an orientation-altering
       SV (INV, BND, or TRA) — hallmark of FoSTeS/MMBIR events

  C5 — Chromoanagenesis signature: ≥3 unique SVs packed within a
       10–50 kb window, consistent with localised catastrophic repair
       (FoSTeS/MMBIR or chromothripsis)

Important: genomic context (intergenic, intronic, repeat, segdup) is
NOT used as a criterion. CxSVs are expected to be enriched in the
'grey zone' of the genome and must not be filtered by annotation.

Input
-----
clusters.bed produced by bedtools cluster on merged_sv.bed.
Expected columns (tab-separated, no header):
  CHROM  START  END  SVTYPE  ID  IS_INTER  CLUSTER_ID

Output TSV (tab-separated, with header)
---------------------------------------
  ClusterID  Chr  Start  End  Size
  N_Breakpoints  N_SVtypes  SV_Types
  Has_Inter  IsCxSV  CxSV_Criteria
"""

import argparse
import sys
from collections import defaultdict


# ─── Variant record ───────────────────────────────────────────────────────────

class Variant:
    __slots__ = ("chrom", "start", "end", "svtype", "vid", "is_inter", "cluster")

    def __init__(self, chrom, start, end, svtype, vid, is_inter, cluster):
        self.chrom    = chrom
        self.start    = int(start)
        self.end      = int(end)
        self.svtype   = svtype.upper()
        self.vid      = vid
        self.is_inter = is_inter
        self.cluster  = cluster

    @property
    def canonical_id(self):
        """Strip _mate/_dest suffixes added by vcf_to_bed.py."""
        cid = self.vid
        if cid.endswith("_mate"):
            cid = cid[:-5]
        elif cid.endswith("_dest"):
            cid = cid[:-5]
        return cid

    @property
    def is_mate_record(self):
        return self.vid.endswith("_mate") or self.vid.endswith("_dest")

    @property
    def is_span_like(self):
        """True for DEL/DUP/INV/INS — records that have a meaningful span."""
        return (self.svtype not in ("BND", "TRA")
                and not self.is_mate_record
                and self.end > self.start)

    @property
    def is_cn(self):
        return self.svtype in ("DEL", "DUP", "DUP:TANDEM")

    @property
    def is_orient(self):
        return self.svtype in ("INV", "BND", "TRA")


# ─── Five criteria ────────────────────────────────────────────────────────────

def criterion1(variants: list, min_bp: int) -> bool:
    """C1: ≥min_bp unique canonical breakpoints in the cluster."""
    return len({v.canonical_id for v in variants}) >= min_bp


def criterion2(variants: list) -> bool:
    """
    C2: Nested SV — any span-like variant is wholly contained
    inside a different span-like variant (different canonical ID).
    """
    spans = [v for v in variants if v.is_span_like]
    for i, inner in enumerate(spans):
        for j, outer in enumerate(spans):
            if i == j or inner.canonical_id == outer.canonical_id:
                continue
            if outer.start <= inner.start and inner.end <= outer.end:
                return True
    return False


def criterion3(variants: list) -> bool:
    """
    C3: Overlapping SV types — ≥2 distinct types whose intervals
    physically share at least 1 bp.
    """
    spans = [v for v in variants if v.is_span_like]
    n = len(spans)
    for i in range(n):
        for j in range(i + 1, n):
            a, b = spans[i], spans[j]
            if (a.svtype != b.svtype
                    and a.canonical_id != b.canonical_id
                    and a.start < b.end
                    and b.start < a.end):
                return True
    return False


def criterion4(variants: list) -> bool:
    """C4: Copy-number shift (DEL/DUP) AND orientation change (INV/BND/TRA)."""
    return any(v.is_cn for v in variants) and any(v.is_orient for v in variants)


def criterion5(variants: list, min_window: int, max_window: int,
               min_svs: int) -> bool:
    """
    C5: ≥min_svs unique SVs within a [min_window, max_window] bp span.
    Uses a sliding window over position-sorted unique canonical variants.
    """
    # One representative record per canonical ID, sorted by start position
    seen: dict[str, Variant] = {}
    for v in variants:
        cid = v.canonical_id
        if cid not in seen or v.start < seen[cid].start:
            seen[cid] = v
    uniq = sorted(seen.values(), key=lambda v: v.start)

    n = len(uniq)
    if n < min_svs:
        return False

    for i in range(n):
        span_start = uniq[i].start
        span_end   = uniq[i].end
        count      = 1
        for j in range(i + 1, n):
            span_end  = max(span_end, uniq[j].end)
            span_size = span_end - span_start
            if span_size > max_window:
                break
            count += 1
            if count >= min_svs and span_size >= min_window:
                return True
    return False


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(
        description="Classify SV clusters as CxSV using all 5 project criteria."
    )
    p.add_argument("--bed",               required=True,  help="clusters.bed from bedtools")
    p.add_argument("--output",            required=True,  help="Output TSV")
    p.add_argument("--min-bp",            type=int, default=3)
    p.add_argument("--max-window",        type=int, default=50000)
    p.add_argument("--chromoanag-min",    type=int, default=10000)
    p.add_argument("--chromoanag-max",    type=int, default=50000)
    p.add_argument("--chromoanag-min-sv", type=int, default=3)
    p.add_argument("--max-cluster-size",  type=int, default=500000,
                   help=(
                       "Maximum total cluster span in bp to be considered a CxSV. "
                       "Clusters larger than this are the result of bedtools chaining "
                       "(SV A -> B -> C -> ... spanning Mb) rather than true localised "
                       "complex rearrangements. Default: 500000 (500 kb). "
                       "Set to 0 to disable this filter."
                   ))
    args = p.parse_args()

    # ── Parse clusters.bed ────────────────────────────────────────────────────
    clusters: dict[str, list[Variant]] = defaultdict(list)
    with open(args.bed) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 7:
                sys.stderr.write(f"WARN: short line skipped: {line!r}\n")
                continue
            chrom, start, end, svtype, vid, is_inter, cluster_id = parts[:7]
            clusters[cluster_id].append(
                Variant(chrom, start, end, svtype, vid, is_inter, cluster_id)
            )

    sys.stderr.write(f"classify_cxsv: {len(clusters)} clusters loaded\n")

    # ── Header ────────────────────────────────────────────────────────────────
    header = [
        "ClusterID", "Chr", "Start", "End", "Size",
        "N_Breakpoints", "N_SVtypes", "SV_Types",
        "Has_Inter", "IsCxSV", "CxSV_Criteria"
    ]

    rows = []
    n_cxsv = 0

    for cid, variants in clusters.items():
        chrom     = variants[0].chrom
        start     = min(v.start for v in variants)
        end       = max(v.end   for v in variants)
        size      = end - start
        unique_ids = {v.canonical_id for v in variants}
        n_bp      = len(unique_ids)
        type_set  = {v.svtype for v in variants}
        sv_types  = ",".join(sorted(type_set))
        n_types   = len(type_set)
        has_inter = "YES" if any(v.is_inter == "inter" for v in variants) else "NO"

        # ── Max cluster size filter ───────────────────────────────────────────
        # Clusters spanning more than max_cluster_size bp are bedtools chaining
        # artifacts (SV A->B->C->... connected within 50kb windows but spanning
        # hundreds of Mb in total). These are NOT true CxSV loci.
        # We record them in the output as IsCxSV=NO so they are traceable.
        if args.max_cluster_size > 0 and size > args.max_cluster_size:
            rows.append([
                cid, chrom, str(start), str(end), str(size),
                str(n_bp), str(n_types), sv_types,
                has_inter, "NO", "FILTERED_oversized_cluster"
            ])
            continue

        # ── Apply five criteria ───────────────────────────────────────────────
        fired = []
        if criterion1(variants, args.min_bp):
            fired.append("C1_min_breakpoints")
        if criterion2(variants):
            fired.append("C2_nested_SV")
        if criterion3(variants):
            fired.append("C3_overlapping_types")
        if criterion4(variants):
            fired.append("C4_CN_plus_orient")
        if criterion5(variants, args.chromoanag_min, args.chromoanag_max,
                      args.chromoanag_min_sv):
            fired.append("C5_chromoanagenesis_cluster")

        is_cxsv  = "YES" if fired else "NO"
        criteria = ",".join(fired) if fired else "NONE"
        if is_cxsv == "YES":
            n_cxsv += 1

        rows.append([
            cid, chrom, str(start), str(end), str(size),
            str(n_bp), str(n_types), sv_types,
            has_inter, is_cxsv, criteria
        ])

    # ── Sort by chromosome then position ─────────────────────────────────────
    def chr_sort_key(row):
        c = row[1].lstrip("chr")
        return (0, int(c)) if c.isdigit() else (1, c)

    rows.sort(key=lambda r: (chr_sort_key(r), int(r[2])))

    # ── Write ─────────────────────────────────────────────────────────────────
    with open(args.output, "w") as out:
        out.write("\t".join(header) + "\n")
        for row in rows:
            out.write("\t".join(row) + "\n")

    sys.stderr.write(
        f"classify_cxsv: {n_cxsv} CxSV / {len(clusters)} clusters\n"
    )
    for criterion_name in [
        "C1_min_breakpoints", "C2_nested_SV", "C3_overlapping_types",
        "C4_CN_plus_orient",  "C5_chromoanagenesis_cluster"
    ]:
        n = sum(1 for r in rows if criterion_name in r[10])
        sys.stderr.write(f"  {criterion_name}: {n}\n")


if __name__ == "__main__":
    main()
