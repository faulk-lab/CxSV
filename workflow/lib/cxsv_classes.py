"""CxSV locus structure: distinct junctions, breakends, interleaving, criteria C1-C5, class.

One SV cluster (bedtools cluster on a sample's merged SV calls, 50 kb) becomes
a *locus*. Its structure is computed on distinct junctions, not on VCF records
(roadmap item 5; positive control E12):

1. **Collapse duplicate and reciprocal records.** Records with the same svtype
   and contigs whose two ends both lie within `collapse_distance` bp are one
   junction (`cxsv_junctions.collapse`, the same rule as the junction tables:
   Sniffles2 writes each inversion breakpoint as two reciprocal records, a BND
   and its mate are two records). A simple balanced inversion is then ONE
   junction with two breakends, and no longer "nested" in itself.
2. **Drop composite records.** SVIM reports two offset inversion junctions
   (e.g. the two junctions of a delINVdel) as one extra inversion placed
   between them. A single-caller junction whose two ends each lie at the
   midpoint of the corresponding ends of two other, overlapping junctions of
   the same type (within max(collapse_distance, composite_fraction x their
   distance)) adds
   no breakend of its own and is not counted. The record stays in the VCF and
   in the junction tables.
3. **Breakends** are the distinct break positions of the remaining junctions
   (both ends, partner ends on other chromosomes included), merged within
   `collapse_distance` bp: two segment ends joined at one DNA break count once.
   This is the breakpoint count X of the PI's criteria (criterion 2):
   intermediate (3-5) vs extreme (> 5, i.e. >= extreme_min_breakends).
4. **Interleaving** (criterion 1, ShatterSeek's defining relation): two
   intrachromosomal junctions interleave when their intervals overlap without
   one containing the other (a1 < b1 < a2 < b2, every gap > collapse_distance).
   Max_Interleaved is the size of the largest connected component.

Criteria (any one makes the locus a CxSV; genomic context is never used):

  C1  >= min_breakpoints distinct junctions, or >= 2 junctions with
      >= min_breakpoints distinct breakends of which at least one changes
      orientation (INV or BND): an inversion or translocation whose junctions
      do not share their breakpoints (sequence lost or gained at the breaks,
      e.g. delINVdel, an unbalanced reciprocal translocation). Two nearby
      deletions are not a CxSV by C1.
  C2  nested: one span-like junction lies inside another
  C3  two span-like junctions of different types overlap by >= 1 bp
  C4  a copy-number junction (DEL/DUP) and an orientation junction (INV/BND)
  C5  >= chromoanag_min_svs junctions within a 10-50 kb window

Class (Cmplx_Class), first match wins; `Chromothripsis_like` no longer comes
from C5 alone (it labelled compact 2-4 junction events):

  Chromothripsis_like         extreme tier AND interleaving (component >=
                              min_interleaved, or >= alt_min_intra interleaved
                              plus >= alt_min_inter interchromosomal junctions)
  DEL+DUP+INV, DEL+INV, DUP+INV, INS+DEL   called junction types
  BND_cluster                 >= 2 distinct BND (translocation) junctions
  Multi_breakpoint_cluster    extreme tier without interleaving
  Nested_SV / Overlapping_types / CN_plus_orientation_change   C2 / C3 / C4
  Complex_inversion           >= 2 distinct inversion junctions, nothing else
  Complex_other               anything else that passed C1-C5

Types are the called junction types: a deletion created by an inversion
junction (delINVdel) has no DEL record of its own and is not a "DEL" here.
Copy number is assessed separately from read depth in scoring (criterion 5).

This module decides which loci exist and how they are labelled, so it is part
of the calling fingerprint (cxsv_fingerprint.CALLING_CODE).
"""
from __future__ import annotations

import cxsv_junctions as cj

CN_TYPES = ("DEL", "DUP")
ORIENT_TYPES = ("INV", "BND")
SPAN_TYPES = ("DEL", "DUP", "INV")          # intervals that can nest / overlap / interleave

CRITERIA = ("C1_min_breakpoints", "C2_nested_SV", "C3_overlapping_types",
            "C4_CN_plus_orient", "C5_chromoanagenesis_cluster")
CLASSES = ("Chromothripsis_like", "DEL+DUP+INV", "DEL+INV", "DUP+INV", "INS+DEL", "BND_cluster",
           "Multi_breakpoint_cluster", "Nested_SV", "Overlapping_types", "CN_plus_orientation_change",
           "Complex_inversion", "Complex_other")
TIERS = ("simple", "intermediate", "extreme")

DEFAULTS = {
    "collapse_distance": 300,
    "ins_len_ratio": 0.5,
    "composite_fraction": 0.25,
    "min_breakpoints": 3,
    "chromoanag_min_window": 10000,
    "chromoanag_max_window": 50000,
    "chromoanag_min_svs": 3,
    "extreme_min_breakends": 6,
    "min_interleaved": 3,
    "interleave_alt_min_intra": 3,
    "interleave_alt_min_inter": 4,
}


def n_callers(j):
    return j.supp_vec.count("1") if j.supp_vec not in (".", "") else 0


def is_span(j):
    return j.svtype in SPAN_TYPES and j.chrom1 == j.chrom2 and j.pos2 > j.pos1


# ─── 1-2. distinct junctions ─────────────────────────────────────────────────
def _near_mid(x, a, b, tol, frac):
    return abs(x - (a + b) / 2.0) <= max(tol, frac * abs(a - b))


def composite_of(x, others, tol, frac):
    """(a, b) when single-caller junction x sits at the midpoint of two other,
    overlapping same-type junctions a and b at both ends, else None."""
    if n_callers(x) != 1 or x.chrom1 != x.chrom2:
        return None
    cand = [o for o in others if o is not x and o.svtype == x.svtype
            and o.chrom1 == x.chrom1 and o.chrom2 == x.chrom2]
    for i, a in enumerate(cand):
        for b in cand[i + 1:]:
            if (abs(a.pos1 - b.pos1) > tol or abs(a.pos2 - b.pos2) > tol) \
                    and a.pos1 < b.pos2 and b.pos1 < a.pos2 \
                    and _near_mid(x.pos1, a.pos1, b.pos1, tol, frac) \
                    and _near_mid(x.pos2, a.pos2, b.pos2, tol, frac):
                return a, b
    return None


def distinct_junctions(junctions, order, p):
    """(kept junctions, composite junctions) for one locus."""
    merged = cj.collapse(junctions, p["collapse_distance"], order, p["ins_len_ratio"])
    composite = [j for j in merged
                 if composite_of(j, merged, p["collapse_distance"], p["composite_fraction"])]
    kept = [j for j in merged if all(j is not c for c in composite)]
    return kept, composite


# ─── 3. breakends ────────────────────────────────────────────────────────────
def breakends(junctions, tol):
    """Distinct break positions [(chrom, pos)], ends within tol merged (leader-based)."""
    ends = {}
    for j in junctions:
        ends.setdefault(j.chrom1, []).append(j.pos1)
        if not (j.svtype == "INS" and j.chrom2 == j.chrom1 and j.pos2 == j.pos1):
            ends.setdefault(j.chrom2, []).append(j.pos2)
    out = []
    for chrom in sorted(ends):
        lead = None
        for p in sorted(ends[chrom]):
            if lead is None or p - lead > tol:
                lead = p
                out.append((chrom, p))
    return out


def tier(n_breakends, extreme_min):
    if n_breakends >= extreme_min:
        return "extreme"
    return "intermediate" if n_breakends >= 3 else "simple"


# ─── 4. interleaving ─────────────────────────────────────────────────────────
def interleaved(a, b, tol):
    """Intervals overlap, neither contains the other, every gap > tol."""
    if a.chrom1 != b.chrom1:
        return False
    if (b.pos1, b.pos2) < (a.pos1, a.pos2):
        a, b = b, a
    return (b.pos1 - a.pos1 > tol and a.pos2 - b.pos1 > tol and b.pos2 - a.pos2 > tol)


def interleaved_components(junctions, tol):
    """Connected components (lists of junctions) of the interleaving graph of
    intrachromosomal span junctions, largest first."""
    spans = [j for j in junctions if is_span(j)]
    parent = list(range(len(spans)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(len(spans)):
        for k in range(i + 1, len(spans)):
            if interleaved(spans[i], spans[k], tol):
                parent[find(i)] = find(k)
    comp = {}
    for i in range(len(spans)):
        comp.setdefault(find(i), []).append(spans[i])
    return sorted(comp.values(), key=lambda c: (-len(c), c[0].pos1))


# ─── criteria ────────────────────────────────────────────────────────────────
def c1(js, n_bkends, p):
    if len(js) >= p["min_breakpoints"]:
        return True
    return (len(js) >= 2 and n_bkends >= p["min_breakpoints"]
            and any(j.svtype in ORIENT_TYPES for j in js))


def c2(js):
    spans = [j for j in js if is_span(j)]
    for i, inner in enumerate(spans):
        for k, outer in enumerate(spans):
            if i != k and outer.chrom1 == inner.chrom1 \
                    and outer.pos1 <= inner.pos1 and inner.pos2 <= outer.pos2:
                return True
    return False


def c3(js):
    spans = [j for j in js if is_span(j)]
    for i, a in enumerate(spans):
        for b in spans[i + 1:]:
            if a.svtype != b.svtype and a.chrom1 == b.chrom1 and a.pos1 < b.pos2 and b.pos1 < a.pos2:
                return True
    return False


def c4(js):
    return any(j.svtype in CN_TYPES for j in js) and any(j.svtype in ORIENT_TYPES for j in js)


def c5(js, chrom, p):
    """>= chromoanag_min_svs junctions within a [min_window, max_window] span (positions on
    the locus chromosome; a translocation counts at its end on this chromosome)."""
    iv = []
    for j in js:
        if j.chrom1 == j.chrom2 == chrom:
            iv.append((j.pos1, max(j.pos1, j.pos2)))
        elif j.chrom1 == chrom:
            iv.append((j.pos1, j.pos1))
        elif j.chrom2 == chrom:
            iv.append((j.pos2, j.pos2))
    iv.sort()
    n, k = len(iv), p["chromoanag_min_svs"]
    if n < k:
        return False
    for i in range(n):
        start, end, count = iv[i][0], iv[i][1], 1
        for t in range(i + 1, n):
            end = max(end, iv[t][1])
            if end - start > p["chromoanag_max_window"]:
                break
            count += 1
            if count >= k and end - start >= p["chromoanag_min_window"]:
                return True
    return False


# ─── class ───────────────────────────────────────────────────────────────────
def interleave_ok(max_interleaved, n_inter, p):
    return (max_interleaved >= p["min_interleaved"]
            or (max_interleaved >= p["interleave_alt_min_intra"]
                and n_inter >= p["interleave_alt_min_inter"]))


def complexity_class(types, n_breakends, criteria, max_interleaved, n_inter, n_bnd, n_inv, p=None):
    """Cmplx_Class of a CxSV locus (see module docstring for the precedence)."""
    p = {**DEFAULTS, **(p or {})}
    types = set(types)
    crit = set(criteria)
    extreme = n_breakends >= p["extreme_min_breakends"]
    if extreme and interleave_ok(max_interleaved, n_inter, p):
        return "Chromothripsis_like"
    if {"DEL", "INV", "DUP"} <= types:
        return "DEL+DUP+INV"
    if {"DEL", "INV"} <= types:
        return "DEL+INV"
    if {"DUP", "INV"} <= types:
        return "DUP+INV"
    if {"INS", "DEL"} <= types:
        return "INS+DEL"
    if n_bnd >= 2:
        return "BND_cluster"
    if extreme:
        return "Multi_breakpoint_cluster"
    if "C2_nested_SV" in crit:
        return "Nested_SV"
    if "C3_overlapping_types" in crit:
        return "Overlapping_types"
    if "C4_CN_plus_orient" in crit:
        return "CN_plus_orientation_change"
    if types == {"INV"} and n_inv >= 2:
        return "Complex_inversion"
    return "Complex_other"


def analyse_locus(junctions, chrom, order, params=None):
    """Structure, criteria and class of one locus from its records' junctions."""
    p = {**DEFAULTS, **(params or {})}
    tol = p["collapse_distance"]
    kept, composite = distinct_junctions(junctions, order, p)
    bk = breakends(kept, tol)
    comps = interleaved_components(kept, tol)
    max_il = len(comps[0]) if comps else 0
    n_inter = sum(1 for j in kept if j.chrom1 != j.chrom2)
    fired = []
    if c1(kept, len(bk), p):
        fired.append(CRITERIA[0])
    if c2(kept):
        fired.append(CRITERIA[1])
    if c3(kept):
        fired.append(CRITERIA[2])
    if c4(kept):
        fired.append(CRITERIA[3])
    if c5(kept, chrom, p):
        fired.append(CRITERIA[4])
    types = sorted({j.svtype for j in kept})
    n_bnd = sum(1 for j in kept if j.svtype == "BND")
    n_inv = sum(1 for j in kept if j.svtype == "INV")
    cls = complexity_class(types, len(bk), fired, max_il, n_inter, n_bnd, n_inv, p) if fired else "NA"
    return {
        "junctions": kept, "composite": composite, "breakends": bk,
        "n_records": len(junctions), "n_junctions": len(kept), "n_composite": len(composite),
        "n_breakends": len(bk), "max_interleaved": max_il, "n_inter": n_inter,
        "tier": tier(len(bk), p["extreme_min_breakends"]),
        "types": types, "criteria": fired, "is_cxsv": bool(fired), "class": cls,
    }
