"""Structural reconstruction of one CxSV locus from a run's own outputs (clinician report).

For a scored locus (normally each Tier 1 locus of a patient run) this module derives

* the **reference segments**: the intervals between the locus's junction breakpoints
  (junctions/junctions.tsv, distinct junctions as in scoring), with a flank on each side;
* the **novel adjacencies (junctions)**, each with its join orientation and support,
  read from the callers' own records (callers/*_filtered.vcf: BND ALT brackets, symbolic
  DEL/DUP/INV, cuteSV INV STRAND) — the merged call set keeps one record per merged group and
  can lose an orientation (SURVIVOR types a BND pair as INV);
* the **copy number** of each segment, from the locus read depth (scoring/locus_depth.tsv,
  the same per-locus depth the scorer uses) or, for segments the locus depth does not
  measure, the run's copy-number segments (cn/segments.tsv);
* the **derivative**: a walk through the adjacency graph (reference adjacencies plus the
  novel junctions) from a flank to a flank that visits every segment as many times as its
  copy number on the rearranged haplotype and uses every junction;
* **gene names** overlapping each segment (the bundle's genes BED).

Nothing is read from alignments; the reconstruction works offline inside a bundle.

Rules (README "structural reconstruction" and docs/SCORING.md section 6f):

1. Heterozygous model: the derivative haplotype carries ``state - (ploidy - 1)`` copies of
   a segment (CN 1 of 2 = deleted, CN 3 = two copies). Segments without a usable depth
   estimate have unknown copy number and may be visited 0-2 times at a small cost.
2. The walk with the lowest cost is drawn. Costs: 10 per junction not on the derivative,
   10 per open end (the walk had to stop: a junction is probably missing), 5 per copy that
   disagrees with read depth, 1 per copy away from one for a segment without depth.
   Ties are broken by the read support of the joins used.
3. **Ambiguous** when the walk is not unique (several orders use the same joins equally
   well), a junction is left over, an end stays open, read depth disagrees with the path,
   or the search limit is hit. The panel says so, the alternatives are listed, and nothing
   is drawn that the junctions and depth do not show.
4. Caller type labels are not trusted over join orientation: a junction called DEL by one
   caller and INV by another is resolved by the orientation of the records and by read
   depth (a loss supports a deletion; an inversion keeps the copy number).

Not fingerprinted: reporting only (cxsv_fingerprint.CALLING_CODE does not list it).
"""
from __future__ import annotations

import csv
import gzip
import textwrap
import os
import re
from collections import defaultdict
from dataclasses import dataclass, field

import cxsv_junctions as cj
import cxsv_classes as cc
import cxsv_scoring as sc

CALLER_FILES = {"Sniffles2": "sniffles_filtered.vcf", "cuteSV": "cutesv_filtered.vcf", "SVIM": "svim_filtered.vcf"}
READ_KEYS = {"Sniffles2": ("SUPPORT",), "cuteSV": ("RE",), "SVIM": ("SUPPORT",)}
JOIN_NAMES = {("+", "-"): "deletion-type", ("-", "+"): "duplication-type",
              ("+", "+"): "inversion-type (head-to-head)", ("-", "-"): "inversion-type (tail-to-tail)"}

JOIN_KIND = {("+", "-"): "del", ("-", "+"): "dup", ("+", "+"): "inv", ("-", "-"): "inv"}
COST = {"unused_junction": 10, "open_end": 5, "cn_mismatch": 5, "unknown_cn": 1}
CN_SLACK = 0.25                    # depth noise allowed before a copy-number cost accrues
CN_FLAG = 0.6                      # |depth CN - path CN| above this is reported (scoring min_cn_delta)
SEARCH_LIMIT = 400000
MAX_SOLUTIONS = 200
CN_BOUNDARY_SLACK = 20000          # a CN segment boundary this close to a junction is the same break
BREAKPOINT_TOL = 500               # junction ends this close are one breakpoint (caller-agreement tolerance)
FLANK_BP = 20000                   # nominal flank length written to the segment table (display shortens)

_BND = (
    (re.compile(r"^[A-Za-z.*]+\[([^:\[\]]+):(\d+)\[$"), "+", "-"),   # t[p[ : piece ending at pos, then piece starting at p
    (re.compile(r"^[A-Za-z.*]+\]([^:\[\]]+):(\d+)\]$"), "+", "+"),   # t]p] : piece ending at pos, then rc of piece ending at p
    (re.compile(r"^\]([^:\[\]]+):(\d+)\][A-Za-z.*]+$"), "-", "+"),   # ]p]t : piece ending at p, then piece starting at pos
    (re.compile(r"^\[([^:\[\]]+):(\d+)\[[A-Za-z.*]+$"), "-", "-"),   # [p[t : rc of piece starting at p, then piece starting at pos
)


# ─── data model ──────────────────────────────────────────────────────────────
@dataclass
class Record:
    """One caller record near the locus: its two ends with their join sides.

    side '+' = the sequence ending at the break is kept (joined at its 3' end);
    side '-' = the sequence starting at the break. A symbolic INV without a strand stands
    for both of its junctions (orients has two entries)."""
    caller: str
    rid: str
    svtype: str
    chrom1: str
    pos1: int
    chrom2: str
    pos2: int
    orients: list            # [(side at end 1, side at end 2)]
    reads: int | None
    precise: bool            # orientation stated by the record (BND / DEL / DUP / stranded INV)


@dataclass
class Group:
    """A junction: one pair of breakpoints, with every join orientation seen there."""
    jid: str
    end1: tuple              # (chrom, breakpoint index on that chrom)
    end2: tuple
    adjacencies: dict = field(default_factory=dict)   # (side1, side2) -> {caller: reads}
    precise: dict = field(default_factory=dict)       # (side1, side2) -> n precise records
    types: dict = field(default_factory=dict)         # caller or 'junction table' -> set of svtypes
    rows: list = field(default_factory=list)          # junction-table rows (chrom1,pos1,chrom2,pos2,svtype)
    optional: bool = False                            # only composite records here: no cost if unused


# ─── input readers ───────────────────────────────────────────────────────────
def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def _info(s):
    d = {}
    for kv in s.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
        else:
            d[kv] = True
    return d


def record_ends(chrom, pos, alt, info):
    """[(chrom1, pos1, chrom2, pos2, [(s1, s2), ...], precise)] for one VCF record, or [].

    Positions are the break coordinates the callers write (POS / END / BND mate);
    within a few bp of the 0-based segment boundary, which is all clustering needs."""
    t = str(info.get("SVTYPE", ".")).upper()
    for rx, s_pos, s_mate in _BND:
        m = rx.match(alt)
        if m:
            return [(chrom, pos, m.group(1), int(m.group(2)), [(s_pos, s_mate)], True)]
    try:
        end = int(info.get("END", pos))
    except ValueError:
        end = pos
    c2 = info.get("CHR2", chrom) if isinstance(info.get("CHR2"), str) else chrom
    if t == "DEL":
        return [(chrom, pos, c2, end, [("+", "-")], True)]
    if t in ("DUP", "DUP:TANDEM"):
        return [(chrom, pos, c2, end, [("-", "+")], True)]
    if t == "INV":
        strand = str(info.get("STRAND", info.get("STRANDS", "")))
        if strand in ("++", "--") and "SVIM" not in str(info.get("_caller", "")):
            return [(chrom, pos, c2, end, [(strand[0], strand[1])], True)]
        return [(chrom, pos, c2, end, [("+", "+"), ("-", "-")], False)]
    return []


def read_caller_records(paths, windows, tol):
    """Caller records with both ends inside the locus windows.

    paths: {caller: filtered VCF}; windows: {chrom: [(lo, hi)]}. Returns [Record]."""
    out = []

    def inside(c, p):
        return any(lo - tol <= p <= hi + tol for lo, hi in windows.get(c, ()))

    for caller, path in paths.items():
        if not path or not os.path.exists(path):
            continue
        with _open(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.split("\t", 8)
                if len(f) < 8 or f[0] not in windows:
                    continue
                pos = int(f[1])
                if not inside(f[0], pos):
                    continue
                info = _info(f[7])
                info["_caller"] = caller
                ends = record_ends(f[0], pos, f[4], info)
                reads = None
                for k in READ_KEYS.get(caller, ("SUPPORT",)):
                    if k in info:
                        try:
                            reads = int(float(info[k]))
                        except (TypeError, ValueError):
                            pass
                        break
                for c1, p1, c2, p2, orients, precise in ends:
                    if inside(c1, p1) and inside(c2, p2):
                        out.append(Record(caller, f[2], str(info.get("SVTYPE", ".")), c1, p1, c2, p2,
                                          orients, reads, precise))
    return out


def read_merged_orientations(vcf_path, ids):
    """{record id: [(chrom1,pos1,chrom2,pos2,orients,precise)]} from the merged VCF (fallback)."""
    out = {}
    if not vcf_path or not os.path.exists(vcf_path) or not ids:
        return out
    with _open(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if len(f) < 8 or f[2] not in ids:
                continue
            info = _info(f[7])
            out[f[2]] = record_ends(f[0], int(f[1]), f[4], info)
    return out


def read_locus_depth(path, locus):
    """[(chrom, start, end, kind, mean depth)] of one locus from scoring/locus_depth.tsv."""
    rows = []
    if not path or not os.path.exists(path):
        return rows
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 7 and f[3] == locus:
                s, e = int(f[1]), int(f[2])
                rows.append((f[0], s, e, f[4], float(f[6]) / max(1, e - s)))
    return rows


def read_cn_segments(path):
    """cn/segments.tsv rows as dicts (chrom, start, end, state, ploidy, cn)."""
    out = []
    if not path or not os.path.exists(path):
        return out
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                out.append({"chrom": r["chrom"], "start": int(r["start"]), "end": int(r["end"]),
                            "state": int(r["state"]), "ploidy": int(r["ploidy"]), "cn": float(r["cn"])})
            except (KeyError, ValueError):
                continue
    return out


class GeneIndex:
    """Genes BED -> overlap queries. The bundle's genes BED merges overlapping genes into one
    interval with a comma-separated name (e.g. NF1 with EVI2A, EVI2B and OMG, which lie in an
    NF1 intron); such an interval is reported as one unit, never split into invented extents."""

    def __init__(self, path):
        self.rows = defaultdict(list)
        if path and os.path.exists(path):
            for c, s, e, g in sc.read_bed(path, split_labels=False):
                if g and g != ".":
                    self.rows[c].append((s, e, g))
        for c in self.rows:
            self.rows[c].sort()

    def overlapping(self, c, s, e):
        return [(gs, ge, g) for gs, ge, g in self.rows.get(c, ()) if gs < e and ge > s]


# ─── locus model ─────────────────────────────────────────────────────────────
def _natural(c):
    m = re.match(r"^(chr)?(\d+)$", c)
    if m:
        return (0, int(m.group(2)), c)
    return (1, 0, {"chrX": "X", "chrY": "Y", "chrM": "Z"}.get(c, c))


def cluster_breakpoints(ends, tol):
    """{chrom: [break positions]} merging ends within tol (leader-based, as cxsv_classes.breakends)."""
    by = defaultdict(list)
    for c, p in ends:
        by[c].append(p)
    out = {}
    for c, ps in by.items():
        lead, keep = None, []
        for p in sorted(ps):
            if lead is None or p - lead > tol:
                lead = p
                keep.append(p)
        out[c] = keep
    return out


def _nearest(bps, c, p, tol):
    best = None
    for i, b in enumerate(bps.get(c, ())):
        d = abs(b - p)
        if d <= tol and (best is None or d < best[1]):
            best = (i, d)
    return None if best is None else best[0]


def locus_junction_rows(junction_rows, locus):
    return [r for r in junction_rows if locus in r["loci"]]


def build_groups(kept, records, merged_fallback, bps, tol, optional=()):
    """Junction groups (one per breakpoint pair) with orientations and support.

    optional: junctions dropped as composite by the scorer whose two ends both lie at
    breakpoints of the kept junctions; they may be used by the walk but cost nothing unused."""
    groups = {}
    opt_ids = {id(j) for j in optional}
    kept = list(kept) + list(optional)

    def key_of(c1, i1, c2, i2):
        a, b = (c1, i1), (c2, i2)
        return (a, b) if (_natural(a[0]), a[1]) <= (_natural(b[0]), b[1]) else (b, a)

    def add(g, e1, e2, s1, s2, caller, reads, precise):
        # orient the side pair to the group's end order
        if (e1, e2) != (g.end1, g.end2):
            s1, s2 = s2, s1
        adj = g.adjacencies.setdefault((s1, s2), {})
        if caller is not None:
            adj[caller] = max(adj.get(caller) or 0, reads or 0) if reads is not None else adj.get(caller)
        g.precise[(s1, s2)] = g.precise.get((s1, s2), 0) + (1 if precise else 0)

    for j in kept:
        i1 = _nearest(bps, j.chrom1, j.pos1, tol)
        i2 = _nearest(bps, j.chrom2, j.pos2, tol)
        if i1 is None or i2 is None:
            continue
        e1, e2 = key_of(j.chrom1, i1, j.chrom2, i2)
        g = groups.get((e1, e2))
        if g is None:
            g = groups[(e1, e2)] = Group("", e1, e2, optional=id(j) in opt_ids)
        elif id(j) not in opt_ids:
            g.optional = False
        g.rows.append((j.chrom1, j.pos1, j.chrom2, j.pos2, j.svtype))
        g.types.setdefault("junction table", set()).add(j.svtype)
        matched = False
        for r in records:
            for (a_c, a_p, b_c, b_p) in ((r.chrom1, r.pos1, r.chrom2, r.pos2), (r.chrom2, r.pos2, r.chrom1, r.pos1)):
                if a_c == j.chrom1 and b_c == j.chrom2 and abs(a_p - j.pos1) <= tol and abs(b_p - j.pos2) <= tol:
                    swap = (a_c, a_p) != (r.chrom1, r.pos1)
                    for s1, s2 in r.orients:
                        if swap:
                            s1, s2 = s2, s1
                        add(g, (j.chrom1, i1), (j.chrom2, i2), s1, s2, r.caller, r.reads, r.precise)
                    g.types.setdefault(r.caller, set()).add(r.svtype)
                    matched = True
                    break
        if not matched:
            for rid in j.record_ids:
                for c1, p1, c2, p2, orients, precise in merged_fallback.get(rid, []):
                    swap = not (c1 == j.chrom1 and abs(p1 - j.pos1) <= tol)
                    for s1, s2 in orients:
                        if swap:
                            s1, s2 = s2, s1
                        add(g, (j.chrom1, i1), (j.chrom2, i2), s1, s2, None, None, precise)
                    matched = True
            if not matched:
                t = j.svtype.upper()
                default = {"DEL": [("+", "-")], "DUP": [("-", "+")], "INV": [("+", "+"), ("-", "-")]}.get(t, [])
                lo_first = j.chrom1 == j.chrom2 and j.pos1 <= j.pos2
                for s1, s2 in default:
                    if not lo_first:
                        s1, s2 = s2, s1
                    add(g, (j.chrom1, i1), (j.chrom2, i2), s1, s2, None, None, False)
    out = sorted(groups.values(), key=lambda g: ((_natural(g.end1[0]), g.end1[1]), (_natural(g.end2[0]), g.end2[1])))
    out = [g for g in out if g.adjacencies]
    for k, g in enumerate(out, 1):
        g.jid = f"J{k}"
    return out


def make_segments(bps, flank_bp=FLANK_BP):
    """Segments per chromosome: [{name, chrom, start, end, kind}] in chromosome then position order.

    Single chromosome: flanks p / q, internal A, B, C ... ; several chromosomes: flanks
    p_<chrom> / q_<chrom>, internal letters continue across chromosomes."""
    chroms = sorted(bps, key=_natural)
    multi = len(chroms) > 1
    segs, k = [], 0
    for c in chroms:
        b = bps[c]
        segs.append({"name": f"p_{c}" if multi else "p", "chrom": c, "start": max(0, b[0] - flank_bp),
                     "end": b[0], "kind": "p"})
        for i in range(len(b) - 1):
            segs.append({"name": _letter(k), "chrom": c, "start": b[i], "end": b[i + 1], "kind": "internal"})
            k += 1
        segs.append({"name": f"q_{c}" if multi else "q", "chrom": c, "start": b[-1], "end": b[-1] + flank_bp,
                     "kind": "q"})
    return segs


def _letter(k):
    s = ""
    k += 1
    while k:
        k, r = divmod(k - 1, 26)
        s = chr(65 + r) + s
    return s


# ─── copy number ─────────────────────────────────────────────────────────────
def segment_copy_numbers(segs, depth_rows, cn_rows, sex, par, cn_params):
    """Per internal segment: {cn, state, ploidy, source, m} (m = copies on the derivative)."""
    depth_cn = sc.cn_segments(depth_rows, sex, par, cn_params) if depth_rows else None
    out = {}
    for s in segs:
        if s["kind"] != "internal":
            continue
        L = s["end"] - s["start"]
        info = {"cn": None, "state": None, "ploidy": sc.ploidy(s["chrom"], (s["start"] + s["end"]) // 2, sex, par) or 1,
                "source": "not measurable", "m": None}
        for c, a, b, cn, state, pl, changed in depth_cn or ():
            ov = min(b, s["end"]) - max(a, s["start"])
            if c == s["chrom"] and ov >= 0.5 * min(L, b - a) and ov > 0:
                info.update(cn=cn, state=state, ploidy=pl, source="locus read depth")
                break
        if info["state"] is None and L >= int(cn_params.get("min_segment_bp", 1000)):
            for r in cn_rows:
                ov = min(r["end"], s["end"]) - max(r["start"], s["start"])
                if r["chrom"] == s["chrom"] and ov >= 0.5 * L:
                    info.update(cn=r["cn"], state=r["state"], ploidy=r["ploidy"], source="copy-number segment")
                    break
        if info["state"] is not None:
            info["m"] = max(0, int(info["state"]) - (int(info["ploidy"]) - 1))
        out[s["name"]] = info
    return out


def cn_boundary_notes(segs, cn_rows, bps, slack=CN_BOUNDARY_SLACK):
    """Copy-number segment boundaries inside / at the edge of the locus with no junction nearby."""
    notes = []
    chroms = {s["chrom"] for s in segs}
    for r in cn_rows:
        if r["chrom"] not in chroms or r["state"] == r["ploidy"]:
            continue
        b = bps[r["chrom"]]
        lo, hi = b[0] - slack, b[-1] + slack
        if r["end"] < lo or r["start"] > hi:
            continue
        for edge in (r["start"], r["end"]):
            d = min(abs(edge - x) for x in b)
            if d > slack:
                kind = "loss" if r["state"] < r["ploidy"] else "gain"
                notes.append(f"the read-depth copy-number {kind} {r['chrom']}:{r['start']:,}-{r['end']:,} "
                             f"(CN {r['cn']:.2f}) has a boundary at {edge:,}, {d / 1000:.1f} kb from the nearest "
                             f"junction; no junction was called there, so that stretch is not drawn as "
                             f"{'deleted' if kind == 'loss' else 'duplicated'}")
    return notes


# ─── graph walk ──────────────────────────────────────────────────────────────
class Graph:
    """Segment ends and their adjacencies. End (i, 0) = start of segment i, (i, 1) = its end."""

    def __init__(self, segs, groups, bps):
        self.segs = segs
        self.edges = defaultdict(list)            # end -> [(other end, edge id)]
        self.edge_info = {}                       # edge id -> ('ref', i) | ('adj', group index, sides)
        idx = {}
        for i, s in enumerate(segs):
            idx[(s["chrom"], s["start"], "start")] = i
            idx[(s["chrom"], s["end"], "end")] = i
        for i in range(len(segs) - 1):
            a, b = segs[i], segs[i + 1]
            if a["chrom"] == b["chrom"] and a["end"] == b["start"]:
                eid = f"r{i}"
                self.edge_info[eid] = ("ref", i)
                self._add((i, 1), (i + 1, 0), eid)
        for gi, g in enumerate(groups):
            (c1, i1), (c2, i2) = g.end1, g.end2
            x1, x2 = bps[c1][i1], bps[c2][i2]
            for (s1, s2) in sorted(g.adjacencies):
                e1 = self._end(idx, c1, x1, s1)
                e2 = self._end(idx, c2, x2, s2)
                if e1 is None or e2 is None:
                    continue
                eid = f"a{gi}{s1}{s2}"
                self.edge_info[eid] = ("adj", gi, (s1, s2))
                self._add(e1, e2, eid)
        for e in self.edges:
            self.edges[e].sort(key=lambda t: (t[1].startswith("r"), t[1], t[0]))

    @staticmethod
    def _end(idx, c, x, side):
        if side == "+":
            i = idx.get((c, x, "end"))
            return None if i is None else (i, 1)
        i = idx.get((c, x, "start"))
        return None if i is None else (i, 0)

    def _add(self, e1, e2, eid):
        self.edges[e1].append((e2, eid))
        if e1 != e2:
            self.edges[e2].append((e1, eid))


def cn_cost(cn, ploidy, v):
    """Cost of visiting a segment v times when read depth gives copy number cn (heterozygous model)."""
    return COST["cn_mismatch"] * max(0.0, abs(cn - (v + ploidy - 1)) - CN_SLACK)


def walk(segs, groups, bps, cn_of, support_of, allow_open=False, limit=SEARCH_LIMIT):
    """All lowest-cost derivative walks. Returns (solutions, n_expanded, hit_limit).

    cn_of: {segment name: (copy number, ploidy) or None}. A solution = dict(paths=[[(seg index,
    strand, edge id into it or None)]], cost, open_ends, visits, used_edges)."""
    G = Graph(segs, groups, bps)
    n = len(segs)
    flank_order = [i for i, s in enumerate(segs) if s["kind"] in ("p", "q")]
    cap, cnp, vbest = [], [], []
    for i, s in enumerate(segs):
        x = cn_of.get(s["name"]) if s["kind"] == "internal" else None
        cnp.append(x)
        if s["kind"] != "internal":
            cap.append(1)
            vbest.append(1)
        elif x is None:
            cap.append(2)
            vbest.append(1)
        else:
            vb = max(0, int(round(x[0] - x[1] + 1)))
            vbest.append(vb)
            cap.append(vb + 1)
    visits = [0] * n
    used = defaultdict(int)
    best = {"cost": None, "sols": []}
    counter = [0]
    hit = [False]

    def seg_cost(i, v):
        if segs[i]["kind"] != "internal":
            return 0.0
        if cnp[i] is None:
            return COST["unknown_cn"] * abs(v - 1)
        return cn_cost(cnp[i][0], cnp[i][1], v)

    def over_cost(i):
        # cost already certain: visits beyond the best count only add
        v = visits[i]
        return seg_cost(i, v) - seg_cost(i, v - 1) if v > vbest[i] else 0.0

    def final_cost(partial):
        cost = partial
        ug = {G.edge_info[e][1] for e, k in used.items() if k > 0 and G.edge_info[e][0] == "adj"}
        cost += COST["unused_junction"] * sum(1 for gi, g in enumerate(groups) if not g.optional and gi not in ug)
        for i in range(n):
            v = visits[i]
            if v <= vbest[i]:
                cost += seg_cost(i, v)
            else:
                cost += seg_cost(i, vbest[i])
        return round(cost, 6), ug

    def record(paths, partial, open_ends):
        cost, ug = final_cost(partial)
        if best["cost"] is not None and cost > best["cost"] + 1e-9:
            return
        sol = {"paths": [list(p) for p in paths], "cost": cost, "open_ends": list(open_ends),
               "visits": list(visits), "used_groups": sorted(ug),
               "used_edges": sorted(e for e, k in used.items() if k > 0 and G.edge_info[e][0] == "adj")}
        if best["cost"] is None or cost < best["cost"] - 1e-9:
            best["cost"], best["sols"] = cost, [sol]
        elif len(best["sols"]) < MAX_SOLUTIONS:
            best["sols"].append(sol)

    def start_next(paths, partial, open_ends):
        nxt = next((i for i in flank_order if visits[i] == 0), None)
        if nxt is None:
            record(paths, partial, open_ends)
            return
        strand = "+" if segs[nxt]["kind"] == "p" else "-"
        enter(paths + [[]], nxt, strand, None, partial, open_ends)

    def enter(paths, i, strand, eid, partial, open_ends):
        counter[0] += 1
        if counter[0] > limit:
            hit[0] = True
            return
        if visits[i] >= cap[i]:
            return
        visits[i] += 1
        partial2 = partial + over_cost(i)
        if best["cost"] is None or partial2 <= best["cost"] + 1e-9:
            paths[-1].append((i, strand, eid))
            exit_end = (i, 1) if strand == "+" else (i, 0)
            telomeric = (segs[i]["kind"] == "p" and exit_end == (i, 0)) or \
                        (segs[i]["kind"] == "q" and exit_end == (i, 1))
            if telomeric:
                if len(paths[-1]) > 1:
                    start_next(paths, partial2, open_ends)
            else:
                for (f, e) in G.edges.get(exit_end, ()):
                    j = f[0]
                    if visits[j] >= cap[j]:
                        continue
                    used[e] += 1
                    enter(paths, j, "+" if f[1] == 0 else "-", e, partial2, open_ends)
                    used[e] -= 1
                    if hit[0]:
                        break
                if allow_open and not hit[0]:
                    p2 = partial2 + COST["open_end"]
                    if best["cost"] is None or p2 <= best["cost"] + 1e-9:
                        start_next(paths, p2, open_ends + [(i, strand)])
            paths[-1].pop()
        visits[i] -= 1

    start_next([], 0.0, [])
    return best["sols"], counter[0], hit[0], G


def reverse_path(p):
    """The same derivative read from its other end (strands flipped, edge labels moved)."""
    n = len(p) - 1
    out = []
    for k in range(n + 1):
        i, s, _ = p[n - k]
        e = p[n - k + 1][2] if k > 0 else None
        out.append((i, "-" if s == "+" else "+", e))
    return out


def orient_path(segs, p):
    """Read a derivative from a p flank in reference orientation when possible."""
    def good(q):
        f = segs[q[0][0]]
        return f["kind"] == "p" and q[0][1] == "+"
    if good(p):
        return p
    r = reverse_path(p)
    if good(r):
        return r
    if segs[p[0][0]]["kind"] == "q" and p[0][1] == "-":
        return r
    return p


def canonical(segs, paths):
    """Paths as tuples of (name, strand), each read from its p flank (or the smaller end), sorted."""
    out = []
    for p in paths:
        fw = tuple((segs[i]["name"], s) for i, s, _ in p)
        rv = tuple((segs[i]["name"], "-" if s == "+" else "+") for i, s, _ in reversed(p))
        def key(t):
            first = t[0] if t else ("", "+")
            return (not (first[0].startswith("p") and first[1] == "+"), t)
        out.append(min(fw, rv, key=key))
    return tuple(sorted(out))


def path_strings(segs, paths, open_ends=()):
    return [" ".join(f"{segs[i]['name']}{s}" for i, s, _ in p) for p in paths]


# ─── main entry ──────────────────────────────────────────────────────────────
def reconstruct(locus_row, junction_rows, records, merged_fallback, depth_rows, cn_rows, genes, params,
                sex="unknown", par=()):
    """Reconstruction of one scored locus. Returns a dict (see module docstring)."""
    L = locus_row["Locus_ID"]
    bp_tol = max(int(params.get("collapse_distance", 300)), BREAKPOINT_TOL)
    key_genes = set(re.findall(r"([^;\s(]+):\w+\(HI [123],", str(locus_row.get("gene_dosage_value", ""))))
    tol = int(params.get("collapse_distance", 300))
    cls_p = {**cc.DEFAULTS, "collapse_distance": tol,
             "composite_fraction": float(params.get("composite_fraction", cc.DEFAULTS["composite_fraction"]))}
    jr = locus_junction_rows(junction_rows, L)
    pool = [cj.Junction(r["chrom1"], r["pos1"], r["chrom2"], r["pos2"], r["svtype"], r["svlen"], r["supp_vec"],
                        r["dosage"], r["record_ids"].split(",")) for r in jr]
    order = {c: k for k, c in enumerate(sorted({x for j in pool for x in (j.chrom1, j.chrom2)}, key=_natural))}
    kept, composite = cc.distinct_junctions(pool, order, cls_p) if pool else ([], [])
    kept = [j for j in kept if not (j.svtype == "INS")]
    res = {"locus": L, "status": "resolved", "reasons": [], "notes": [], "alternatives": [],
           "segments": [], "junctions": [], "paths": {}, "segment_info": {}, "inter": False,
           "composite": [f"{j.svtype} {j.chrom1}:{j.pos1:,}-{j.pos2:,} ({j.supp_vec})" for j in composite]}
    if not kept:
        res["status"] = "ambiguous"
        res["reasons"].append("no junction of this locus could be placed (copy-number-only or insertion-only call); "
                              "no structure is drawn")
        return res
    bps = cluster_breakpoints([(j.chrom1, j.pos1) for j in kept] + [(j.chrom2, j.pos2) for j in kept], bp_tol)
    at_bps = [j for j in composite if _nearest(bps, j.chrom1, j.pos1, bp_tol) is not None
              and _nearest(bps, j.chrom2, j.pos2, bp_tol) is not None]
    groups = build_groups(kept, records, merged_fallback, bps, bp_tol, optional=at_bps)
    segs = make_segments(bps)
    res["inter"] = len(bps) > 1
    cninfo = segment_copy_numbers(segs, depth_rows, cn_rows, sex, par, params.get("copy_number", {}))
    m_of = {k: v["m"] for k, v in cninfo.items()}
    cn_of = {k: (None if v["cn"] is None else (float(v["cn"]), int(v["ploidy"]))) for k, v in cninfo.items()}

    sols, n_exp, hit, G = walk(segs, groups, bps, cn_of, None)
    if not sols or min(s["cost"] for s in sols) >= COST["unused_junction"]:
        sols2, n2, hit2, G = walk(segs, groups, bps, cn_of, None, allow_open=True)
        if sols2 and (not sols or min(s["cost"] for s in sols2) < min(s["cost"] for s in sols) - 1e-9):
            sols, hit = sols2, hit or hit2
    ranked = rank_solutions(segs, groups, G, sols)
    res["segments"] = segs
    res["segment_info"] = cninfo
    res["groups"] = groups
    if not ranked:
        res["status"] = "ambiguous"
        res["reasons"].append("no derivative walk was found through the called junctions")
        _describe(res, segs, groups, bps, genes, None, key_genes)
        return res
    top = ranked[0]
    ties = [s for s in ranked if s["tie_key"] == top["tie_key"]]
    res["best"] = top
    res["paths"] = _name_paths(segs, top["paths"])
    res["alternatives"] = [list(_name_paths(segs, s["paths"]).values()) for s in ties[1:6]]
    res["n_equivalent"] = len(ties)
    if len(ties) > 1:
        res["status"] = "ambiguous"
        res["reasons"].append(f"{len(ties)} derivative orders are equally consistent with the junctions and read depth "
                              "(the same joins taken in a different order); the first is drawn, the others are listed")
    flips = [s for s in ranked if s["tie_key"][:3] == top["tie_key"][:3] and s["arm_flips"] > top["arm_flips"]]
    if flips:
        res["notes"].append(f"{len(flips)} further order(s) fit equally well but join two left (or two right) "
                            "chromosome ends, which gives a dicentric or acentric chromosome when the breaks lie on the "
                            "same side of both centromeres; not drawn (the bundle has no centromere positions)")
    unused = [g for gi, g in enumerate(groups) if gi not in top["used_groups"] and not g.optional]
    bridged = {g.jid: br for g in unused for br in [_bridged(g, segs, bps, top)] if br}
    for jid, br in bridged.items():
        res["notes"].append(f"junction {jid} (deletion-type, {br[0]}) spans the inserted {br[1]} "
                            f"({br[2] / 1000:.1f} kb) on the derivative: a caller's summary of a deletion with an "
                            "insertion, not drawn as a separate join")
    unused = [g for g in unused if g.jid not in bridged]
    for g in unused:
        res["status"] = "ambiguous"
        res["reasons"].append(f"junction {g.jid} is not on the drawn derivative (no walk uses it together with the "
                              "others; a partner junction may be missing, it may lie on the other haplotype, or it "
                              "is a caller's summary of a more complex join)")
    if top["open_ends"]:
        res["status"] = "ambiguous"
        ends = ", ".join(f"after {segs[i]['name']}{s}" for i, s in top["open_ends"])
        res["reasons"].append(f"the derivative could not be completed ({ends}): a junction is probably missing; "
                              "the pieces are drawn joined at a hollow marker (no junction called there)")
    for i, s in enumerate(segs):
        if s["kind"] != "internal" or cn_of.get(s["name"]) is None:
            continue
        cn, pl = cn_of[s["name"]]
        v = top["visits"][i]
        if abs(cn - (v + pl - 1)) > CN_FLAG:
            res["status"] = "ambiguous"
            res["reasons"].append(f"segment {s['name']}: read depth gives CN {cn:.2f}, but the drawn derivative "
                                  f"carries it {v} time{'s' if v != 1 else ''} (CN {v + pl - 1} expected)")
    for gi, g in enumerate(groups):
        kinds = {JOIN_KIND[a] for a in g.adjacencies}
        if {"del", "inv"} <= kinds or {"dup", "inv"} <= kinds:
            (c1, i1), (c2, i2) = g.end1, g.end2
            spanned = [s for s in segs if s["kind"] == "internal" and s["chrom"] == c1 == c2
                       and bps[c1][i1] <= s["start"] and s["end"] <= bps[c2][i2]]
            if spanned and all(cn_of.get(s["name"]) is None for s in spanned):
                res["status"] = "ambiguous"
                res["reasons"].append(f"junction {g.jid}: the callers disagree on the join (deletion- or duplication-"
                                      "type vs inversion-type) and read depth cannot decide (spanned sequence too short "
                                      "to measure); the drawn choice rests on parsimony only")
    if hit:
        res["status"] = "ambiguous"
        res["reasons"].append("walk search limit reached; the drawn derivative may not be the best one")
    res["notes"] += cn_boundary_notes(segs, cn_rows, bps)
    _describe(res, segs, groups, bps, genes, top, key_genes)
    return res


def rank_solutions(segs, groups, G, sols):
    """Distinct solutions (oriented, fragments joined) sorted by
    (cost, -support, inferred joins used, arm flips, inverted bp, canonical form).

    support = per junction used, the callers seeing the orientation used (+0.5 per precise
    record); an *inferred* join is an orientation only implied by a symbolic INV record (no
    record states it), so a walk that needs fewer of them is preferred. Solutions equal in
    (cost, support, inferred joins, arm flips) are ties: the locus is then ambiguous."""
    seen = {}
    for s in sols:
        s = dict(s)
        if s["open_ends"]:
            s = _join_fragments(segs, s)
        s["paths"] = [orient_path(segs, p) for p in s["paths"]]
        key = canonical(segs, s["paths"])
        sup, inferred = 0.0, 0
        best_per_group = {}
        for e in s["used_edges"]:
            _, gi, sides = G.edge_info[e]
            g = groups[gi]
            w = len(g.adjacencies.get(sides, {})) + 0.5 * g.precise.get(sides, 0)
            best_per_group[gi] = max(best_per_group.get(gi, 0), w)
            if g.precise.get(sides, 0) == 0:
                inferred += 1
        sup = sum(best_per_group.values())
        flips = 0
        for p in s["paths"]:
            a, b = segs[p[0][0]], segs[p[-1][0]]
            if a["kind"] in ("p", "q") and b["kind"] in ("p", "q") and a["chrom"] != b["chrom"]:
                flips += int(a["kind"] == b["kind"])          # two left or two right chromosome ends
        inv_bp = sum(segs[i]["end"] - segs[i]["start"] for p in s["paths"] for i, st, _ in p
                     if st == "-" and segs[i]["kind"] == "internal")
        s["support"], s["inferred"], s["arm_flips"] = sup, inferred, flips
        s["tie_key"] = (round(s["cost"], 3), -round(sup, 3), inferred, flips)
        rank = s["tie_key"] + (inv_bp, key)
        if key not in seen or rank < seen[key][0]:
            seen[key] = (rank, s)
    return [v[1] for v in sorted(seen.values(), key=lambda v: v[0])]


def _join_fragments(segs, sol):
    """Derivatives from the fragments of an open-ended walk: each fragment that stopped at an
    open end is joined to the next such fragment read backwards, at a 'gap' (no junction
    called). Complete fragments are kept as they are."""
    opened = {(i, s) for i, s in sol["open_ends"]}
    done, frags = [], []
    for p in sol["paths"]:
        (frags if (p[-1][0], p[-1][1]) in opened else done).append(list(p))
    while len(frags) >= 2:
        a, b = frags.pop(0), frags.pop(0)
        rb = reverse_path(b)
        done.append(a + [(rb[0][0], rb[0][1], "gap")] + rb[1:])
    done += frags
    out = dict(sol)
    out["paths"] = done
    return out


BRIDGE_MAX_BP = 10000


def _bridged(g, segs, bps, sol):
    """(span text, inserted names, inserted bp) when an unused deletion-type junction joins two
    segment ends that the drawn derivative connects through a short inserted piece, else None."""
    if set(g.adjacencies) != {("+", "-")}:
        return None
    (c1, i1), (c2, i2) = g.end1, g.end2
    x1, x2 = bps[c1][i1], bps[c2][i2]
    for p in sol["paths"]:
        for a in range(len(p)):
            sa = segs[p[a][0]]
            exit_pos = (sa["chrom"], sa["end"]) if p[a][1] == "+" else None
            if exit_pos != (c1, x1):
                continue
            ins = 0
            for b in range(a + 1, len(p)):
                sb = segs[p[b][0]]
                if p[b][1] == "+" and (sb["chrom"], sb["start"]) == (c2, x2):
                    if b > a + 1 and ins <= BRIDGE_MAX_BP:
                        names = " ".join(f"{segs[i]['name']}{s}" for i, s, _ in p[a + 1:b])
                        return (f"{c1}:{x1:,}-{x2:,}", names, ins)
                    break
                ins += sb["end"] - sb["start"]
                if ins > BRIDGE_MAX_BP:
                    break
    return None


def _name_paths(segs, paths):
    """{derivative name: ["p+", "A-", ...]} (PC events.yaml notation)."""
    out = {}
    chroms = sorted({s["chrom"] for s in segs}, key=_natural)
    for p in paths:
        names = [f"{segs[i]['name']}{s}" for i, s, _ in p]
        if len(chroms) == 1:
            key = "der" if "der" not in out else f"der{len(out) + 1}"
        else:
            first = segs[p[0][0]]
            c = first["chrom"].replace("chr", "")
            ref = [f"{s['name']}+" for s in segs if s["chrom"] == first["chrom"]]
            key = f"normal({c})" if names == ref else f"der({c})"
            if key in out:
                key = f"{key}'"
        out[key] = names
    return out


def _describe(res, segs, groups, bps, genes, top, key_genes=()):
    """Junction table, segment genes, caller-type notes and the plain-language inference."""
    used_adj = defaultdict(set)
    if top:
        for p in top["paths"]:
            for i, s, e in p:
                if e and e.startswith("a"):
                    gi = int(re.match(r"a(\d+)", e).group(1))
                    used_adj[gi].add(e[len(f"a{gi}"):])
    J = []
    for gi, g in enumerate(groups):
        (c1, i1), (c2, i2) = g.end1, g.end2
        x1, x2 = bps[c1][i1], bps[c2][i2]
        callers = defaultdict(int)
        for sides, adj in g.adjacencies.items():
            for k, v in adj.items():
                callers[k] = max(callers[k], v or 0)
        orients = sorted(g.adjacencies)
        used = sorted(used_adj.get(gi, ()))
        J.append({"id": g.jid, "end1": f"{c1}:{x1:,}", "end2": f"{c2}:{x2:,}",
                  "orientations": [f"{a}{b} {JOIN_NAMES.get((a, b), '')}" for a, b in orients],
                  "used": [f"{u} {JOIN_NAMES.get((u[0], u[1]), '')}" for u in used],
                  "callers": dict(callers), "types": {k: sorted(v) for k, v in g.types.items()},
                  "max_reads": max(callers.values()) if callers else None})
        if label_conflict(g, used):
            res["notes"].append(_type_note(g, J[-1], segs, res["segment_info"], bps))
    res["junctions"] = J
    # genes per internal segment, and genes broken at a breakpoint
    seg_genes = {}
    for s in segs:
        if s["kind"] != "internal":
            continue
        gl = []
        for gs, ge, gname in genes.overlapping(s["chrom"], s["start"], s["end"]) if genes else ():
            whole = gs >= s["start"] and ge <= s["end"]
            orig = gname.split(",")
            names = sorted(orig, key=lambda x: (x not in key_genes, orig.index(x)))
            ov = min(ge, s["end"]) - max(gs, s["start"])
            gl.append(("/".join(names), "whole" if whole else f"part: {ov / 1000:.1f} of {(ge - gs) / 1000:.1f} kb"))
        seg_genes[s["name"]] = gl
    res["segment_genes"] = seg_genes
    broken = []
    for c, ps in bps.items():
        for x in ps:
            for gs, ge, gname in genes.overlapping(c, x - 1, x) if genes else ():
                if gs < x < ge:
                    orig = gname.split(",")
                    nm = "/".join(sorted(orig, key=lambda g: (g not in key_genes, orig.index(g))))
                    if nm not in broken:
                        broken.append(nm)
    res["broken_genes"] = broken
    res["summary"] = summarise(res, segs, top)


EXPECTED_JOINS = {"DEL": {"+-"}, "DUP": {"-+"}, "INV": {"++", "--"}}


def label_conflict(g, used):
    """True when the callers' type labels disagree with each other or with the join used."""
    labels = {t.upper() for v in g.types.values() for t in v} & set(EXPECTED_JOINS)
    if len(labels) > 1:
        return True
    return bool(used) and any(not (EXPECTED_JOINS[t] & set(used)) for t in labels)


def _type_note(g, j, segs, cninfo, bps):
    types = "; ".join(f"{k} {'/'.join(v)}" for k, v in g.types.items())
    used = ", ".join(j["used"]) or "none"
    (c1, i1), (c2, i2) = g.end1, g.end2
    spanned = [s for s in segs if s["kind"] == "internal" and s["chrom"] == c1 == c2
               and bps[c1][i1] <= s["start"] and s["end"] <= bps[c2][i2]]
    cn_txt = ", ".join(f"{s['name']} CN {cninfo[s['name']]['cn']:.2f}" for s in spanned
                       if cninfo.get(s["name"], {}).get("cn") is not None)
    txt = (f"{g.jid}: the callers label this junction differently ({types}); the orientation written in the "
           f"records is {', '.join(j['orientations'])}; the derivative uses it as {used}")
    if cn_txt:
        txt += f"; read depth of the spanned segment(s): {cn_txt}"
    del_like = any(u.startswith("+-") for u in j["used"])
    lost = spanned and all((cninfo.get(s["name"], {}).get("m") == 0) for s in spanned)
    if del_like and lost:
        txt += (". Read depth shows the spanned sequence lost on one copy, which supports a deletion: an inversion "
                "would keep the copy number")
    elif any(u[:2] in ("++", "--") for u in j["used"]) and spanned and all(
            cninfo.get(s["name"], {}).get("m") == 1 for s in spanned):
        txt += ". Read depth shows no copy-number change across it, as an inversion would"
    return txt


def summarise(res, segs, top):
    """One plain-language line: which segments are deleted, duplicated, inverted."""
    if not top:
        return "structure not determined"
    visits = top["visits"]
    inv = {segs[i]["name"] for p in top["paths"] for i, s, _ in p if s == "-" and segs[i]["kind"] == "internal"}
    parts = []
    def sz(s):
        n = s["end"] - s["start"]
        return f"{n / 1e6:.2f} Mb" if n >= 1e6 else (f"{n / 1e3:.1f} kb" if n >= 1e3 else f"{n} bp")
    dele = [f"{s['name']} ({sz(s)})" for i, s in enumerate(segs) if s["kind"] == "internal" and visits[i] == 0]
    dup = [f"{s['name']} x{visits[i]}" for i, s in enumerate(segs) if s["kind"] == "internal" and visits[i] > 1]
    if dele:
        parts.append(("deletion of " if len(dele) == 1 else "deletions of ") + ", ".join(dele))
    if dup:
        parts.append("extra copies of " + ", ".join(dup))
    if inv:
        parts.append("inverted " + ", ".join(sorted(inv)))
    if res.get("inter"):
        parts.append("joins between chromosomes " + "/".join(sorted({s["chrom"] for s in segs}, key=_natural)))
    if not parts:
        parts.append("no change of order, orientation or copy number")
    return "; ".join(parts)


def layout(res, title=None, subtitle=None):
    """The reconstruction as a drawing dict for cxsv_recon_figure.draw_event (PC events.yaml shape)."""
    segs = res["segments"]
    ev = {"id": res["locus"], "name": "", "title": title, "segments": {s["name"]: [s["chrom"], s["start"], s["end"]] for s in segs},
          "paths": res["paths"], "subtitle": subtitle}
    marks = {}
    if res.get("best"):
        names = list(res["paths"])
        gid = {f"a{gi}": g.jid for gi, g in enumerate(res["groups"])}
        for name, p in zip(names, res["best"]["paths"]):
            labs = []
            for k, (i, s, e) in enumerate(p):
                if e and e.startswith("a"):
                    labs.append(gid.get(re.match(r"(a\d+)", e).group(1)))
                else:
                    labs.append(None)
            marks[name] = labs
        ev["junction_labels"] = marks
        ev["gap_marks"] = {name: [k for k, (i, s, e) in enumerate(p) if e == "gap"]
                           for name, p in zip(names, res["best"]["paths"])}
    ev["segment_notes"] = {}
    for n, info in res["segment_info"].items():
        g = res["segment_genes"].get(n, [])
        gt = ", ".join(x for x, _ in g[:3]) + (f" +{len(g) - 3}" if len(g) > 3 else "")
        ev["segment_notes"][n] = {"cn": None if info["cn"] is None else f"CN {info['cn']:.1f}", "genes": gt}
    ev["status"] = None if res["status"] == "resolved" else "ambiguous"
    ev["title_line2"] = "inferred: " + textwrap.shorten(res.get("summary", ""), 110)
    return ev


def flat_rows(sample, rank, res):
    """TSV rows (one per segment and one per junction) for report/<sample>.tier1_reconstruction.tsv."""
    rows = []
    for s in res.get("segments", []):
        info = res["segment_info"].get(s["name"], {})
        rows.append({"sample": sample, "rank": rank, "locus": res["locus"], "status": res["status"],
                     "item": "segment", "name": s["name"], "where": f"{s['chrom']}:{s['start']}-{s['end']}",
                     "detail": (f"{s['kind']}; CN {info['cn']:.2f} ({info['source']}); copies on derivative "
                                f"{info['m'] if info.get('m') is not None else 'unknown'}"
                                if s["kind"] == "internal" and info.get("cn") is not None
                                else s["kind"] + ("; CN not measurable" if s["kind"] == "internal" else "")),
                     "genes": ",".join(f"{g}({w})" for g, w in res.get("segment_genes", {}).get(s["name"], []))})
    for j in res.get("junctions", []):
        rows.append({"sample": sample, "rank": rank, "locus": res["locus"], "status": res["status"], "item": "junction",
                     "name": j["id"], "where": f"{j['end1']}|{j['end2']}",
                     "detail": f"orientations {'/'.join(j['orientations'])}; used {'/'.join(j['used']) or 'no'}; "
                               + "; ".join(f"{k} {v} reads" for k, v in sorted(j["callers"].items())),
                     "genes": ""})
    for name, p in res.get("paths", {}).items():
        rows.append({"sample": sample, "rank": rank, "locus": res["locus"], "status": res["status"], "item": "path",
                     "name": name, "where": "", "detail": " ".join(p), "genes": ""})
    for k, alt in enumerate(res.get("alternatives", []), 2):
        rows.append({"sample": sample, "rank": rank, "locus": res["locus"], "status": res["status"],
                     "item": "alternative", "name": f"order {k}", "where": "", "detail": " | ".join(" ".join(p) for p in alt),
                     "genes": ""})
    for r in res.get("reasons", []):
        rows.append({"sample": sample, "rank": rank, "locus": res["locus"], "status": res["status"], "item": "ambiguity",
                     "name": "", "where": "", "detail": r, "genes": ""})
    for r in res.get("notes", []):
        rows.append({"sample": sample, "rank": rank, "locus": res["locus"], "status": res["status"], "item": "note",
                     "name": "", "where": "", "detail": r, "genes": ""})
    return rows


TSV_COLUMNS = ("sample", "rank", "locus", "status", "item", "name", "where", "detail", "genes")


def locus_windows(junction_rows, locus, pad=2000):
    """{chrom: [(lo, hi)]} around the locus's junction ends (for reading caller records)."""
    w = defaultdict(list)
    for r in locus_junction_rows(junction_rows, locus):
        w[r["chrom1"]].append((r["pos1"] - pad, r["pos1"] + pad))
        w[r["chrom2"]].append((r["pos2"] - pad, r["pos2"] + pad))
    return dict(w)


def reconstruct_run(t1_rows, junctions_path, caller_vcfs, merged_vcf, locus_depth_path, cn_segments_path,
                    genes_path, params, sex="unknown", par=()):
    """Reconstructions of several loci of one run (reads each input once). Returns [res]."""
    jrows = cj.read_table(junctions_path) if junctions_path and os.path.exists(junctions_path) else []
    cn_rows = read_cn_segments(cn_segments_path)
    genes = GeneIndex(genes_path)
    tol = max(int(params.get("collapse_distance", 300)), 500)
    windows = defaultdict(list)
    for r in t1_rows:
        for c, ivs in locus_windows(jrows, r["Locus_ID"]).items():
            windows[c] += ivs
    records = read_caller_records(caller_vcfs, dict(windows), tol)
    ids = {i for r in t1_rows for j in locus_junction_rows(jrows, r["Locus_ID"]) for i in j["record_ids"].split(",")}
    merged = read_merged_orientations(merged_vcf, ids)
    out = []
    for r in t1_rows:
        w = locus_windows(jrows, r["Locus_ID"])
        rec = [x for x in records if any(lo - tol <= x.pos1 <= hi + tol for lo, hi in w.get(x.chrom1, ()))
               and any(lo - tol <= x.pos2 <= hi + tol for lo, hi in w.get(x.chrom2, ()))]
        depth = read_locus_depth(locus_depth_path, r["Locus_ID"])
        try:
            out.append(reconstruct(r, jrows, rec, merged, depth, cn_rows, genes, params, sex, par))
        except Exception as exc:                         # never let one locus stop the report
            out.append({"locus": r["Locus_ID"], "status": "failed", "reasons": [f"reconstruction failed: {exc}"],
                        "notes": [], "segments": [], "junctions": [], "paths": {}, "segment_info": {},
                        "alternatives": []})
    return out
