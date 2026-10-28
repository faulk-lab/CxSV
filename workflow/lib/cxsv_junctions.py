"""Junction model: turn a sample's merged SV calls into breakpoint junctions.

A *junction* is one novel adjacency in a genome, written as two reference
positions joined by a structural variant:

    (chrom1, pos1) -- svtype --> (chrom2, pos2)      1-based, hs1 contig names

with the conventions of the shared junction-site schema (docs/JUNCTIONS.md):

    DEL / DUP / INV   chrom2 = chrom1, pos2 = END
    INS               pos2 = pos1, svlen = inserted length
    BND               (chrom2, pos2) = the mate; SURVIVOR's TRA is written BND
    order             (chrom1, pos1) <= (chrom2, pos2) in .fai contig order, then position

Per sample, every record of merged_sv.vcf becomes one junction (`parse_vcf`).
`collapse` then merges records that describe the same junction twice (for
matching only; the CxSV criteria in classify_cxsv.py are untouched):

  * Sniffles2 writes each inversion breakpoint as two reciprocal breakend
    records; SURVIVOR keeps both, e.g. INV 51213665->51223040 and INV
    51223041->51213667. After ordering they are the same junction.
  * a BND and its mate record (A->B and B->A).

Two records are merged when they have the same svtype and contigs and both ends
lie within `collapse_distance` bp (and, for INS, comparable lengths). Merging is
leader-based (each record is compared with the first record of a group, never
with a neighbour of a neighbour), so it cannot chain.

This module defines which junctions a sample has, so it is part of the calling
fingerprint (workflow/lib/cxsv_fingerprint.py: CALLING_CODE). Clustering across
individuals and patient matching live in cxsv_matching.py (not fingerprinted).
"""
from __future__ import annotations

import bisect
import re
from dataclasses import dataclass, field

SVTYPES = ("DEL", "DUP", "INV", "INS", "BND")
_TYPE_MAP = {"TRA": "BND", "DUP:TANDEM": "DUP", "DUP:INT": "DUP", "CNV": "DUP"}
_BND_RE = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")

# SURVIVOR writes one sample column per input VCF, in file-list order
# (rules/per_sample.smk: cuteSV, SVIM, Sniffles2). Genotype dosage is read from
# the first caller in this priority that has a called GT for the record.
GT_PRIORITY = (2, 0, 1)          # Sniffles2, cuteSV, SVIM
CALLER_NAMES = ("cuteSV", "SVIM", "Sniffles2")

JUNCTION_COLUMNS = ("chrom1", "pos1", "chrom2", "pos2", "svtype", "svlen",
                    "locus_id", "supp_vec", "dosage", "n_records", "record_ids")


@dataclass
class Junction:
    chrom1: str
    pos1: int
    chrom2: str
    pos2: int
    svtype: str
    svlen: int
    supp_vec: str = "."
    dosage: int = 1
    record_ids: list = field(default_factory=list)
    locus_ids: set = field(default_factory=set)

    def key(self, order):
        return (self.svtype, order.get(self.chrom1, 1 << 30), self.chrom1,
                order.get(self.chrom2, 1 << 30), self.chrom2, self.pos1, self.pos2, self.svlen)

    def as_row(self):
        return {"chrom1": self.chrom1, "pos1": self.pos1, "chrom2": self.chrom2, "pos2": self.pos2,
                "svtype": self.svtype, "svlen": self.svlen,
                "locus_id": ",".join(sorted(self.locus_ids)) or ".",
                "supp_vec": self.supp_vec, "dosage": self.dosage,
                "n_records": len(self.record_ids), "record_ids": ",".join(self.record_ids)}


# ─── contig order ────────────────────────────────────────────────────────────
def contig_order(fai_path=None, fallback=()):
    """{contig: rank} from a .fai (column 1, file order). Contigs not in the .fai
    (only possible with a foreign VCF) sort after all known ones, by name."""
    order = {}
    if fai_path:
        with open(fai_path) as fh:
            for line in fh:
                if line.strip():
                    order.setdefault(line.split("\t", 1)[0], len(order))
    for c in fallback:
        order.setdefault(c, len(order))
    return order


def ordered_ends(c1, p1, c2, p2, order):
    """Return the two ends with (chrom1, pos1) <= (chrom2, pos2)."""
    k1 = (order.get(c1, 1 << 30), c1, p1)
    k2 = (order.get(c2, 1 << 30), c2, p2)
    return (c1, p1, c2, p2) if k1 <= k2 else (c2, p2, c1, p1)


# ─── VCF parsing ─────────────────────────────────────────────────────────────
def _info(s):
    d = {}
    for kv in s.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
        elif kv:
            d[kv] = True
    return d


def _int(x, default=None):
    try:
        return int(str(x).split(",")[0])
    except (TypeError, ValueError):
        return default


def dosage_from_samples(fmt, samples):
    """Allele dosage (1 or 2) from the per-caller genotype columns, or 1 when none
    is called (a merged record exists, so the sample carries at least one allele)."""
    keys = fmt.split(":") if fmt else []
    if "GT" not in keys:
        return 1
    gi = keys.index("GT")
    for idx in GT_PRIORITY:
        if idx >= len(samples):
            continue
        f = samples[idx].split(":")
        gt = f[gi] if gi < len(f) else "./."
        alleles = [a for a in re.split(r"[/|]", gt) if a not in (".", "")]
        if not alleles:
            continue
        alt = sum(a != "0" for a in alleles)
        if alt == 0:
            continue
        return 2 if alt >= 2 else 1
    return 1


def record_to_junction(fields, order):
    """One VCF record (split on tabs) -> Junction, or None if it cannot be resolved."""
    chrom, pos, vid, alt, info_s = fields[0], int(fields[1]), fields[2], fields[4], fields[7]
    info = _info(info_s)
    raw = str(info.get("SVTYPE", "")).upper()
    svtype = _TYPE_MAP.get(raw, raw)
    if svtype not in SVTYPES:
        return None
    svlen = abs(_int(info.get("SVLEN"), 0) or 0)
    if svtype == "BND":
        m = _BND_RE.search(alt)
        if m:
            c2, p2 = m.group(1), int(m.group(2))
        else:
            c2, p2 = info.get("CHR2"), _int(info.get("END2", info.get("END")))
            if not c2 or p2 is None:
                return None
        svlen = 0
    elif svtype == "INS":
        c2, p2 = chrom, pos
    else:
        c2 = chrom
        p2 = _int(info.get("END"))
        if p2 is None or p2 == pos:
            p2 = pos + svlen if svlen else None
        if p2 is None:
            return None
        if not svlen:
            svlen = abs(p2 - pos)
    c1, p1, c2, p2 = ordered_ends(chrom, pos, c2, p2, order)
    supp = str(info.get("SUPP_VEC", "."))
    dose = dosage_from_samples(fields[8], fields[9:]) if len(fields) > 9 else 1
    return Junction(c1, p1, c2, p2, svtype, svlen, supp, dose, [vid])


def parse_vcf(path, order):
    """All resolvable records of a (merged) SV VCF as Junctions, in file order.
    Contigs from ##contig lines extend `order` when they are missing from it."""
    out = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("##contig=<ID="):
                name = line[len("##contig=<ID="):].split(",", 1)[0].rstrip(">\n")
                order.setdefault(name, len(order))
                continue
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            j = record_to_junction(f, order)
            if j is not None:
                out.append(j)
    return out


# ─── collapsing duplicate / reciprocal records ───────────────────────────────
def same_junction(a, b, distance, ins_len_ratio=0.5):
    """True when a and b are one junction: same type and contigs, both ends within
    `distance` bp and, for INS, lengths within `ins_len_ratio` (shorter/longer)."""
    if a.svtype != b.svtype or a.chrom1 != b.chrom1 or a.chrom2 != b.chrom2:
        return False
    if abs(a.pos1 - b.pos1) > distance or abs(a.pos2 - b.pos2) > distance:
        return False
    if a.svtype == "INS":
        lo, hi = sorted((a.svlen, b.svlen))
        if hi > 0 and lo / hi < ins_len_ratio:
            return False
    return True


def _or_supp(a, b):
    if a == "." or b == "." or len(a) != len(b):
        return a if b == "." else (b if a == "." else a)
    return "".join("1" if x == "1" or y == "1" else "0" for x, y in zip(a, b))


def _n_callers(j):
    return j.supp_vec.count("1") if j.supp_vec != "." else 0


def collapse(junctions, distance, order, ins_len_ratio=0.5):
    """Merge records describing the same junction (leader-based, never chains).

    The leader of each group is its record with the most supporting callers
    (ties: smaller pos1, pos2, then record id), so the reported coordinates are
    those of the best-supported call. SUPP_VEC is OR-ed, dosage is the maximum,
    record ids and loci are pooled. Output is sorted by contig order and position.
    """
    groups = {}
    for j in junctions:
        groups.setdefault((j.svtype, j.chrom1, j.chrom2), []).append(j)
    merged = []
    for _, js in sorted(groups.items(), key=lambda kv: (kv[0][0], order.get(kv[0][1], 1 << 30),
                                                        order.get(kv[0][2], 1 << 30))):
        # leaders first: best-supported records claim their neighbours
        rank = sorted(js, key=lambda j: (-_n_callers(j), j.pos1, j.pos2, j.record_ids[0]))
        by_pos = sorted(range(len(rank)), key=lambda i: (rank[i].pos1, i))
        pos_sorted = [rank[i].pos1 for i in by_pos]
        taken = [False] * len(rank)
        for i, lead in enumerate(rank):
            if taken[i]:
                continue
            taken[i] = True
            g = Junction(lead.chrom1, lead.pos1, lead.chrom2, lead.pos2, lead.svtype, lead.svlen,
                         lead.supp_vec, lead.dosage, list(lead.record_ids), set(lead.locus_ids))
            lo = bisect.bisect_left(pos_sorted, lead.pos1 - distance)
            hi = bisect.bisect_right(pos_sorted, lead.pos1 + distance)
            for k in sorted(by_pos[lo:hi]):          # rank order, so the result is order-independent
                if not taken[k] and same_junction(lead, rank[k], distance, ins_len_ratio):
                    taken[k] = True
                    o = rank[k]
                    g.supp_vec = _or_supp(g.supp_vec, o.supp_vec)
                    g.dosage = max(g.dosage, o.dosage)
                    g.record_ids += o.record_ids
                    g.locus_ids |= o.locus_ids
            g.record_ids = sorted(set(g.record_ids))
            merged.append(g)
    merged.sort(key=lambda j: (order.get(j.chrom1, 1 << 30), j.pos1, order.get(j.chrom2, 1 << 30),
                               j.pos2, j.svtype, j.svlen))
    return merged


# ─── CxSV locus membership ───────────────────────────────────────────────────
def canonical_id(vid):
    """Strip the _mate/_dest suffixes vcf_to_bed.py adds to breakend rows."""
    for suf in ("_mate", "_dest"):
        if vid.endswith(suf):
            return vid[: -len(suf)]
    return vid


def cxsv_loci_of_records(clusters_bed, summary_tsv):
    """{record id: {locus_id}} for records inside the sample's CxSV clusters.

    clusters.bed (bedtools cluster on merged_sv.bed: CHROM START END SVTYPE ID
    IS_INTER CLUSTER) says which cluster each record (or breakend) is in;
    cxsv_summary.tsv says which clusters are CxSVs. A locus id is
    "Chr:Start-End" exactly as in cxsv_master_table.tsv. A translocation with
    one breakend in each of two CxSV loci belongs to both.
    """
    cxsv = {}
    with open(summary_tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > ix["IsCxSV"] and f[ix["IsCxSV"]] == "YES":
                cxsv[f[ix["ClusterID"]]] = f"{f[ix['Chr']]}:{f[ix['Start']]}-{f[ix['End']]}"
    members = {}
    with open(clusters_bed) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 7 and f[6] in cxsv:
                members.setdefault(canonical_id(f[4]), set()).add(cxsv[f[6]])
    return members


def sample_junctions(vcf, order, collapse_distance, clusters_bed=None, summary_tsv=None,
                     ins_len_ratio=0.5):
    """The collapsed junction table of one sample, with CxSV locus membership."""
    js = parse_vcf(vcf, order)
    if clusters_bed and summary_tsv:
        members = cxsv_loci_of_records(clusters_bed, summary_tsv)
        for j in js:
            j.locus_ids = set(members.get(j.record_ids[0], ()))
    return collapse(js, collapse_distance, order, ins_len_ratio)


def write_table(path, junctions):
    with open(path, "w") as fh:
        fh.write("\t".join(JUNCTION_COLUMNS) + "\n")
        for j in junctions:
            r = j.as_row()
            fh.write("\t".join(str(r[c]) for c in JUNCTION_COLUMNS) + "\n")


def read_table(path):
    """Junction rows as dicts with int positions and a set of locus ids."""
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            for k in ("pos1", "pos2", "svlen", "dosage", "n_records"):
                if k in r:
                    r[k] = int(r[k])
            r["loci"] = set() if r.get("locus_id", ".") == "." else set(r["locus_id"].split(","))
            rows.append(r)
    return rows
