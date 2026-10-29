"""Shared junction-site schema for external SV catalogues (roadmap item 4).

Every catalogue is converted to one tab-separated table, one row per SV site:

    chrom1 pos1 chrom2 pos2 svtype svlen ac an af source version [af_<group> ...]

Conventions (agreed for both item-4 tracks):

* 1-based positions, UCSC hs1 contig names (chr1..chr22, chrX, chrY, chrM).
* svtype is one of DEL, DUP, INV, INS, BND. Translocations are BND.
* DEL/DUP/INV: chrom2 = chrom1 and pos2 = END (VCF END: last affected base).
  INS: pos2 = pos1 (the anchor base) and svlen = inserted length.
  BND: (chrom2, pos2) is the mate breakend.
* Each site is written once, with (chrom1, pos1) <= (chrom2, pos2) in hs1
  .fai contig order, then position.
* ac / an / af are '.' when the source does not provide them.

Choices this module makes where the shared schema is silent (documented in
acquisition/catalogues/README.md and docs/EXTERNAL_CATALOGUES.md):

* svlen is UNSIGNED (absolute length; DEL svlen > 0). For BND it is the
  source's SVLEN when that is a positive span on one chromosome (e.g. the
  outer span of a complex SV), otherwise '.'.
* af is taken from the source when it is given, otherwise ac/an when both
  are known, otherwise '.'.
* Rows are sorted by (chrom1, pos1, chrom2, pos2, svtype, svlen) in hs1 .fai
  contig order, so the bgzipped table can be tabix-indexed on chrom1/pos1.

Standard library only.
"""
import csv
import gzip
import io
import os
from collections import Counter

SCHEMA_VERSION = "1"
CORE_COLUMNS = ["chrom1", "pos1", "chrom2", "pos2", "svtype", "svlen",
                "ac", "an", "af", "source", "version"]
SVTYPES = ("DEL", "DUP", "INV", "INS", "BND")
MISSING = "."
HS1_CONTIGS = tuple([f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"])

# Source-specific SV type spellings -> schema svtype. Types absent from this
# map are not representable as a junction site and are counted as dropped.
SVTYPE_MAP = {
    "DEL": "DEL", "DUP": "DUP", "INV": "INV", "INS": "INS", "BND": "BND",
    "TRA": "BND", "CTX": "BND",
    "DUP:TANDEM": "DUP", "DEL:ME": "DEL",
    "INS:ME": "INS", "INS:ME:ALU": "INS", "INS:ME:LINE1": "INS",
    "INS:ME:SVA": "INS", "INS:ME:HERV": "INS", "MEI": "INS",
}


def normalise_contig(name):
    """Map common contig spellings to hs1 names; return None if not hs1."""
    if name is None:
        return None
    n = name
    for prefix in ("CHM13#0#", "GRCh38#0#"):
        if n.startswith(prefix):
            n = n[len(prefix):]
    if not n.startswith("chr"):
        n = "chr" + n
    if n in ("chrMT", "chrM"):
        n = "chrM"
    if n == "chr23":
        n = "chrX"
    if n == "chr24":
        n = "chrY"
    return n if n in HS1_CONTIGS else None


def normalise_svtype(svtype):
    if svtype is None:
        return None
    s = svtype.strip().upper().strip("<>")
    if s in SVTYPE_MAP:
        return SVTYPE_MAP[s]
    if s.startswith("INS:"):
        return "INS"
    if s.startswith("DUP:"):
        return "DUP"
    if s.startswith("DEL:"):
        return "DEL"
    return None


class ContigOrder:
    """Contig rank in hs1 .fai order (the order of the reference index)."""

    def __init__(self, names):
        self.names = list(names)
        self.rank = {n: i for i, n in enumerate(self.names)}

    @classmethod
    def from_fai(cls, path):
        with open(path) as fh:
            return cls([line.split("\t", 1)[0] for line in fh if line.strip()])

    def key(self, chrom, pos):
        return (self.rank[chrom], int(pos))

    def __contains__(self, chrom):
        return chrom in self.rank


def make_site(chrom1, pos1, chrom2, pos2, svtype, svlen=None,
              ac=None, an=None, af=None, pop_af=None):
    """Return a site dict (not yet oriented). Missing values become None."""
    return {"chrom1": chrom1, "pos1": int(pos1),
            "chrom2": chrom2, "pos2": int(pos2), "svtype": svtype,
            "svlen": svlen, "ac": ac, "an": an, "af": af,
            "pop_af": dict(pop_af or {})}


def orient(site, order):
    """Swap the two ends if needed so that (chrom1,pos1) <= (chrom2,pos2)."""
    k1 = order.key(site["chrom1"], site["pos1"])
    k2 = order.key(site["chrom2"], site["pos2"])
    if k2 < k1:
        site = dict(site)
        site["chrom1"], site["chrom2"] = site["chrom2"], site["chrom1"]
        site["pos1"], site["pos2"] = site["pos2"], site["pos1"]
    return site


def validate(site, order):
    """Raise ValueError if a site breaks the schema rules."""
    if site["svtype"] not in SVTYPES:
        raise ValueError(f"bad svtype {site['svtype']!r}")
    for c in (site["chrom1"], site["chrom2"]):
        if c not in order:
            raise ValueError(f"contig {c!r} not in hs1 order")
    if site["pos1"] < 1 or site["pos2"] < 1:
        raise ValueError("positions must be 1-based (>= 1)")
    if order.key(site["chrom1"], site["pos1"]) > order.key(site["chrom2"], site["pos2"]):
        raise ValueError("site not oriented")
    if site["svtype"] != "BND" and site["chrom1"] != site["chrom2"]:
        raise ValueError(f"{site['svtype']} with two contigs")
    if site["svtype"] == "INS" and site["pos1"] != site["pos2"]:
        raise ValueError("INS must have pos2 == pos1")


def fmt_int(v):
    return MISSING if v is None or v == "" else str(int(v))


def fmt_af(v):
    if v is None or v == "":
        return MISSING
    return f"{float(v):.6g}"


def resolve_af(ac, an, af):
    """Source AF if given, else AC/AN when both known, else None."""
    if af not in (None, ""):
        return float(af)
    if ac not in (None, "") and an not in (None, "") and int(an) > 0:
        return int(ac) / int(an)
    return None


def site_row(site, source, version, af_groups):
    row = [site["chrom1"], str(site["pos1"]), site["chrom2"], str(site["pos2"]),
           site["svtype"], fmt_int(site["svlen"]), fmt_int(site["ac"]),
           fmt_int(site["an"]), fmt_af(resolve_af(site["ac"], site["an"], site["af"])),
           source, version]
    row += [fmt_af(site["pop_af"].get(g)) for g in af_groups]
    return row


def sort_key(order):
    def k(s):
        return (order.rank[s["chrom1"]], s["pos1"], order.rank[s["chrom2"]],
                s["pos2"], s["svtype"], -1 if s["svlen"] is None else int(s["svlen"]))
    return k


def write_sites(fh, sites, source, version, af_groups, order):
    """Validate, orient, sort and write sites as TSV to an open text handle.

    Returns a Counter of sites by svtype.
    """
    oriented = [orient(s, order) for s in sites]
    for s in oriented:
        validate(s, order)
    oriented.sort(key=sort_key(order))
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(CORE_COLUMNS + [f"af_{g}" for g in af_groups])
    counts = Counter()
    for s in oriented:
        w.writerow(site_row(s, source, version, af_groups))
        counts[s["svtype"]] += 1
    return counts


def read_sites(path):
    """Read a junction-site table (plain or gzip/bgzip) as a list of dicts."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def classify_sequence_allele(pos, ref, alt, min_len=50, inv_flag=False):
    """Classify a sequence-resolved VCF allele as a junction site.

    Returns (svtype, pos1, pos2, svlen) or None when the allele is not an
    SV of at least ``min_len`` bp net length (and is not flagged inversion).

    * Net loss  >= min_len -> DEL, pos1 = POS + shared-prefix - 1 (the last
      retained base; POS for a padded allele), pos2 = POS + len(REF) - 1 (END).
    * Net gain  >= min_len -> INS at the anchor base, svlen = net gain.
    * ``inv_flag`` and len(REF) >= min_len -> INV over the REF span.
    """
    if alt.startswith("<") or "[" in alt or "]" in alt or alt in ("*", "."):
        return None
    lr, la = len(ref), len(alt)
    p = 0
    while p < min(lr, la) and ref[p].upper() == alt[p].upper():
        p += 1
    anchor = pos + p - 1 if p >= 1 else pos
    d = la - lr
    if inv_flag and lr >= min_len and abs(d) < min_len:
        return ("INV", pos, pos + lr - 1, lr)
    if d <= -min_len:
        return ("DEL", anchor, pos + lr - 1, -d)
    if d >= min_len:
        return ("INS", anchor, anchor, d)
    return None


def parse_bnd_alt(alt):
    """Mate (chrom, pos) from a VCF breakend ALT such as N[chr2:100[ ."""
    for br in ("[", "]"):
        if br in alt:
            inner = alt.split(br)[1]
            chrom, pos = inner.rsplit(":", 1)
            return chrom, int(pos)
    return None


def open_text_out(path):
    """Text handle writing to ``path`` atomically (temp file + rename)."""
    tmp = path + ".tmp"
    return tmp, io.open(tmp, "w", newline="")


def finish_atomic(tmp, path):
    os.replace(tmp, path)
