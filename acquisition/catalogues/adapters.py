"""Per-source adapters: one catalogue record -> zero or more junction sites.

Pure functions over already-parsed fields (no I/O), so they are unit-tested
in tests/unit/test_catalogue_junction_sites.py. convert_catalogue.py feeds
them from ``bcftools query`` output or from the HGSVC3 complex-SV BED.

Each adapter returns ``(sites, reason)``: a list of site dicts from
junction_sites.make_site, and, when the list is empty, a short reason string
that convert_catalogue.py counts into the manifest (``dropped_<reason>``).

Contig names are passed through junction_sites.normalise_contig. For a
GRCh38 source (gnomAD-SV) the GRCh38 names (chr1..chrY, chrM) are kept as
they are, and liftover.py maps both ends to hs1 afterwards.
"""
from junction_sites import (make_site, normalise_contig, normalise_svtype,
                            classify_sequence_allele, parse_bnd_alt)

MIN_SV_LEN = 50


def _int(v):
    if v in (None, "", "."):
        return None
    return int(float(v))


def _float(v):
    if v in (None, "", "."):
        return None
    return float(v)


def vcf_sv_sites(chrom, pos, ref, alt, svtype=None, end=None, svlen=None,
                 chr2=None, end2=None, inv_flag=False, ac=None, an=None,
                 af=None, pop_af=None, min_len=MIN_SV_LEN):
    """Generic biallelic VCF SV record -> sites.

    Handles symbolic ALTs (<DEL>, <DUP>, <INV>, <INS>, <INS:ME:ALU>, <BND>,
    <CTX>), breakend ALTs (N[chr2:100[) and sequence-resolved alleles.
    """
    c1 = normalise_contig(chrom)
    if c1 is None:
        return [], "non_hs1_contig"
    pos = int(pos)
    ac, an, af = _int(ac), _int(an), _float(af)
    if ac == 0:
        return [], "ac_zero"
    end, svlen = _int(end), _int(svlen)
    kw = dict(ac=ac, an=an, af=af, pop_af=pop_af)

    if "[" in alt or "]" in alt:
        mate = parse_bnd_alt(alt)
        c2 = normalise_contig(mate[0]) if mate else None
        if c2 is None:
            return [], "bnd_mate_unparsed"
        return [make_site(c1, pos, c2, mate[1], "BND", None, **kw)], None

    if alt.startswith("<"):
        t = normalise_svtype(svtype) or normalise_svtype(alt)
        if t is None:
            return [], "svtype_" + (svtype or alt.strip("<>")).replace(":", "_")
        if t == "BND":
            c2 = normalise_contig(chr2)
            p2 = _int(end2)
            if c2 is None or p2 is None:
                return [], "bnd_mate_missing"
            span = abs(svlen) if (svlen and svlen > 0 and c2 == c1) else None
            return [make_site(c1, pos, c2, p2, "BND", span, **kw)], None
        if t == "INS":
            ln = abs(svlen) if svlen is not None else None
            if ln is not None and ln < min_len:
                return [], "below_min_len"
            return [make_site(c1, pos, c1, pos, "INS", ln, **kw)], None
        # DEL / DUP / INV: pos2 = END
        if end is None or end == pos:
            if svlen is None:
                return [], "no_end"
            end = pos + abs(svlen)
        ln = abs(svlen) if svlen is not None else end - pos
        if ln < min_len:
            return [], "below_min_len"
        return [make_site(c1, pos, c1, end, t, ln, **kw)], None

    cls = classify_sequence_allele(pos, ref, alt, min_len=min_len,
                                   inv_flag=inv_flag)
    if cls is None:
        return [], "below_min_len"
    t, p1, p2, ln = cls
    return [make_site(c1, p1, c1, p2, t, ln, **kw)], None


# --------------------------------------------------------------- gnomAD-SV
GNOMAD_POPS = ["afr", "ami", "amr", "asj", "eas", "fin", "mid", "nfe", "rmi", "sas"]

# CPX_INTERVALS component prefix -> schema svtype
_GNOMAD_CPX_COMPONENT = {"DEL": "DEL", "DUP": "DUP", "dDUP": "DUP",
                         "INV": "INV", "INS": "INS", "dDEL": "DEL"}


def parse_gnomad_cpx_intervals(text):
    """'DUP_chr1:100-200,INV_chr1:150-300' -> [('DUP','chr1',100,200), ...]."""
    out = []
    if text in (None, "", "."):
        return out
    for tok in text.split(","):
        tok = tok.strip()
        if "_" not in tok or ":" not in tok:
            continue
        kind, loc = tok.split("_", 1)
        chrom, rng = loc.rsplit(":", 1)
        a, b = rng.split("-")
        out.append((kind, chrom, int(a), int(b)))
    return out


def gnomad_sites(chrom, pos, alt, filt, svtype, end, svlen, chr2, end2,
                 cpx_intervals, ac, an, af, pop_af, pass_only=True,
                 min_len=MIN_SV_LEN):
    """gnomAD-SV v4 record -> sites (in GRCh38 coordinates).

    * PASS only by default (``pass_only``); other FILTERs -> dropped_filter_*.
    * CNV (multiallelic copy-number sites) are not representable -> dropped.
    * CTX (reciprocal translocation) and BND -> BND to (CHR2, END2).
    * CPX -> one BND for the outer span (POS..END on one contig, svlen=SVLEN)
      plus one site per CPX_INTERVALS component (DEL/DUP/INV; dDUP -> DUP);
      all inherit the CPX record's AC/AN/AF.
    """
    if pass_only and filt not in ("PASS", "."):
        return [], "filter_" + filt.split(";")[0]
    st = (svtype or "").upper()
    if st in ("CNV", "MCNV"):
        return [], "svtype_CNV"
    if st == "CPX":
        c1 = normalise_contig(chrom)
        if c1 is None:
            return [], "non_hs1_contig"
        acv, anv, afv = _int(ac), _int(an), _float(af)
        if acv == 0:
            return [], "ac_zero"
        kw = dict(ac=acv, an=anv, af=afv, pop_af=pop_af)
        pos, endv = int(pos), _int(end) or int(pos)
        sites = [make_site(c1, pos, c1, endv, "BND", _int(svlen) or None, **kw)]
        for kind, cc, a, b in parse_gnomad_cpx_intervals(cpx_intervals):
            t = _GNOMAD_CPX_COMPONENT.get(kind)
            cn = normalise_contig(cc)
            if t is None or cn is None:
                continue
            if t == "INS":
                sites.append(make_site(cn, a, cn, a, "INS", b - a, **kw))
            elif b - a >= min_len:
                sites.append(make_site(cn, a, cn, b, t, b - a, **kw))
        return sites, None
    if st in ("CTX", "BND"):
        return vcf_sv_sites(chrom, pos, "N", "<BND>", "BND", end, svlen,
                            chr2=chr2, end2=end2, ac=ac, an=an, af=af,
                            pop_af=pop_af, min_len=min_len)
    return vcf_sv_sites(chrom, pos, "N", alt, svtype, end, svlen, ac=ac,
                        an=an, af=af, pop_af=pop_af, min_len=min_len)


# ------------------------------------------------------------------ HGSVC3
# STRUCT_REF segment token -> schema svtypes emitted for that segment
_HGSVC_SEGMENT = {
    "DEL": ("DEL",), "DUP": ("DUP",), "TRP": ("DUP",), "QUAD": ("DUP",),
    "MIXTRP": ("DUP",), "MIXQUAD": ("DUP",), "INV": ("INV",),
    "INVDUP": ("INV", "DUP"), "INVTRP": ("INV", "DUP"), "INVQUAD": ("INV", "DUP"),
    "NML": (),
}


def parse_struct_ref(text):
    """'DEL(4259):INV(10):NML(5)' -> [('DEL', 4259), ('INV', 10), ('NML', 5)]."""
    out = []
    for tok in (text or "").split(":"):
        tok = tok.strip()
        if not tok or "(" not in tok:
            continue
        name, n = tok.split("(", 1)
        out.append((name, int(n.rstrip(")"))))
    return out


def hgsvc3_cpx_sites(chrom, start0, end0, svlen, struct_ref, ac, an, af=None,
                     min_len=MIN_SV_LEN):
    """One HGSVC3 PAV complex-SV BED row -> sites.

    The BED start is 0-based, so as a 1-based coordinate it is the last
    unaffected base left of the event; the BED end is the 1-based last
    affected base. The outer span becomes one BND (chrom:start..chrom:end,
    svlen = SVLEN). When the STRUCT_REF segment lengths tile the span exactly
    (sum == end - start), each segment is also emitted as typed sites
    (DEL; DUP for DUP/TRP/QUAD/MIX*; INV; INV+DUP for INVDUP/INVTRP/INVQUAD;
    NML skipped). Returns (sites, reason, tiled: bool).
    """
    c1 = normalise_contig(chrom)
    if c1 is None:
        return [], "non_hs1_contig", False
    s, e = int(start0), int(end0)
    acv, anv, afv = _int(ac), _int(an), _float(af)
    if acv == 0:
        return [], "ac_zero", False
    kw = dict(ac=acv, an=anv, af=afv)
    sites = [make_site(c1, max(s, 1), c1, e, "BND", _int(svlen), **kw)]
    segs = parse_struct_ref(struct_ref)
    tiled = bool(segs) and sum(n for _, n in segs) == e - s
    if tiled:
        cur = s
        for name, n in segs:
            for t in _HGSVC_SEGMENT.get(name, ()):
                if n >= min_len:
                    sites.append(make_site(c1, max(cur, 1), c1, cur + n, t, n, **kw))
            cur += n
    return sites, None, tiled
