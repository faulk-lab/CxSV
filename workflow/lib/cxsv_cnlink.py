"""Copy-number-linked candidates of one genome (roadmap item 15, Phase 2; docs/CN_LINKAGE.md).
Loci-defining: fingerprinted. Phase 1 prototype: calibration/cn_linkage/cnlink.py.
Candidate structures, distance-free junction association and the structure model
(roadmap item 15, Phase 1; docs/CN_LINKAGE_PHASE1.md, docs/DUP_TRP_PLAN.md section 3).

Input per genome: abnormal copy-number segments (cnlib.segment_genome) and every
record of merged_sv.vcf, including the records discovery sets aside
(cxsv/long_span_records.tsv: DEL/DUP/INV > span_link_max; cxsv/single_caller_filtered.tsv:
svim_lowq).

  group_segments      abnormal segments on one contig within max_gap of each other form
                      a candidate; candidates joined by an attached translocation merge
  attach_records      a record joins a candidate when a breakend lies within tol of a
                      segment boundary (no distance clustering), and its type fits:
                        DEL  reciprocal overlap >= 0.5 with a loss segment of the candidate
                        DUP  reciprocal overlap >= 0.5 with a gain segment
                        INV  other end at a boundary or inside the candidate span
                        BND/TRA  any other end (it may link a second candidate)
                      INS and DEL/DUP < min_cn_record_bp are never attached
  classify            DUP-TRP/INV-DUP: a >= +2 segment adjacent to a +1 segment (nested
                      gain) AND an inverted junction (INV, or BND with ++ / -- sides);
                      other complex patterns: CN_nested_gain, CN_linked_gain,
                      CN_linked_loss, CN_linked_gain_loss
  is_complex          >= 1 attached junction and (>= 2 abnormal segments, or an
                      orientation-changing junction, or >= 2 junctions); a single gain or
                      loss explained by one DUP / DEL is a simple CNV, not a locus;
                      candidates with no attached junction are CN-only (reported, no locus)
  cn_public           segment-level panel frequency: every abnormal segment of the
                      candidate seen (same direction, reciprocal overlap >= 0.5) in at
                      least public_min_genomes (autosomes, 1) or public_min_genomes_sexchrom
                      (chrX / chrY, 2; panel-v1.3 D2) other individuals -> dropped as public.
                      Counts are distinct individuals. On chrX / chrY only individuals of the
                      sample's sex call count, and libraries excluded by sexchrom_discordant
                      count for no one.
  locus_span          a "locus" candidate becomes one CxSV locus per contig that carries
                      one of its segments: Start-End = the union of those segments and the
                      attached breakends on that contig (panel-v1.3 item g). Breakends on a
                      contig without a segment (a translocation's mate) stay records of the
                      locus but make no locus of their own; check_locus rejects a zero-length,
                      off-contig or out-of-segment locus.
  sexchrom_discordant an individual with > 1 panel library whose libraries disagree on
                      sex (different sex_call, or chrX / chrY ratios differing by more than
                      sexchrom_discordance_x / _y): all its libraries are left out of
                      chrX / chrY public counting (panel-v1.3 D2; NA19238: LRSC chrX 1.17 vs
                      ONT-1019 1.02). Recorded in pon_members.tsv (sexchrom_public).
"""
import re
from collections import defaultdict

CN_DEFAULTS = dict(
    max_gap=500_000,          # abnormal segments closer than this on a contig: one candidate
    tol=10_000,               # breakend to segment boundary (bp)
    min_cn_record_bp=10_000,  # DEL/DUP shorter than this are not CN records
    rec_overlap=0.5,          # DEL/DUP vs loss/gain segment, reciprocal
    public_overlap=0.5,       # segment vs another genome's segment, reciprocal
    public_min_genomes=1,     # other individuals needed to call an autosomal segment public
    public_min_genomes_sexchrom=2,   # ... a chrX / chrY segment (sex-matched; panel-v1.3 D2)
    sexchrom_discordance_x=0.10,     # chrX ratio difference between two libraries of one
                                     # individual above which both leave chrX/chrY counting
    sexchrom_discordance_y=0.10,     # the same for the chrY ratio
    adjacent_bp=20_000,       # nested-gain adjacency
    span_link_max=500_000,    # long records (calling.yaml cxsv_params.span_link_max)
    inv_other_slack=100_000,  # INV with one end at a boundary: other end within the candidate
                              # span +/- this (the outer DUP of a DUP-TRP can be below detection)
)
INV_ALT = re.compile(r"^[ACGTNacgtn]*\][^\]]+\]$|^\[[^\[]+\[[ACGTNacgtn]*$")   # N]c:p] (++) or [c:p[N (--)
BND_ALT = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")


def info_dict(s):
    out = {}
    for kv in s.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
        else:
            out[kv] = True
    return out


def read_records(vcf):
    """merged_sv.vcf -> list of dicts: id, chrom, pos, chr2, end, svtype, svlen, supp, qual, inverted, line."""
    recs = []
    for line in open(vcf):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 8:
            continue
        info = info_dict(f[7])
        st = str(info.get("SVTYPE", "")).upper()
        c, p = f[0], int(f[1])
        c2, e = info.get("CHR2", c), info.get("END")
        m = BND_ALT.search(f[4])
        if st in ("BND", "TRA") and m:
            c2, e = m.group(1), m.group(2)
        try:
            e = int(e) if e not in (None, True) else p
        except ValueError:
            e = p
        try:
            q = float(f[5])
        except ValueError:
            q = None
        sl = info.get("SVLEN")
        try:
            sl = abs(int(sl))
        except (TypeError, ValueError):
            sl = abs(e - p) if c2 == c else 0
        inverted = (st == "INV") or (st in ("BND", "TRA") and bool(INV_ALT.match(f[4])))
        recs.append(dict(id=f[2], chrom=c, pos=p, chr2=c2, end=e, svtype=st, svlen=sl,
                         supp=str(info.get("SUPP_VEC", "")), qual=q, inverted=inverted, line=line))
    return recs


def _recip(a0, a1, b0, b1):
    ov = min(a1, b1) - max(a0, b0)
    if ov <= 0:
        return 0.0
    return min(ov / max(1, a1 - a0), ov / max(1, b1 - b0))


def group_segments(segs, P=CN_DEFAULTS):
    by = defaultdict(list)
    for s in segs:
        by[s["chrom"]].append(s)
    groups = []
    for c, ss in by.items():
        ss.sort(key=lambda s: s["start"])
        cur = [ss[0]]
        for s in ss[1:]:
            if s["start"] - max(x["end"] for x in cur) <= P["max_gap"]:
                cur.append(s)
            else:
                groups.append(cur)
                cur = [s]
        groups.append(cur)
    return [dict(id=i, segs=g, chrom=g[0]["chrom"], start=min(s["start"] for s in g),
                 end=max(s["end"] for s in g), records=[]) for i, g in enumerate(groups)]


def _boundaries(g):
    """(contig, position, masked bp before, masked bp after) of every segment boundary."""
    return [(s["chrom"], s[k], int(s.get(f"{k}_mask_before", 0) or 0), int(s.get(f"{k}_mask_after", 0) or 0))
            for s in g["segs"] for k in ("start", "end")]


def attach_records(groups, recs, P=CN_DEFAULTS):
    """Attach records to candidate groups; returns {group id: [record]} and translocation links."""
    bidx = defaultdict(list)
    for g in groups:
        for c, p, mb, ma in _boundaries(g):
            bidx[c].append((p - mb - P["tol"], p + ma + P["tol"], g["id"]))
    for v in bidx.values():
        v.sort()
    gid = {g["id"]: g for g in groups}

    def at_boundary(c, p):
        """groups with a boundary within tol of p (widened over adjacent masked bins)."""
        return {i for lo, hi, i in bidx.get(c, ()) if lo <= p <= hi}

    links = []
    for r in recs:
        st = r["svtype"]
        if st == "INS":
            continue
        h1, h2 = at_boundary(r["chrom"], r["pos"]), at_boundary(r["chr2"], r["end"])
        if not h1 and not h2:
            continue
        if st in ("DEL", "DUP"):
            if r["svlen"] < P["min_cn_record_bp"] or r["chr2"] != r["chrom"]:
                continue
            want = "loss" if st == "DEL" else "gain"
            for i in h1 | h2:
                g = gid[i]
                if any(_recip(r["pos"], r["end"], s["start"], s["end"]) >= P["rec_overlap"]
                       and (s["state"] < s["ploidy"]) == (want == "loss") for s in g["segs"]):
                    g["records"].append(r)
        elif st == "INV":
            for i in h1 | h2:
                g = gid[i]
                other = r["end"] if i in h1 else r["pos"]
                if (i in h1 and i in h2) or (g["start"] - P["inv_other_slack"] <= other <= g["end"] + P["inv_other_slack"]):
                    g["records"].append(r)
        else:                                       # BND / TRA
            for i in h1 | h2:
                gid[i]["records"].append(r)
            if h1 and h2 and not (h1 & h2):
                links.append((min(h1), min(h2)))
    return links


def merge_linked(groups, links):
    """Union candidates joined by a translocation attached at both ends."""
    parent = {g["id"]: g["id"] for g in groups}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in links:
        parent[find(a)] = find(b)
    out = defaultdict(lambda: dict(segs=[], records=[], parts=[]))
    for g in groups:
        o = out[find(g["id"])]
        o["segs"] += g["segs"]
        o["parts"].append((g["chrom"], g["start"], g["end"]))
        seen = {r["id"] for r in o["records"]}
        o["records"] += [r for r in g["records"] if r["id"] not in seen]
    res = []
    for k, o in out.items():
        o["id"] = k
        o["chrom"], o["start"], o["end"] = o["parts"][0]
        res.append(o)
    return res


def classify(c, P=CN_DEFAULTS):
    segs = sorted(c["segs"], key=lambda s: (s["chrom"], s["start"]))
    delta = [s["state"] - s["ploidy"] for s in segs]
    gains, losses = any(d > 0 for d in delta), any(d < 0 for d in delta)
    nested = False
    for s in segs:
        if s["state"] - s["ploidy"] >= 2:
            for t in segs:
                if (t is not s and t["chrom"] == s["chrom"] and t["state"] - t["ploidy"] == 1
                        and (abs(t["end"] - s["start"]) <= P["adjacent_bp"] or abs(s["end"] - t["start"]) <= P["adjacent_bp"])):
                    nested = True
    inverted = any(r["inverted"] for r in c["records"])
    if nested and inverted:
        cls = "DUP-TRP/INV-DUP"
    elif nested:
        cls = "CN_nested_gain"
    elif gains and losses:
        cls = "CN_linked_gain_loss"
    elif gains:
        cls = "CN_linked_gain"
    else:
        cls = "CN_linked_loss"
    return cls, nested, inverted


def is_complex(c):
    recs = c["records"]
    if not recs:
        return False, "cn_only"
    orient = any(r["svtype"] in ("INV", "BND", "TRA") for r in recs)
    if len(c["segs"]) >= 2 or orient or len(recs) >= 2:
        return True, "complex"
    return False, "simple_cnv"


SEX_CHROMS = ("chrX", "chrY")


def _ratio(v):
    try:
        v = float(v)
    except (TypeError, ValueError):
        return None
    return v if v == v else None              # NaN -> None


def sexchrom_discordant(libs, P=CN_DEFAULTS):
    """{key: reason} of panel libraries left out of chrX / chrY public counting.

    libs: iterable of dicts with key, individual_id, sex_call, x_ratio, y_ratio (the
    columns of pon_members.tsv). An individual with more than one library is discordant
    when two of its libraries have different sex calls, or chrX ratios differing by more
    than sexchrom_discordance_x, or chrY ratios by more than sexchrom_discordance_y (a
    missing ratio is not compared). Then every library of that individual is excluded:
    which one is wrong cannot be told from depth alone."""
    by = defaultdict(list)
    for r in libs:
        by[str(r["individual_id"])].append(r)
    out = {}
    for ind, rs in sorted(by.items()):
        if len(rs) < 2:
            continue
        why = []
        for i, a in enumerate(rs):
            for b in rs[i + 1:]:
                ka, kb = a["key"], b["key"]
                ca, cb = (a.get("sex_call") or "").strip(), (b.get("sex_call") or "").strip()
                if ca != cb:
                    why.append(f"sex call {ca} ({ka}) vs {cb} ({kb})")
                for ax, tol in (("x", P["sexchrom_discordance_x"]), ("y", P["sexchrom_discordance_y"])):
                    va, vb = _ratio(a.get(f"{ax}_ratio")), _ratio(b.get(f"{ax}_ratio"))
                    if va is not None and vb is not None and abs(va - vb) > tol:
                        why.append(f"chr{ax.upper()} ratio {va:.3f} ({ka}) vs {vb:.3f} ({kb}), "
                                   f"|diff| {abs(va - vb):.3f} > {tol}")
        if why:
            reason = "discordant libraries of " + ind + ": " + "; ".join(why)
            for r in rs:
                out[r["key"]] = reason
    return out


def min_genomes(chrom, P=CN_DEFAULTS):
    """Other individuals needed to call a segment on `chrom` public."""
    if chrom in SEX_CHROMS:
        return P.get("public_min_genomes_sexchrom", P["public_min_genomes"])
    return P["public_min_genomes"]


def cn_public(c, others, P=CN_DEFAULTS, sex=None):
    """(public?, per-segment counts of other individuals carrying it).
    others: list of (genome, individual, [segments]), (..., sex_call) or
    (..., sex_call, sexchrom_ok) of the reference libraries. Counts are distinct
    individuals (two libraries of one individual count once). A segment on chrX / chrY is
    compared only with libraries of the same sex call as `sex` (the sample's ploidy sex;
    panel-v1.2-r2) and never with a library whose sexchrom_ok is False (sexchrom_discordant;
    panel-v1.3 D2); it is public at >= public_min_genomes_sexchrom individuals, an
    autosomal segment at >= public_min_genomes (every sex). sex None (or others without a
    sex): no sex matching."""
    counts = []
    for s in c["segs"]:
        sign = s["state"] > s["ploidy"]
        sexchrom = s["chrom"] in SEX_CHROMS
        inds = set()
        for o in others:
            g, ind, ss = o[0], o[1], o[2]
            if sexchrom:
                if sex and len(o) > 3 and o[3] != sex:
                    continue
                if len(o) > 4 and not o[4]:
                    continue
            if ind in inds:
                continue
            if any(t["chrom"] == s["chrom"] and (t["state"] > t["ploidy"]) == sign
                   and _recip(s["start"], s["end"], t["start"], t["end"]) >= P["public_overlap"] for t in ss):
                inds.add(ind)
        counts.append(len(inds))
    pub = all(n >= min_genomes(s["chrom"], P) for n, s in zip(counts, c["segs"])) if counts else False
    return pub, counts


def segment_contigs(c):
    """Contigs carrying at least one abnormal segment of candidate c (in segment order)."""
    out = []
    for s in sorted(c["segs"], key=lambda s: (s["chrom"], s["start"])):
        if s["chrom"] not in out:
            out.append(s["chrom"])
    return out


def locus_span(c, chrom, rows=(), length=None):
    """(start, end) of the copy-number locus of candidate c on `chrom`: the union of its
    segments there and the breakend intervals `rows` ((start, end) of the records' BED
    rows on that contig, vcf_to_bed.py). panel-v1.2 took the span from the rows alone, so
    a loss whose only junction was a telomeric translocation end became chr17:0-2
    (BOC_5883, 17p13.3). With `length` (the contig length) the end is clamped to it: segments
    are built on 1 kb bins, so one that reaches the last, partial bin ends past the contig
    (HG03667 chr2, design check panel-v1.3)."""
    segs = [s for s in c["segs"] if s["chrom"] == chrom]
    if not segs:
        raise ValueError(f"candidate {c.get('id')} has no copy-number segment on {chrom}")
    lo = min([s["start"] for s in segs] + [int(r[0]) for r in rows])
    hi = max([s["end"] for s in segs] + [int(r[1]) for r in rows])
    if length is not None:
        hi = min(hi, int(length))
    return max(0, lo), hi


def check_locus(chrom, start, end, c, length=None):
    """Raise ValueError unless chrom:start-end is a well-formed copy-number locus of
    candidate c: end > start >= 0, inside the contig (length, when known), and overlapping
    at least one of c's segments on chrom by >= 1 bp."""
    start, end = int(start), int(end)
    if start < 0 or end <= start:
        raise ValueError(f"malformed copy-number locus {chrom}:{start}-{end} (zero or negative length)")
    if length is not None and end > length:
        raise ValueError(f"copy-number locus {chrom}:{start}-{end} extends past the contig end ({length})")
    segs = [s for s in c["segs"] if s["chrom"] == chrom]
    if not any(min(end, s["end"]) > max(start, s["start"]) for s in segs):
        raise ValueError(f"copy-number locus {chrom}:{start}-{end} overlaps none of candidate "
                         f"{c.get('id')}'s segments ({len(segs)} on {chrom})")
    return start, end


def long_record(r, P=CN_DEFAULTS):
    return r["svtype"] in ("DEL", "DUP", "INV") and r["chr2"] == r["chrom"] and abs(r["end"] - r["pos"]) > P["span_link_max"]


def candidates(segs, recs, others=None, P=CN_DEFAULTS, single_filtered=(), sex=None):
    """Every candidate of one genome with its attached records, class and decision.
    sex: the genome's ploidy sex, for sex-matched panel frequency on chrX / chrY."""
    groups = group_segments(segs, P)
    links = attach_records(groups, recs, P)
    cands = merge_linked(groups, links)
    out = []
    for c in cands:
        cls, nested, inverted = classify(c, P)
        cx, why = is_complex(c)
        pub, counts = cn_public(c, others or [], P, sex=sex)
        decision = "locus" if cx and not pub else ("public_cn" if cx else why)
        n_long = sum(long_record(r, P) for r in c["records"])
        n_single = sum(r["id"] in single_filtered for r in c["records"])
        out.append(dict(c, cls=cls, nested=nested, inverted=inverted, decision=decision,
                        public_counts=counts, n_long=n_long, n_single_filtered=n_single,
                        span=sum(e - s for _, s, e in c["parts"])))
    return out
