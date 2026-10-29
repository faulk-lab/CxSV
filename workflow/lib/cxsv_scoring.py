"""Configurable, explainable CxSV scoring (roadmap item 5; docs/SCORING.md).

An explicit points system modelled on the ACMG/ClinGen CNV scoring framework
(Riggs et al. 2020): every feature contributes signed points, the points add
up to a total, and the total maps to a tier. Weights and cut-offs come from
config/scoring.yaml; nothing is learned. Each feature function returns a
`Feature(name, group, value, points, rule)`, so a report can show the raw value,
the points and the rule that produced them.

Feature groups
--------------
privacy        panel carriers weighted by panel callability; external catalogues
               (max AF per catalogue, the most common one counts, never summed)
genes          genes broken by a breakend or wholly lost / gained (copy number
               from read depth), weighted by ClinGen HI / TS scores, gnomAD
               LOEUF as fallback; exons; gene count; no gene at all
event          complexity class and size
evidence       caller agreement (SUPP_VEC) and supporting reads; repeat /
               segmental-duplication context as an artefact penalty
chromothripsis the PI's criteria (germline-chromothripsis-pipeline-criteria):
               1 interleaving, 2 breakpoint-count tier, 5 copy-number
               oscillation (supporting, never a gate; sex-aware, PAR excluded),
               7 breakpoint clustering vs a chromosome-intrinsic simulated null
               (descriptive), 8 fragment-join balance / fold-back (BFB) as a
               modifier of the chromothripsis points only
phenotype      optional: HPO match of the patient's terms to the affected genes;
               optional OMIM (licence key, off by default)
inheritance    optional: de novo vs inherited when parents are sequenced

All thresholds are operational, cohort-specific choices, not validated cutoffs.
Standard library + numpy.
"""
from __future__ import annotations

import bisect
import hashlib
import math
import re
from collections import namedtuple



Feature = namedtuple("Feature", "name group value points rule")

FEATURE_ORDER = ("panel_privacy", "catalogue_frequency", "gene_dosage", "exons", "gene_count", "no_genes",
                 "complexity_class", "size", "caller_agreement", "read_support", "artefact_context",
                 "ct_interleaving", "ct_breakpoint_tier", "ct_cn_oscillation", "ct_clustering",
                 "ct_bfb_modifier", "hpo_match", "omim", "inheritance")
CT_FEATURES = ("ct_interleaving", "ct_breakpoint_tier", "ct_cn_oscillation", "ct_clustering")


# ─── small helpers ───────────────────────────────────────────────────────────
def r2(x):
    return round(float(x) + 0.0, 4)


def fnum(x, nd=3):
    return "NA" if x is None else (f"{x:.{nd}g}" if isinstance(x, float) else str(x))


def tier_of(total, tiers):
    """tiers: list of {name, min} sorted by min descending (config)."""
    for t in sorted(tiers, key=lambda t: -t["min"]):
        if total >= t["min"] - 1e-9:
            return t["name"]
    return tiers[-1]["name"] if tiers else "NA"


class IntervalIndex:
    """Sorted intervals per chromosome (0-based half-open BED) with labels."""

    def __init__(self, rows):
        self.by = {}
        for chrom, s, e, label in rows:
            self.by.setdefault(chrom, []).append((s, e, label))
        self.starts, self.maxlen = {}, {}
        for c, v in self.by.items():
            v.sort()
            self.starts[c] = [x[0] for x in v]
            self.maxlen[c] = max(x[1] - x[0] for x in v)

    def overlapping(self, chrom, s, e):
        """Intervals overlapping [s, e) (0-based half-open)."""
        if chrom not in self.by:
            return []
        v, st = self.by[chrom], self.starts[chrom]
        i = bisect.bisect_left(st, s - self.maxlen[chrom])
        out = []
        while i < len(v) and v[i][0] < e:
            if v[i][1] > s:
                out.append(v[i])
            i += 1
        return out

    def covered_bp(self, chrom, s, e):
        iv = sorted((max(a, s), min(b, e)) for a, b, _ in self.overlapping(chrom, s, e))
        tot, cur_s, cur_e = 0, None, None
        for a, b in iv:
            if cur_e is None or a > cur_e:
                if cur_e is not None:
                    tot += cur_e - cur_s
                cur_s, cur_e = a, b
            else:
                cur_e = max(cur_e, b)
        if cur_e is not None:
            tot += cur_e - cur_s
        return tot


def read_bed(path, label_col=3, split_labels=False):
    rows = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            lab = f[label_col] if len(f) > label_col else "."
            if split_labels:
                for g in lab.split(","):
                    rows.append((f[0], int(f[1]), int(f[2]), g))
            else:
                rows.append((f[0], int(f[1]), int(f[2]), lab))
    return rows


# ─── orientation of junctions (criterion 8) ─────────────────────────────────
_ALT_FORMS = (
    (re.compile(r"^[A-Za-z.]+\[([^:\[\]]+):(\d+)\[$"), "+", "-"),   # t[p[
    (re.compile(r"^[A-Za-z.]+\]([^:\[\]]+):(\d+)\]$"), "+", "+"),   # t]p]
    (re.compile(r"^\]([^:\[\]]+):(\d+)\][A-Za-z.]+$"), "-", "+"),   # ]p]t
    (re.compile(r"^\[([^:\[\]]+):(\d+)\[[A-Za-z.]+$"), "-", "-"),   # [p[t
)


def record_orientation(chrom, pos, alt, svtype):
    """Join orientation of one VCF record as a set of (side at lower end, side at
    higher end) pairs: '+' = the segment ending at the break is kept (joined at its
    3' end), '-' = the segment starting at it. DEL-like (+,-); DUP-like (-,+);
    inversion-type joins (+,+) and (-,-) (ShatterSeek's h2hINV / t2tINV). A
    symbolic or sequence-resolved INV stands for both of its junctions; INS and
    unresolvable records: empty set."""
    t = svtype.upper()
    for rx, s_pos, s_mate in _ALT_FORMS:
        m = rx.match(alt)
        if m:
            c2, p2 = m.group(1), int(m.group(2))
            if c2 == chrom and p2 < pos:
                return {(s_mate, s_pos)}
            return {(s_pos, s_mate)}
    if t == "DEL":
        return {("+", "-")}
    if t in ("DUP", "DUP:TANDEM"):
        return {("-", "+")}
    if t == "INV":
        return {("+", "+"), ("-", "-")}
    return set()


JOIN_NAMES = {("+", "-"): "DEL-like", ("-", "+"): "DUP-like", ("+", "+"): "h2hINV", ("-", "-"): "t2tINV"}


def join_balance(junction_orients):
    """Counts of the four ShatterSeek fragment-join types over intrachromosomal
    junctions, and the chi-square goodness-of-fit p-value against equal
    proportions (df = 3; NA below 4 junctions)."""
    counts = {v: 0 for v in JOIN_NAMES.values()}
    for o in junction_orients:
        for pair in o:
            if pair in JOIN_NAMES:
                counts[JOIN_NAMES[pair]] += 1
    n = sum(counts.values())
    if n < 4:
        return counts, None
    e = n / 4.0
    x = sum((v - e) ** 2 / e for v in counts.values())
    p = math.erfc(math.sqrt(x / 2)) + math.sqrt(2 * x / math.pi) * math.exp(-x / 2)   # chi2 sf, df=3
    return counts, min(1.0, p)


def foldbacks(junctions, orients, tol, max_span, cn_gain_intervals):
    """Fold-back junctions (Xing et al. 2019, as adapted in the PI's criterion 8): an
    inversion-orientation junction (+,+ or -,-) with no reciprocal junction (the
    opposite orientation with both ends within tol), ends < max_span apart; with_cn
    = those whose span overlaps a copy-number gain."""
    out = []
    for j, o in zip(junctions, orients):
        if j.chrom1 != j.chrom2 or j.pos2 - j.pos1 >= max_span:
            continue
        inv = {p for p in o if p in (("+", "+"), ("-", "-"))}
        if len(inv) != 1:
            continue
        opp = ("-", "-") if inv == {("+", "+")} else ("+", "+")
        reciprocal = any(k is not j and opp in ok and abs(k.pos1 - j.pos1) <= tol and abs(k.pos2 - j.pos2) <= tol
                         for k, ok in zip(junctions, orients))
        if reciprocal:
            continue
        gain = any(c == j.chrom1 and s < j.pos2 and e > j.pos1 for c, s, e in cn_gain_intervals)
        out.append((j, gain))
    return out


# ─── copy number from read depth (criterion 5) ──────────────────────────────
def ploidy(chrom, pos, sex, par):
    """Copies of (chrom, pos): autosomes 2; chrX 1 in males outside the PAR, chrY
    1 in males (0 in females); PAR positions count 2. Unknown sex: chrX 2."""
    c = chrom.lower().removeprefix("chr")
    in_par = any(pc == chrom and s <= pos < e for pc, s, e in par)
    if c == "x":
        return 2 if (in_par or sex != "male") else 1
    if c == "y":
        return 0 if sex == "female" else (2 if in_par else 1)
    if c in ("m", "mt"):
        return 1
    return 2


def infer_sex(flank_depths, par, ratio_male=0.75, ratio_female=0.85):
    """('male'|'female'|'unknown', chrX/autosome depth ratio) from locus flank depths
    [(chrom, start, end, mean_depth)]."""
    auto = sorted(d for c, s, e, d in flank_depths if c.lower().removeprefix("chr").isdigit() and d > 0)
    xs = sorted(d for c, s, e, d in flank_depths if c == "chrX" and d > 0
                and not any(pc == c and ps <= s < pe for pc, ps, pe in par))
    if not auto or not xs:
        return "unknown", None
    ra = xs[len(xs) // 2] / auto[len(auto) // 2]
    if ra < ratio_male:
        return "male", ra
    if ra > ratio_female:
        return "female", ra
    return "unknown", ra


def cn_segments(depth_rows, sex, par, p):
    """Copy-number states of one locus's segments.

    depth_rows: [(chrom, start, end, kind, mean_depth)], kind flank_left /
    segment / flank_right. Baseline = length-weighted mean of the flanks on the
    same chromosome. Returns [(chrom, start, end, cn, state, ploidy, changed)] for
    segments >= min_segment_bp, or None when no usable baseline."""
    out = []
    for chrom in sorted({r[0] for r in depth_rows}):
        rows = [r for r in depth_rows if r[0] == chrom]
        fl = [(r[2] - r[1], r[4]) for r in rows if r[3].startswith("flank") and r[2] > r[1]]
        if not fl or sum(w for w, _ in fl) == 0:
            continue
        base = sum(w * d for w, d in fl) / sum(w for w, _ in fl)
        if base < p["min_baseline_depth"]:
            continue
        for c, s, e, kind, d in sorted(rows, key=lambda r: r[1]):
            if kind != "segment" or e - s < p["min_segment_bp"]:
                continue
            pl = ploidy(c, (s + e) // 2, sex, par)
            # the flanks carry the locus chromosome's ploidy; chrX flanks of a male are 1 copy
            base_pl = ploidy(c, s, sex, par) or 1
            cn = d / base * base_pl
            state = max(0, int(round(cn)))
            changed = abs(cn - pl) >= p["min_cn_delta"]
            out.append((c, s, e, cn, state if changed else pl, pl, changed))
    return out or None


def oscillation(segments):
    """Longest run of consecutive segments alternating between exactly two
    copy-number states (A B A B ...) that includes a changed state; and the
    number of distinct states."""
    if not segments:
        return 0, 0
    states = [s[4] for s in segments]
    best, cur = 1, 1
    for i in range(1, len(states)):
        if states[i] != states[i - 1] and (cur == 1 or states[i] == states[i - 2]):
            cur += 1
        elif states[i] != states[i - 1]:
            cur = 2
        else:
            cur = 1
        best = max(best, cur)
    if not any(s[6] for s in segments):
        best = 0
    return (best if best >= 2 else 0), len(set(states))


# ─── breakpoint clustering (criterion 7) ────────────────────────────────────
class EffectiveCoords:
    """Chromosome coordinates with masked intervals removed (segdups / hotspots)."""

    def __init__(self, chrom_len, masked):
        merged = []
        for s, e in sorted(masked):
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        self.iv = merged
        self.starts = [s for s, _ in merged]
        self.cum = [0]
        for s, e in merged:
            self.cum.append(self.cum[-1] + (e - s))
        self.length = chrom_len - self.cum[-1]

    def masked(self, pos):
        i = bisect.bisect_right(self.starts, pos) - 1
        return i >= 0 and self.iv[i][0] <= pos < self.iv[i][1]

    def eff(self, pos):
        i = bisect.bisect_right(self.starts, pos)
        before = self.cum[i]                     # masked bp in intervals starting at or before pos
        if i > 0 and pos < self.iv[i - 1][1]:
            before -= self.iv[i - 1][1] - pos
        return pos - before


def clustering_test(locus_pos, chrom_pos, eff, n_sim, seed):
    """Scan statistic: is the span of the locus's k breakends shorter than the
    shortest span of k consecutive breakends among n uniformly placed ones (n =
    the chromosome's non-public, unmasked breakends), on the masked chromosome?
    Also the KS distance of the locus's breakend spacings to the exponential
    expected under uniform breakage (rate n / effective length). Returns a dict;
    p is a simulated, descriptive value."""
    import numpy as np
    k, n = len(locus_pos), max(len(chrom_pos), len(locus_pos))
    if k < 3 or eff.length <= 0:
        return None
    x = sorted(eff.eff(p) for p in locus_pos)
    span = x[-1] - x[0]
    rate = n / eff.length
    gaps = sorted(b - a for a, b in zip(x, x[1:]))
    m = len(gaps)
    D = max(max(abs((i + 1) / m - (1 - math.exp(-rate * g))), abs(i / m - (1 - math.exp(-rate * g))))
            for i, g in enumerate(gaps))
    rng = np.random.default_rng(seed)
    sims = np.sort(rng.uniform(0, eff.length, size=(n_sim, n)), axis=1)
    min_span = (sims[:, k - 1:] - sims[:, :n - k + 1]).min(axis=1)
    p = (1 + int((min_span <= span).sum())) / (n_sim + 1)
    return {"k": k, "n_chrom": n, "span_eff": span, "expected_spacing": eff.length / n,
            "observed_mean_spacing": span / (k - 1), "ks_D": D, "p_scan": p}


def bh_qvalues(pvals):
    """Benjamini-Hochberg q-values; None entries stay None."""
    idx = [i for i, p in enumerate(pvals) if p is not None]
    q = [None] * len(pvals)
    if not idx:
        return q
    order = sorted(idx, key=lambda i: pvals[i])
    m, prev = len(order), 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        prev = min(prev, pvals[i] * m / rank)
        q[i] = prev
    return q


def seed_of(text, base=0):
    return (int(hashlib.sha256(text.encode()).hexdigest()[:12], 16) + base) % (2 ** 32)


# ─── HPO ────────────────────────────────────────────────────────────────────
def read_obo_parents(path):
    parents, cur = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line == "[Term]":
                cur = None
            elif line.startswith("id: HP:"):
                cur = line[4:]
                parents.setdefault(cur, set())
            elif cur and line.startswith("is_a: "):
                parents[cur].add(line[6:].split(" ", 1)[0])
            elif cur and line.startswith("alt_id: "):
                parents.setdefault(line[8:], set()).add(cur)
    return parents


def ancestors(term, parents, _cache={}):
    key = (id(parents), term)
    if key in _cache:
        return _cache[key]
    out, stack = set(), [term]
    while stack:
        t = stack.pop()
        for p in parents.get(t, ()):
            if p not in out:
                out.add(p)
                stack.append(p)
    _cache[key] = out
    return out


def hpo_match(patient_terms, gene_terms, parents):
    """Patient terms matched by the genes' annotations: a patient term T matches
    when a gene annotation A equals T, is a descendant of T, or is an ancestor of
    T other than the root ('All' / 'Phenotypic abnormality'). Returns the matched
    patient terms."""
    roots = {"HP:0000001", "HP:0000118"}
    matched = set()
    for T in patient_terms:
        anc_T = ancestors(T, parents) - roots
        for A in gene_terms:
            if A == T or T in ancestors(A, parents) or A in anc_T:
                matched.add(T)
                break
    return matched


# ─── gene points (ClinGen-style section 2) ───────────────────────────────────
def gene_points(mode, hi, ts, loeuf, w):
    """Points of one affected gene. mode: broken | lost | gained."""
    def score(table, s):
        return float(table.get(str(s), 0.0)) if s not in (None, "", "NA") else 0.0
    if mode == "gained":
        return score(w["ts_gained"], ts), f"TS {ts or 'NA'} (whole gene gained)"
    table = w["hi_lost"] if mode == "lost" else w["hi_broken"]
    pts = score(table, hi)
    why = f"HI {hi or 'NA'} ({'whole gene lost' if mode == 'lost' else 'gene broken by a breakend'})"
    # LOEUF only where ClinGen has no HI evidence (not curated, 0 or 30); never over a curated 40
    if pts == 0.0 and loeuf is not None and str(hi) not in ("1", "2", "3", "40"):
        for cut, pt in sorted(((float(k), float(v)) for k, v in w["loeuf"].items())):
            if loeuf < cut:
                return pt, why + f"; LOEUF {loeuf:.3g} < {cut:g}"
    return pts, why


# ─── feature functions (each returns a Feature) ──────────────────────────────
def f_panel_privacy(status, carriers, n_panel, n_callable, w):
    """Panel status weighted by callability: a locus absent from the panel only
    counts as private evidence where the panel could have seen it."""
    if status == "public":
        return Feature("panel_privacy", "privacy", f"public ({carriers}/{n_panel} carriers)",
                       float(w["public"]), "rearrangement carried by panel individuals")
    if status == "shared_junctions":
        return Feature("panel_privacy", "privacy", f"shared_junctions (0 carriers; callable {n_callable}/{n_panel})",
                       float(w["shared_junctions"]), "some junctions in the panel, not the rearrangement")
    try:
        frac = int(n_callable) / int(n_panel) if int(n_panel) else 0.0
        why = f"private, callable in {n_callable}/{n_panel} panel individuals"
    except (TypeError, ValueError):
        frac, why = 0.5, "private, panel callability unknown (x 0.5)"
    return Feature("panel_privacy", "privacy", f"private (0/{n_panel} carriers; callable {n_callable}/{n_panel})",
                   r2(float(w["private_callable"]) * frac), why)


def f_catalogue_frequency(cat_rows, panel_points, w):
    """cat_rows: [(name, matched 'k/n', max_af or None, contains bool)]. The most
    severe containing catalogue sets the penalty; privacy = min(panel, catalogue)
    when either is negative, so points here are that minimum minus the panel's."""
    contained = [(n, af) for n, k, af, c in cat_rows if c]
    shown = "; ".join(f"{n} {k}" + (f" AF {af:.3g}" if af is not None else "") for n, k, af, c in cat_rows) or "none configured"
    if not contained:
        return Feature("catalogue_frequency", "privacy", shown, 0.0, "no catalogue contains the locus")
    worst, why = 0.0, ""
    for n, af in contained:
        if af is not None and af >= float(w["common_af"]):
            pen, r = float(w["common"]), f"{n} AF {af:.3g} >= {w['common_af']}"
        elif af is not None and af >= float(w["rare_af"]):
            pen, r = float(w["rare"]), f"{n} AF {af:.3g} >= {w['rare_af']}"
        else:
            pen, r = float(w["seen"]), f"{n} contains it (AF {fnum(af)})"
        if pen < worst or not why:
            worst, why = min(worst, pen), r
    combined = min(worst, panel_points) if panel_points < 0 else worst
    return Feature("catalogue_frequency", "privacy", shown, r2(combined - panel_points),
                   why + ("; privacy = min(panel, catalogue)" if panel_points < 0 else "; replaces the private bonus"))


def f_gene_dosage(affected, dosage, constraint, w):
    """affected: {gene: mode} (broken | lost | gained). Max over genes."""
    best, best_why, shown = 0.0, "no dosage-sensitive gene affected", []
    for g, mode in sorted(affected.items()):
        d = dosage.get(g, {})
        loeuf = constraint.get(g)
        pts, why = gene_points(mode, d.get("hi"), d.get("ts"), loeuf, w)
        shown.append(f"{g}:{mode}(HI {d.get('hi', 'NA')},TS {d.get('ts', 'NA')},LOEUF {fnum(loeuf)})")
        if pts > best:
            best, best_why = pts, f"{g}: {why}"
    return Feature("gene_dosage", "genes", "; ".join(shown) or "none", r2(best), best_why)


def f_exons(n_exons, w):
    return Feature("exons", "genes", n_exons, float(w) if n_exons else 0.0,
                   "exons in the rearranged interval or hit by a breakend" if n_exons else "no exon disrupted")


def f_gene_count(n_genes, w):
    pts = 0.0
    for k, v in sorted(((int(k), float(v)) for k, v in w.items())):
        if n_genes >= k:
            pts = v
    return Feature("gene_count", "genes", n_genes, pts, f"{n_genes} genes affected (ClinGen section 3)")


def f_no_genes(n_genes, n_exons, w):
    if n_genes or n_exons:
        return Feature("no_genes", "genes", "no", 0.0, "genes affected")
    return Feature("no_genes", "genes", "yes", float(w), "no gene broken, lost or gained (ClinGen 1B)")


def f_class(cls, w):
    return Feature("complexity_class", "event", cls, float(w.get(cls, w.get("default", 0.0))), "class points")


def f_size(span, w):
    pts = 0.0
    for k, v in sorted(((int(k), float(v)) for k, v in w.items())):
        if span >= k:
            pts = v
    return Feature("size", "event", span, pts, "locus span (bp)")


def f_caller_agreement(supp_vecs, w):
    n = len(supp_vecs)
    multi = sum(1 for s in supp_vecs if s.count("1") >= 2)
    frac = multi / n if n else 0.0
    if n and frac >= float(w["multi_fraction"]):
        pts, why = float(w["multi"]), f">= {w['multi_fraction']:.0%} of junctions with >= 2 callers"
    elif n and multi == 0:
        pts, why = float(w["none"]), "every junction from a single caller"
    else:
        pts, why = 0.0, "mixed caller support"
    return Feature("caller_agreement", "evidence", f"{multi}/{n} junctions with >= 2 callers", pts, why)


def f_read_support(min_reads, w):
    if min_reads is None:
        return Feature("read_support", "evidence", "NA", 0.0, "read counts not available")
    if min_reads < int(w["min_reads"]):
        return Feature("read_support", "evidence", min_reads, float(w["weak"]),
                       f"a junction has < {w['min_reads']} supporting reads")
    return Feature("read_support", "evidence", min_reads, 0.0, f"every junction >= {w['min_reads']} reads")


def f_artefact(frac_segdup, frac_lowc, w):
    pts = float(w["segdup"]) * frac_segdup + float(w["low_complexity"]) * frac_lowc
    return Feature("artefact_context", "evidence", f"segdup {frac_segdup:.2f}, low-complexity {frac_lowc:.2f}",
                   r2(pts), "fraction of breakends in segmental duplications / simple repeats")


def f_interleaving(max_il, w):
    pts = 0.0
    for k, v in sorted(((int(k), float(v)) for k, v in w.items())):
        if max_il >= k:
            pts = v
    return Feature("ct_interleaving", "chromothripsis", max_il, pts,
                   "largest interleaved junction component (criterion 1; floor 3-4, cancer 6)")


def f_breakpoint_tier(tier, n_bk, w):
    return Feature("ct_breakpoint_tier", "chromothripsis", f"{tier} ({n_bk} breakends)", float(w.get(tier, 0.0)),
                   "criterion 2: intermediate < 5, extreme > 5 breakpoints")


def f_cn_oscillation(segs, w, sex_note):
    if segs is None:
        return Feature("ct_cn_oscillation", "chromothripsis", "NA (no depth baseline)", 0.0,
                       "criterion 5: copy number not assessable")
    changed = sum(1 for s in segs if s[6])
    run, n_states = oscillation(segs)
    pts = float(w["change"]) if changed else 0.0
    if run >= int(w["osc_min_segments"]) and n_states >= 2:
        pts += float(w["oscillation"])
    states = ",".join(str(s[4]) for s in segs)
    return Feature("ct_cn_oscillation", "chromothripsis",
                   f"CN states {states}; {changed} changed; oscillating run {run}; {sex_note}", r2(pts),
                   "criterion 5 (supporting, never a gate)")


def f_clustering(res, q, w):
    if res is None:
        return Feature("ct_clustering", "chromothripsis", "NA (< 3 unmasked breakends)", 0.0, "criterion 7")
    val = (f"k={res['k']} of n={res['n_chrom']} chromosome breakends; mean spacing "
           f"{res['observed_mean_spacing']:.0f} vs {res['expected_spacing']:.0f} bp expected; KS D "
           f"{res['ks_D']:.2f}; p_sim {res['p_scan']:.3g}; q {fnum(q)}")
    pts = float(w["points"]) if q is not None and q <= float(w["q_max"]) else 0.0
    return Feature("ct_clustering", "chromothripsis", val, pts,
                   f"criterion 7, chromosome-intrinsic simulated null (descriptive; q <= {w['q_max']})")


def f_bfb(ct_subtotal, counts, p_balance, n_fold, n_fold_cn, n_inv, subtelomeric, p_fold, w):
    frac = n_fold / n_inv if n_inv else 0.0
    bfb = n_fold >= 1 and frac >= float(p_fold["min_fraction"]) and n_fold_cn >= 1
    val = (f"joins {counts}; balance p {fnum(p_balance)}; fold-backs {n_fold}/{n_inv} inversion-type "
           f"({n_fold_cn} with CN gain); subtelomeric {'yes' if subtelomeric else 'no'}")
    if not bfb or ct_subtotal <= 0:
        return Feature("ct_bfb_modifier", "chromothripsis", val, 0.0,
                       "criterion 8: no BFB-like fold-back skew" if not bfb else "criterion 8: nothing to modify")
    f = float(w["subtelomeric_factor"] if subtelomeric else w["factor"])
    return Feature("ct_bfb_modifier", "chromothripsis", val, r2(-f * ct_subtotal),
                   f"criterion 8: BFB-like (fold-back skew with CN gain{', subtelomeric' if subtelomeric else ''}); "
                   f"chromothripsis points x {1 - f:g}")


def f_hpo(patient_terms, matched, w):
    if not patient_terms:
        return Feature("hpo_match", "phenotype", "no patient HPO terms", 0.0, "optional")
    n, k = len(patient_terms), len(matched)
    pts = float(w["half"]) if k and k / n >= 0.5 else (float(w["any"]) if k else 0.0)
    return Feature("hpo_match", "phenotype", f"{k}/{n} patient terms matched ({','.join(sorted(matched)) or '-'})",
                   pts, "HPO annotations of the affected genes")


def f_omim(genes_with_phenotype, enabled, w):
    if not enabled:
        return Feature("omim", "phenotype", "off (no licence key)", 0.0, "optional")
    return Feature("omim", "phenotype", ",".join(sorted(genes_with_phenotype)) or "none",
                   float(w["phenotype_gene"]) if genes_with_phenotype else 0.0, "affected gene with an OMIM phenotype")


def f_inheritance(parent_results, parents_affected, w):
    """parent_results: {role: carries bool}; both parents needed for de novo."""
    if not parent_results:
        return Feature("inheritance", "inheritance", "parents not sequenced", 0.0, "optional")
    carriers = [r for r, c in parent_results.items() if c]
    if carriers:
        pts = float(w["inherited_affected"] if parents_affected else w["inherited_unaffected"])
        return Feature("inheritance", "inheritance", f"inherited ({','.join(carriers)})", pts,
                       "carried by an " + ("affected" if parents_affected else "unaffected") + " parent")
    if len(parent_results) >= 2:
        return Feature("inheritance", "inheritance", "de novo (absent in both parents)", float(w["de_novo"]),
                       "ClinGen 5A-like")
    return Feature("inheritance", "inheritance", f"absent in {','.join(parent_results)}; other parent not sequenced",
                   0.0, "de novo status needs both parents")
