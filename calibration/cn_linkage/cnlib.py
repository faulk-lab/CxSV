"""Copy-number-anchored linkage, Phase 1 prototype (roadmap item 15; docs/CN_LINKAGE_PHASE1.md).

Outside the pipeline: nothing here is read by workflow/ and nothing changes the
calling fingerprint. Pure numpy/pandas; runs in the pipeline image or any env with both.

Pieces
  read_bins          mosdepth --by 1000 regions.bed.gz -> {contig: 1 kb mean depth}
  to_coarse          1 kb -> 10 kb bins (mean of the 1 kb bins, length-weighted)
  genome_median      autosomal median of the 10 kb bins (> 0)
  ploidy_track       expected copies per 10 kb bin from sex (+ PAR on chrX in males)
  PanelOfNormals     per platform stratum: matrix (bins x genomes) of normalised depth
                     (depth / genome median); per-bin median and MAD over the healthy
                     genomes of the stratum, leaving out the genome's own family, same
                     sex only on chrX / chrY
  segment            ratio = depth / (genome median x PoN median); ploidy-aware HMM over
                     copy-number states (Viterbi), >= min_bins per segment, then each
                     boundary refined on the 1 kb bins
  The candidate / junction / structure model is in cnlink.py.
"""
import gzip
import math
from collections import defaultdict

import numpy as np
try:
    import pandas as pd
except ImportError:          # simulation env (badread) has numpy only
    pd = None

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
SEX_CHROMS = ["chrX", "chrY"]
CONTIGS = AUTOSOMES + SEX_CHROMS
COARSE = 10          # 1 kb bins per coarse bin (10 kb)

DEFAULTS = dict(
    min_pon_median=0.25,     # PoN median (normalised, per expected copy) below this: bin masked
    max_pon_relmad=0.30,     # PoN robust relative spread above this: bin masked
    min_pon_genomes=4,       # fewer healthy genomes with data in the bin: masked
    p_change=1e-4,           # HMM per-bin probability of a copy-number change
    outlier=0.01,            # HMM emission mixture weight of a flat outlier density
    max_cn=6,                # highest copy-number state (ploidy 2); ploidy 1: max_cn - 2
    min_bins=3,              # shortest segment (coarse bins)
    min_cn_delta=0.6,        # |mean copy number - ploidy| a segment needs (as config/scoring.yaml)
    refine_kb=15,            # boundary refinement window on the 1 kb bins (+/- kb)
    nested_min_seg_kb=60,    # gain segments at least this long are scanned for a nested level
    nested_min_kb=30,        # shortest nested higher level (1 kb bins)
    nested_min_z=6.0,        # nested level vs the segment's level, in 10 kb-bin standard errors
    nested_min_copy=0.4,     # nested level at least this many copies above the rest of the segment
                             # (observed nested levels are compressed: PC L01 TRP +0.6 copy, not +1)
)


# ─── input ──────────────────────────────────────────────────────────────────
def read_fai(path):
    return {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(path)}


def read_bins(path, contigs=CONTIGS):
    """{contig: float32 array of 1 kb mean depth} from mosdepth regions.bed.gz."""
    if pd is None:
        import gzip
        vals = defaultdict(list)
        with gzip.open(path, "rt") as fh:
            for l in fh:
                c, _, _, d = l.split("\t")
                if c in contigs:
                    vals[c].append(float(d))
        return {c: np.asarray(v, np.float32) for c, v in vals.items()}
    df = pd.read_csv(path, sep="\t", header=None, names=["c", "s", "e", "d"],
                     dtype={"c": str, "s": np.int64, "e": np.int64, "d": np.float32})
    out = {}
    for c, g in df.groupby("c", sort=False):
        if c in contigs:
            out[c] = g["d"].to_numpy(np.float32)
    return out


def to_coarse(fine, k=COARSE):
    """1 kb -> k kb bins (mean; a trailing partial bin is the mean of what it has)."""
    out = {}
    for c, a in fine.items():
        n = len(a)
        m = int(math.ceil(n / k))
        pad = np.full(m * k, np.nan, np.float32)
        pad[:n] = a
        out[c] = np.nanmean(pad.reshape(m, k), axis=1).astype(np.float32)
    return out


def genome_median(coarse):
    v = np.concatenate([coarse[c] for c in AUTOSOMES if c in coarse])
    v = v[v > 0]
    return float(np.median(v)) if len(v) else float("nan")


def read_par(path):
    """PAR intervals {contig: [(start, end)]} (BED; resources/scoring/par.bed)."""
    out = defaultdict(list)
    if path:
        for l in open(path):
            f = l.split()
            if len(f) >= 3 and not l.startswith("#"):
                out[f[0]].append((int(f[1]), int(f[2])))
    return out


def ploidy_track(contig, nbins, sex, par, binsize=COARSE * 1000):
    """Expected copies per coarse bin: autosomes 2; chrX 2 (female) / 1 (male, PAR 2);
    chrY 1 outside PAR in males, 0 in females (PAR on chrY: 0, counted on chrX)."""
    if contig in AUTOSOMES:
        return np.full(nbins, 2, np.int8)
    p = np.zeros(nbins, np.int8)
    if contig == "chrX":
        p[:] = 2 if sex == "female" else 1
        if sex == "male":
            for s, e in par.get("chrX", []):
                p[s // binsize:int(math.ceil(e / binsize))] = 2
    elif contig == "chrY":
        if sex == "male":
            p[:] = 1
            for s, e in par.get("chrY", []):
                p[s // binsize:int(math.ceil(e / binsize))] = 0
    return p


# ─── panel of normals ───────────────────────────────────────────────────────
class PanelOfNormals:
    """Normalised coarse depth of the healthy genomes of one platform stratum.

    Values are depth / genome median, divided by the bin's expected copies / 2, so a
    normal bin is ~1.0 whatever the sex (chrX in a male: 0.5 / (1/2) = 1)."""

    def __init__(self, genomes, coarse_by_genome, sexes, families, par):
        self.genomes = list(genomes)
        self.sex = np.array([sexes[g] for g in self.genomes])
        self.family = np.array([families[g] for g in self.genomes])
        self.mat = {}
        for c in CONTIGS:
            if any(coarse_by_genome[g][0].get(c) is None for g in self.genomes):
                continue                      # contig absent from a genome: no reference
            cols = []
            for g in self.genomes:
                a = coarse_by_genome[g][0].get(c)
                med = coarse_by_genome[g][1]
                pl = ploidy_track(c, len(a), sexes[g], par).astype(np.float32)
                with np.errstate(divide="ignore", invalid="ignore"):
                    cols.append(np.where(pl > 0, a / med / (pl / 2.0), np.nan))
            self.mat[c] = np.vstack(cols).T.astype(np.float32)   # bins x genomes

    @classmethod
    def from_npz(cls, path):
        """The PoN written by segment_all.py (pon/<stratum>.npz)."""
        z = np.load(path, allow_pickle=False)
        self = cls.__new__(cls)
        self.genomes = list(z["genomes"])
        self.sex = z["sex"]
        self.family = z["family"]
        self.mat = {c: z[c] for c in CONTIGS if c in z.files}
        return self

    def reference(self, contig, sex, exclude_family=None):
        """(median, relative MAD, n) per bin over the stratum's genomes, without the
        given family. chrX pools both sexes (values are per expected copy, PAR-aware),
        which doubles the reference at the X-linked loci (a same-sex chrX reference of
        7 R10 females was too noisy to resolve the PC L01 triplication); chrY: males."""
        keep = np.ones(len(self.genomes), bool)
        if exclude_family is not None:
            keep &= self.family != exclude_family
        if contig == "chrY":
            keep &= self.sex == "male"
        M = self.mat[contig][:, keep]
        with np.errstate(all="ignore"):
            n = np.sum(np.isfinite(M) & (M > 0), axis=1)
            med = np.nanmedian(M, axis=1)
            mad = 1.4826 * np.nanmedian(np.abs(M - med[:, None]), axis=1)
            rel = mad / med
        return med, rel, n


# ─── segmentation ───────────────────────────────────────────────────────────
def _viterbi(r, sd, ploidy, kmax, p_change, outlier):
    """Most likely copy-number path; r = ratio (NaN = masked), sd per bin."""
    K = kmax + 1
    cn = np.arange(K, dtype=np.float64)
    mu = cn / ploidy
    scale = np.sqrt(np.maximum(cn, 0.5) / ploidy)
    lstay, lsw = math.log(1 - p_change), math.log(p_change / (K - 1))
    n = len(r)
    back = np.zeros((n, K), np.int8)
    prior = np.full(K, math.log(1e-3 / (K - 1)))
    prior[ploidy] = math.log(1 - 1e-3)
    v = prior.copy()
    lflat = math.log(outlier / 4.0)            # flat density over ratio 0..4
    for t in range(n):
        if t:
            best_other = v.max()
            arg_other = int(v.argmax())
            stay = v + lstay
            sw = best_other + lsw
            take_sw = sw > stay
            back[t] = np.where(take_sw, arg_other, np.arange(K))
            v = np.where(take_sw, sw, stay)
        x = r[t]
        if np.isfinite(x):
            s = sd[t] * scale
            ln = -0.5 * ((x - mu) / s) ** 2 - np.log(s) - 0.9189385
            v = v + np.logaddexp(math.log(1 - outlier) + ln, lflat)
    path = np.empty(n, np.int8)
    path[-1] = int(v.argmax())
    for t in range(n - 1, 0, -1):
        path[t - 1] = back[t, path[t]]
    return path


def runs(states):
    """[(start_bin, end_bin_exclusive, state)] of constant runs."""
    if len(states) == 0:
        return []
    cut = np.flatnonzero(np.diff(states)) + 1
    bounds = np.concatenate([[0], cut, [len(states)]])
    return [(int(a), int(b), int(states[a])) for a, b in zip(bounds[:-1], bounds[1:])]


def segment_contig(contig, coarse_depth, fine_depth, gmed, ploidy, pon_med, pon_rel, pon_n, sigma_g, P):
    """Abnormal segments of one contig: dicts with chrom, start, end (bp, refined),
    state (copies), ploidy, mean_ratio, n_bins, z."""
    n = len(coarse_depth)
    # ratio: depth relative to the healthy panel's depth for this bin (normal = 1;
    # one extra copy on ploidy 2 = 1.5, on ploidy 1 = 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = coarse_depth / gmed / pon_med / (np.maximum(ploidy, 1) / 2.0)
    mask = ((pon_med < P["min_pon_median"]) | ~np.isfinite(pon_med) | (pon_rel > P["max_pon_relmad"])
            | (pon_n < P["min_pon_genomes"]) | (ploidy == 0) | ~np.isfinite(ratio))
    r = np.where(mask, np.nan, ratio)
    sd = np.sqrt(sigma_g ** 2 + np.nan_to_num(pon_rel, nan=1.0) ** 2)
    out = []
    # piecewise by ploidy (PAR in male chrX): run the HMM on each constant-ploidy stretch
    for a, b, pl in runs(ploidy):
        if pl == 0 or b - a < P["min_bins"]:
            continue
        kmax = P["max_cn"] if pl == 2 else P["max_cn"] - 2
        path = _viterbi(r[a:b], sd[a:b], int(pl), kmax, P["p_change"], P["outlier"]).astype(np.int16)
        path = _clean(path, r[a:b], int(pl), P)
        for s, e, st in runs(path):
            if st == pl:
                continue
            seg = r[a + s:a + e]
            ok = np.isfinite(seg)
            if ok.sum() < P["min_bins"]:
                continue
            m = float(np.nanmean(seg))
            if abs(m * pl - pl) < P["min_cn_delta"]:
                continue
            se = float(np.nanmean(sd[a + s:a + e][ok])) / math.sqrt(ok.sum())
            out.append(dict(chrom=contig, b0=a + s, b1=a + e, state=int(st), ploidy=int(pl),
                            mean_ratio=m, cn=m * pl, n_bins=int(ok.sum()), z=(m - 1) / se if se > 0 else 0.0))
    # merge adjacent segments with the same state
    merged = []
    for sg in out:
        if merged and merged[-1]["b1"] == sg["b0"] and merged[-1]["state"] == sg["state"]:
            merged[-1]["b1"] = sg["b1"]
            merged[-1]["n_bins"] += sg["n_bins"]
        else:
            merged.append(sg)
    for sg in merged:
        sg["start"], sg["end"] = sg["b0"] * COARSE * 1000, sg["b1"] * COARSE * 1000
    _refine(merged, fine_depth, gmed, pon_med, ploidy, P)
    merged = _nested(merged, fine_depth, gmed, pon_med, ploidy, mask, P)
    if any(sg.get("nested") for sg in merged):
        _refine(merged, fine_depth, gmed, pon_med, ploidy, P)      # nested boundaries to 1 kb
    for sg in merged:
        sg["start_mask_before"], sg["start_mask_after"] = _masked_ext(mask, sg["start"])
        sg["end_mask_before"], sg["end_mask_after"] = _masked_ext(mask, sg["end"])
    return merged, r



def _nested(segs, fine, gmed, pon_med, ploidy, mask, P):
    """Split a gain segment around a nested higher level (the TRP of a DUP-TRP/INV-DUP
    is often too short for the genome-wide HMM's change penalty).

    On the informative 1 kb bins of the segment, the interval [i, j) with the largest
    two-sample statistic (mean inside vs mean of the rest; circular binary
    segmentation) is found. It becomes a nested level (state + 1) if it is higher, is
    >= nested_min_kb bins with >= nested_min_kb bins outside, rises by >=
    nested_min_copy copies, and the difference is >= nested_min_z standard errors, the
    error taken from the spread of the segment's 10 kb blocks (long reads correlate
    neighbouring 1 kb bins). Observed nested levels are compressed: the PC L01 TRP
    reads +0.5 copy above its DUP, not +1."""
    out = []
    m_min = P["nested_min_kb"]
    for sg in segs:
        pl, st = sg["ploidy"], sg["state"]
        if st <= pl or sg["end"] - sg["start"] < P["nested_min_seg_kb"] * 1000:
            out.append(sg)
            continue
        a, b = sg["start"] // 1000, min(len(fine), sg["end"] // 1000)
        cb = np.minimum(np.arange(a, b) // COARSE, len(pon_med) - 1)
        with np.errstate(divide="ignore", invalid="ignore"):
            rr = fine[a:b] / gmed / pon_med[cb] / (max(pl, 1) / 2.0)
        ok = ~mask[cb] & np.isfinite(rr)
        pos = np.nonzero(ok)[0]
        y = rr[ok].astype(np.float64)
        n = len(y)
        if n < 2 * m_min:
            out.append(sg)
            continue
        cs = np.concatenate([[0.0], np.cumsum(y)])
        tot = cs[-1]
        best, bi, bj = 0.0, 0, 0
        for i in range(0, n - m_min + 1):
            j = np.arange(i + m_min, n + 1)
            k = j - i
            rest = n - k
            good = rest >= m_min
            if not good.any():
                break
            j, k, rest = j[good], k[good], rest[good]
            s_in = cs[j] - cs[i]
            d = s_in / k - (tot - s_in) / rest
            t = d * np.sqrt(k * rest / n)          # proportional to the two-sample z
            q = int(np.argmax(t))
            if t[q] > best:
                best, bi, bj = float(t[q]), i, int(j[q])
        if best <= 0:
            out.append(sg)
            continue
        inside = y[bi:bj]
        rest_y = np.concatenate([y[:bi], y[bj:]])
        diff = float(inside.mean() - rest_y.mean())
        k = len(rest_y) // COARSE
        blocks = rest_y[:k * COARSE].reshape(k, COARSE).mean(axis=1) if k >= 3 else np.array([])
        sd10 = float(1.4826 * np.median(np.abs(blocks - np.median(blocks)))) if len(blocks) >= 3 else float("inf")
        se = max(sd10, 1e-3) * math.sqrt(COARSE / len(inside) + COARSE / len(rest_y))
        z = diff / se
        if diff * pl < P["nested_min_copy"] or z < P["nested_min_z"]:
            out.append(sg)
            continue
        s1 = (a + int(pos[bi])) * 1000
        e1 = (a + int(pos[bj - 1]) + 1) * 1000
        for s0_, e0_, k_ in ((sg["start"], s1, st), (s1, e1, st + 1), (e1, sg["end"], st)):
            if e0_ - s0_ < 1000:
                continue
            i0, i1 = s0_ // 1000 - a, e0_ // 1000 - a
            v = rr[i0:i1][ok[i0:i1]]
            mr = float(v.mean()) if len(v) else float("nan")
            piece = dict(sg, start=s0_, end=e0_, state=k_, mean_ratio=mr, cn=mr * pl,
                         n_bins=int(len(v) // COARSE), nested=(k_ == st + 1))
            piece["b0"], piece["b1"] = s0_ // (COARSE * 1000), -(-e0_ // (COARSE * 1000))
            if k_ == st + 1:
                piece["z"] = round(z, 2)
            out.append(piece)
    return out



def _masked_ext(mask, pos):
    """bp of masked 10 kb bins immediately before / after a boundary: the HMM carries a
    state through masked bins, so the true change point can lie anywhere in them."""
    k = pos // (COARSE * 1000)
    before = after = 0
    i = k - 1
    while i >= 0 and mask[i]:
        before += 1
        i -= 1
    i = k
    while i < len(mask) and mask[i]:
        after += 1
        i += 1
    return before * COARSE * 1000, after * COARSE * 1000

def _clean(path, r, pl, P):
    """Abnormal runs shorter than min_bins informative bins return to ploidy."""
    path = path.copy()
    for s, e, st in runs(path):
        if st != pl and np.isfinite(r[s:e]).sum() < P["min_bins"]:
            path[s:e] = pl
    return path


def _level(sg):
    """Ratio level of a segment for boundary refinement: its observed mean (nested
    levels and some duplications read lower than state / ploidy), else nominal;
    outside any segment 1.0."""
    if sg is None:
        return 1.0
    m = sg.get("mean_ratio")
    return float(m) if m is not None and np.isfinite(m) else sg["state"] / sg["ploidy"]


def _refine(segs, fine, gmed, pon_med, ploidy, P):
    """Move each segment boundary to the 1 kb position that best splits the two levels."""
    w = P["refine_kb"]
    n_f = len(fine)
    # boundaries: (position, segment ending there or None, segment starting there or None);
    # adjacent abnormal segments share one boundary, refined once with both levels
    bmap = {}
    for sg in segs:
        bmap.setdefault(sg["start"], [None, None])[1] = sg
        bmap.setdefault(sg["end"], [None, None])[0] = sg
    for pos, (left_seg, right_seg) in sorted(bmap.items()):
        left, right = _level(left_seg), _level(right_seg)
        pos_kb = pos // 1000
        # the coarse path can be a bin or two off; re-centre while the best split sits on
        # the window edge (at most 3 moves)
        for _ in range(4):
            lo, hi = max(0, pos_kb - w), min(n_f, pos_kb + w)
            if hi - lo < 4:
                break
            idx = np.arange(lo, hi)
            cb = np.minimum(idx // COARSE, len(pon_med) - 1)
            pl = np.maximum(ploidy[cb], 1)
            with np.errstate(divide="ignore", invalid="ignore"):
                rr = fine[lo:hi] / gmed / pon_med[cb] / (pl / 2.0)
            ok = np.isfinite(rr)
            if ok.sum() < 4:
                break
            x = np.where(ok, rr, 0.0)
            wgt = ok.astype(float)
            cl = np.cumsum(wgt * (x - left) ** 2)
            cr = np.cumsum((wgt * (x - right) ** 2)[::-1])[::-1]
            cost = np.concatenate([[0], cl]) + np.concatenate([cr, [0]])
            k = int(np.argmin(cost))
            new_kb = lo + k
            edge = (k <= 1 and lo > 0) or (k >= len(cost) - 2 and hi < n_f)
            if new_kb == pos_kb or not edge:
                pos_kb = new_kb
                break
            pos_kb = new_kb
        newpos = int(pos_kb * 1000)
        if left_seg:
            left_seg["end"] = newpos
        if right_seg:
            right_seg["start"] = newpos
    for sg in segs:
        if sg["end"] <= sg["start"]:
            sg["end"] = sg["start"] + 1000


def genome_noise(coarse, gmed, pon, sex, family, par, P):
    """Robust spread of the ratio over unmasked autosomal bins."""
    vals = []
    for c in AUTOSOMES:
        med, rel, n = pon.reference(c, sex, family)
        with np.errstate(divide="ignore", invalid="ignore"):
            r = coarse[c] / gmed / med
        ok = np.isfinite(r) & (med >= P["min_pon_median"]) & (rel <= P["max_pon_relmad"]) & (n >= P["min_pon_genomes"])
        vals.append(r[ok])
    v = np.concatenate(vals)
    return float(1.4826 * np.median(np.abs(v - np.median(v))))


def segment_genome(fine, sex, par, pon, family=None, P=DEFAULTS, contigs=CONTIGS):
    """All abnormal segments of a genome against a PoN (family left out)."""
    coarse = to_coarse(fine)
    gmed = genome_median(coarse)
    sigma = genome_noise(coarse, gmed, pon, sex, family, par, P)
    segs = []
    for c in contigs:
        if c not in coarse or (c == "chrY" and sex != "male"):
            continue
        med, rel, n = pon.reference(c, sex, family)
        pl = ploidy_track(c, len(coarse[c]), sex, par)
        s, _ = segment_contig(c, coarse[c], fine[c], gmed, pl, med, rel, n, sigma, P)
        segs.extend(s)
    return segs, dict(genome_median=gmed, sigma=sigma)
