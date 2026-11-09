"""Sex-stratified chrX / chrY copy-number reference (roadmap item 15, panel-v1.2-r2).

panel-v1.2 pooled both sexes in the chrX PoN. Where male chrX depth per copy differs from
female depth, males got spurious gains (state 2 at ploidy 1: ONT-1019 males 5.5 chrX gains
per genome vs 0.3 in females), and the segment-level panel frequency then made the positive
control's PLP1 DUP-TRP (L01, a female) public. These tests pin the fix: sex call from chrX and
chrY, one chrX reference per sex, ambiguous calls out of the sex-chromosome reference, and
sex-matched public counting on chrX / chrY. Synthetic depth only; numpy, no pandas.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "workflow" / "lib"))
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))
import cxsv_cn as cn       # noqa: E402
import cxsv_cnlink as cnl  # noqa: E402

NB = 4000
PAR = {"chrX": [(0, 200_000)], "chrY": [(0, 100_000)]}
P = dict(cn.DEFAULTS)


def _fine(rng, cov, sex, events=(), male_extra=None, ychr=True, y_lost=False):
    """chr1, chrX (and chrY for males when ychr) Poisson depth; male_extra = (s, e) 1 kb bins
    where males read twice their copy number on chrX (a male-specific depth artefact)."""
    out = {}
    for c in ("chr1", "chrX") + (("chrY",) if ychr else ()):
        if c == "chrY":
            lam = np.full(NB, cov / 2.0 if sex == "male" and not y_lost else 0.0)
        else:
            pl = 2 if (c != "chrX" or sex == "female") else 1
            lam = np.full(NB, cov * pl / 2.0)
            if c == "chrX" and sex == "male":
                lam[:200] = cov
                if male_extra:
                    lam[male_extra[0]:male_extra[1]] = cov
        for ec, s, e, k in events:
            if ec == c:
                lam[s:e] = cov * k / 2.0
        out[c] = rng.poisson(lam).astype(np.float32)
    return out


def _member(fine):
    coarse = cn.to_coarse(fine)
    gmed = cn.genome_median(coarse)
    sx = cn.call_sex(coarse, gmed, PAR, P)
    return cn.normalised(coarse, gmed, sx["ploidy_sex"], PAR), sx["sex_call"]


@pytest.mark.parametrize("x,y,call,ploidy", [
    (0.49, 0.19, "male", "male"),          # typical panel male (survey medians)
    (0.97, 0.0001, "female", "female"),    # typical panel female
    (0.66, 0.0004, "ambiguous", "male"),    # low-chrX female (X loss?): out of the reference, legacy ploidy
    (0.56, 0.0001, "ambiguous", "male"),   # chrX male-like but no chrY: legacy ploidy (male)
    (0.97, 0.18, "ambiguous", "female"),   # XXY-like
    (0.50, 0.045, "ambiguous", "male"),    # male with low chrY (mosaic loss)
    (1.26, 0.0004, "ambiguous", "female"),  # chrX above the female range
    (1.17, 0.0001, "ambiguous", "female"),  # NA19238 (LRSC R10): chrX-wide gains; matched L01 in the first design check
    (0.61, 0.19, "ambiguous", "male"),      # male with a chrX above the male range
    (0.49, float("nan"), "male", "male"),  # no chrY track: chrX decides
    (float("nan"), 0.2, "unknown", "female"),
])
def test_classify_sex(x, y, call, ploidy):
    assert cn.classify_sex(x, y, P) == (call, ploidy)


def test_sex_ratios_and_call():
    rng = np.random.default_rng(11)
    for sex in ("male", "female"):
        coarse = cn.to_coarse(_fine(rng, 30, sex))
        sx = cn.call_sex(coarse, cn.genome_median(coarse), PAR, P)
        assert sx["sex_call"] == sex
        assert sx["x_ratio"] == pytest.approx(0.5 if sex == "male" else 1.0, abs=0.05)
        assert (sx["y_ratio"] > 0.4) if sex == "male" else (sx["y_ratio"] == 0.0)


def test_chrx_reference_per_sex_and_ambiguous_left_out():
    rng = np.random.default_rng(12)
    members = [_member(_fine(rng, 30, "male" if i % 2 else "female")) for i in range(8)]
    amb = _fine(rng, 30, "male", y_lost=True)          # chrX male-like, chrY without depth: ambiguous
    members.append(_member(amb))
    assert members[-1][1] == "ambiguous"
    pon = cn.Pon.build(members, PAR)
    assert pon.meta["sex_ambiguous"] == 1 and pon.meta["genomes"] == 9
    assert set(pon.med) == {"chr1", "chrX|male", "chrX|female", "chrY|male"}
    assert pon.n["chr1"].min() == 9                     # autosomes: every member
    assert pon.n["chrX|male"].max() == 4 and pon.n["chrX|female"].max() == 4 and pon.n["chrY|male"].max() == 4


def test_male_specific_depth_is_absorbed_by_the_male_reference():
    """The panel-v1.2 failure mode: a chrX stretch where males read 2 copies. Against the
    male-only reference a male shows no gain there; a female still calls her real gain."""
    rng = np.random.default_rng(13)
    extra = (2000, 2300)
    pon = cn.Pon.build([_member(_fine(rng, 30, "male" if i % 2 else "female", male_extra=extra))
                        for i in range(16)], PAR)
    male = _fine(rng, 30, "male", male_extra=extra)
    segs, _ = cn.segment_genome(male, "male", PAR, pon, P)
    assert not [s for s in segs if s["chrom"] == "chrX"]
    fem = _fine(rng, 30, "female", [("chrX", 1000, 1300, 3)])
    segs, _ = cn.segment_genome(fem, "female", PAR, pon, P)
    assert [(s["chrom"], s["state"], s["ploidy"]) for s in segs] == [("chrX", 3, 2)]


def test_legacy_pooled_pon_is_still_read():
    med = {"chr1": np.ones(400, np.float32), "chrX": np.full(400, 0.9, np.float32)}
    rel = {k: np.full(400, 0.05, np.float32) for k in med}
    n = {k: np.full(400, 10, np.int32) for k in med}
    pon = cn.Pon(med, rel, n, {"stratum": "old"})
    m, _, k = pon.reference("chrX", 400, "male")
    assert m[0] == pytest.approx(0.9) and k[0] == 10


def _seg(chrom, s, e, state, ploidy):
    return dict(chrom=chrom, start=s, end=e, state=state, ploidy=ploidy)


def test_public_counting_is_sex_matched_on_chrx_only():
    cand = {"segs": [_seg("chrX", 1_000_000, 1_500_000, 3, 2), _seg("chr1", 5_000_000, 5_400_000, 3, 2)]}
    others = [
        ("m1", "M1", [_seg("chrX", 1_050_000, 1_480_000, 2, 1), _seg("chr1", 5_000_000, 5_400_000, 3, 2)], "male"),
        ("a1", "A1", [_seg("chrX", 1_000_000, 1_500_000, 3, 2)], "ambiguous"),
        ("f1", "F1", [_seg("chrX", 990_000, 1_490_000, 3, 2)], "female"),
    ]
    pub, counts = cnl.cn_public(cand, others, sex="female")
    assert counts == [1, 1] and pub                     # chrX: f1 only; chr1: m1 (autosome, any sex)
    pub, counts = cnl.cn_public(cand, others[:2], sex="female")
    assert counts == [0, 1] and not pub                 # the male's chrX gain no longer makes it public
    assert cnl.cn_public(cand, others[:2], sex=None)[1] == [2, 1]   # old rule, for comparison
    assert cnl.cn_public(cand, [o[:3] for o in others])[1] == [3, 1]  # tables without a sex column


def test_cn_linkage_reads_the_sex_column(tmp_path):
    import cn_linkage
    p = tmp_path / "panel.tsv"
    p.write_text("key\tindividual_id\tstratum\tsex\tchrom\tstart\tend\tstate\tploidy\n"
                 "c/a\tA\tT\tmale\tchrX\t1\t100\t2\t1\n"
                 "c/b\tB\tT\tambiguous\tchr1\t1\t100\t3\t2\n"
                 "c/s\tSELF\tT\tfemale\tchr1\t1\t100\t3\t2\n")
    o = sorted(cn_linkage.others(str(p), "SELF"))
    assert [(x[0], x[3]) for x in o] == [("c/a", "male"), ("c/b", "ambiguous")]
