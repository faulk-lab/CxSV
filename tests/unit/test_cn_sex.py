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
    v12 = {**cnl.CN_DEFAULTS, "public_min_genomes_sexchrom": 1}     # panel-v1.2-r2 rule
    pub, counts = cnl.cn_public(cand, others, v12, sex="female")
    assert counts == [1, 1] and pub                     # chrX: f1 only; chr1: m1 (autosome, any sex)
    pub, counts = cnl.cn_public(cand, others[:2], v12, sex="female")
    assert counts == [0, 1] and not pub                 # the male's chrX gain no longer makes it public
    assert cnl.cn_public(cand, others[:2], v12, sex=None)[1] == [2, 1]   # old rule, for comparison
    assert cnl.cn_public(cand, [o[:3] for o in others], v12)[1] == [3, 1]  # tables without a sex column
    # panel-v1.3 D2 (default): one sex-matched individual is not enough on chrX
    pub, counts = cnl.cn_public(cand, others, sex="female")
    assert counts == [1, 1] and not pub


def test_sexchrom_needs_two_sex_matched_individuals_autosomes_one():
    P = cnl.CN_DEFAULTS
    assert P["public_min_genomes"] == 1 and P["public_min_genomes_sexchrom"] == 2
    assert cnl.min_genomes("chrX") == cnl.min_genomes("chrY") == 2 and cnl.min_genomes("chr7") == 1
    x = {"segs": [_seg("chrX", 1_000_000, 1_500_000, 3, 2)]}
    a = {"segs": [_seg("chr2", 1_000_000, 1_500_000, 3, 2)]}
    one = [("c/f1", "F1", [_seg("chrX", 1_000_000, 1_500_000, 3, 2), _seg("chr2", 1_000_000, 1_500_000, 3, 2)], "female", True)]
    same_ind = one + [("d/f1", "F1", [_seg("chrX", 1_000_000, 1_500_000, 3, 2)], "female", True)]
    two = one + [("c/f2", "F2", [_seg("chrX", 1_010_000, 1_490_000, 3, 2)], "female", True)]
    assert cnl.cn_public(a, one, sex="female") == (True, [1])          # autosome: one individual
    assert cnl.cn_public(x, one, sex="female") == (False, [1])
    assert cnl.cn_public(x, same_ind, sex="female") == (False, [1])    # two libraries, one individual
    assert cnl.cn_public(x, two, sex="female") == (True, [2])
    males = [("c/m%d" % i, "M%d" % i, [_seg("chrX", 1_000_000, 1_500_000, 2, 1)], "male", True) for i in range(3)]
    assert cnl.cn_public(x, males, sex="female") == (False, [0])       # not sex-matched
    y = {"segs": [_seg("chrY", 3_000_000, 3_500_000, 0, 1)]}
    ym = [("c/m%d" % i, "M%d" % i, [_seg("chrY", 3_000_000, 3_500_000, 0, 1)], "male", True) for i in range(2)]
    assert cnl.cn_public(y, ym[:1], sex="male") == (False, [1]) and cnl.cn_public(y, ym, sex="male") == (True, [2])


# pon_members.tsv rows of panel-v1.2-r2 (calibration of D2, docs/CN_LINKAGE.md section 2b)
PANEL_LIBS = [
    dict(key="1kGP_LRSC_ONT_CHM13/NA19238", individual_id="NA19238", sex_call="ambiguous", x_ratio=1.171, y_ratio=0.0001),
    dict(key="1kGP_ONT_1019_CHM13/NA19238", individual_id="NA19238", sex_call="female", x_ratio=1.024, y_ratio=0.0001),
    dict(key="1kGP_LRSC_ONT_CHM13/HG03499", individual_id="HG03499", sex_call="female", x_ratio=0.856, y_ratio=0.0),
    dict(key="1kGP_ONT_1019_CHM13/HG03499", individual_id="HG03499", sex_call="female", x_ratio=0.970, y_ratio=0.0),
    dict(key="1kGP_LRSC_ONT_CHM13/HG03548", individual_id="HG03548", sex_call="female", x_ratio=0.953, y_ratio=0.0),
    dict(key="1kGP_ONT_1019_CHM13/HG03548", individual_id="HG03548", sex_call="female", x_ratio=1.016, y_ratio=0.0),
    dict(key="1kGP_ONT_1019_CHM13/HG00096", individual_id="HG00096", sex_call="ambiguous", x_ratio=1.26, y_ratio=0.0),
    dict(key="a/M1", individual_id="M1", sex_call="male", x_ratio=0.49, y_ratio=0.19),
    dict(key="b/M1", individual_id="M1", sex_call="male", x_ratio=0.50, y_ratio=0.32),
    dict(key="a/M2", individual_id="M2", sex_call="male", x_ratio=0.49, y_ratio="nan"),
    dict(key="b/M2", individual_id="M2", sex_call="male", x_ratio=0.50, y_ratio=0.20),
]


def test_sexchrom_discordant_libraries():
    d = cnl.sexchrom_discordant(PANEL_LIBS)
    assert sorted(d) == ["1kGP_LRSC_ONT_CHM13/HG03499", "1kGP_LRSC_ONT_CHM13/NA19238", "1kGP_ONT_1019_CHM13/HG03499",
                         "1kGP_ONT_1019_CHM13/NA19238", "a/M1", "b/M1"]
    r = d["1kGP_ONT_1019_CHM13/NA19238"]
    assert "sex call ambiguous" in r and "chrX ratio 1.171" in r and "> 0.1" in r
    assert "chrY ratio" in d["a/M1"]                   # chrY 0.19 vs 0.32
    # HG03548 (|diff| 0.063), a single-library ambiguous individual and a missing chrY ratio: kept
    assert not any(k.endswith(("HG03548", "HG00096")) or "M2" in k for k in d)
    loose = {**cnl.CN_DEFAULTS, "sexchrom_discordance_x": 0.15, "sexchrom_discordance_y": 0.2}
    assert sorted(cnl.sexchrom_discordant(PANEL_LIBS, loose)) == ["1kGP_LRSC_ONT_CHM13/NA19238",
                                                                   "1kGP_ONT_1019_CHM13/NA19238"]   # sex call only


def test_discordant_library_counts_on_autosomes_not_on_sex_chromosomes():
    """PC L01 (female chrX DUP-TRP) under D2: NA19238's chrX gain plus one other female."""
    x = {"segs": [_seg("chrX", 1_000_000, 1_500_000, 3, 2), _seg("chr3", 9_000_000, 9_500_000, 1, 2)]}
    others = [("lrsc/NA19238", "NA19238", [_seg("chrX", 900_000, 1_600_000, 3, 2), _seg("chr3", 9_000_000, 9_500_000, 1, 2)],
               "female", False),
              ("c/F2", "F2", [_seg("chrX", 1_000_000, 1_500_000, 3, 2)], "female", True)]
    pub, counts = cnl.cn_public(x, others, sex="female")
    assert counts == [1, 1] and not pub               # chrX: F2 only (< 2); chr3: NA19238 still counts
    others[0] = others[0][:4] + (True,)
    assert cnl.cn_public(x, others, sex="female") == (True, [2, 1])


def test_cn_linkage_reads_the_sex_column(tmp_path):
    import cn_linkage
    p = tmp_path / "panel.tsv"
    p.write_text("key\tindividual_id\tstratum\tsex\tchrom\tstart\tend\tstate\tploidy\n"
                 "c/a\tA\tT\tmale\tchrX\t1\t100\t2\t1\n"
                 "c/b\tB\tT\tambiguous\tchr1\t1\t100\t3\t2\n"
                 "c/s\tSELF\tT\tfemale\tchr1\t1\t100\t3\t2\n")
    o = sorted(cn_linkage.others(str(p), "SELF"))
    assert [(x[0], x[3]) for x in o] == [("c/a", "male"), ("c/b", "ambiguous")]


def test_cn_linkage_sexchrom_exclusion_from_pon_members(tmp_path):
    import cn_linkage
    head = "key\tindividual_id\tstratum\tsex\tsex_call\tx_ratio\ty_ratio\tgenome_median\tin_pon\tin_sex_reference"
    old = tmp_path / "v12.tsv"                          # panel-v1.2-r2 table: recomputed from the ratios
    old.write_text(head + "\n" + "".join(
        f"{r['key']}\t{r['individual_id']}\tT\tfemale\t{r['sex_call']}\t{r['x_ratio']}\t{r['y_ratio']}\t30\tyes\tyes\n"
        for r in PANEL_LIBS[:6]))
    ex = cn_linkage.sexchrom_excluded(str(old))
    assert sorted(ex) == sorted(r["key"] for r in PANEL_LIBS[:4])
    new = tmp_path / "v13.tsv"                          # panel-v1.3 table: the recorded column decides
    new.write_text(head + "\tsexchrom_public\n" + "c/a\tA\tT\tfemale\tfemale\t0.97\t0\t30\tyes\tyes\tyes\n"
                   "c/b\tB\tT\tfemale\tambiguous\t1.2\t0\t30\tyes\tno: sex ambiguous\tno: discordant libraries of B: x\n")
    assert cn_linkage.sexchrom_excluded(str(new)) == {"c/b": "discordant libraries of B: x"}
    assert cn_linkage.sexchrom_excluded("none") == {}
    panel = tmp_path / "panel.tsv"
    panel.write_text("key\tindividual_id\tstratum\tsex\tchrom\tstart\tend\tstate\tploidy\n"
                     "c/a\tA\tT\tfemale\tchrX\t1\t100\t3\t2\nc/b\tB\tT\tfemale\tchrX\t1\t100\t3\t2\n")
    o = sorted(cn_linkage.others(str(panel), "SELF", cn_linkage.sexchrom_excluded(str(new))))
    assert [(x[0], x[4]) for x in o] == [("c/a", True), ("c/b", False)]


def test_cn_pon_records_sexchrom_public(tmp_path):
    import gzip
    import json
    import subprocess
    rng = np.random.default_rng(11)
    rows = ["key\tindividual_id\tstratum\tdepth"]
    libs = [("c/f1", "F1", "female", ()), ("d/f1", "F1", "female", (("chrX", 200, 4000, 3),)),   # chrX ~1.5: discordant
            ("c/f2", "F2", "female", ()), ("d/f2", "F2", "female", ()), ("c/m1", "M1", "male", ())]
    for key, ind, sex, ev in libs:
        p = tmp_path / (key.replace("/", "_") + ".bed.gz")
        with gzip.open(p, "wt") as fh:
            for c, a in _fine(rng, 30, sex, ev).items():
                fh.writelines(f"{c}\t{i * 1000}\t{i * 1000 + 1000}\t{v:.2f}\n" for i, v in enumerate(a))
        rows.append(f"{key}\t{ind}\tT\t{p}")
    (tmp_path / "members.tsv").write_text("\n".join(rows) + "\n")
    (tmp_path / "p.json").write_text(json.dumps({"cn_params": {"segmentation": dict(cn.DEFAULTS),
                                                                "linkage": {"sexchrom_discordance_x": 0.10,
                                                                            "sexchrom_discordance_y": 0.10},
                                                                "par": {"chrX": [[0, 200000]], "chrY": [[0, 100000]]}}}))
    subprocess.run([sys.executable, ROOT / "workflow" / "scripts" / "cn_pon.py", "--members", tmp_path / "members.tsv",
                    "--params", tmp_path / "p.json", "--out-dir", tmp_path / "pon", "--threads", "1"], check=True)
    lines = (tmp_path / "pon" / "pon_members.tsv").read_text().splitlines()
    assert lines[0].split("\t")[-1] == "sexchrom_public"
    col = {l.split("\t")[0]: l.split("\t")[-1] for l in lines[1:]}
    assert col["c/f1"].startswith("no: discordant libraries of F1") and col["d/f1"] == col["c/f1"]
    assert col["c/f2"] == col["d/f2"] == col["c/m1"] == "yes"
