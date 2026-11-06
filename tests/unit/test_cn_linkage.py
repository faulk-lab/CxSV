"""Roadmap item 15 Phase 1 prototype (calibration/cn_linkage): segmentation and
copy-number linkage on synthetic depth. numpy only (CI installs no pandas)."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "calibration" / "cn_linkage"))
import cnlib   # noqa: E402
import cnlink  # noqa: E402

NB = 4000          # 1 kb bins per synthetic contig (4 Mb)
CONTIGS = ["chr1", "chrX"]


def _genome(rng, cov, sex, events=()):
    """1 kb Poisson depth; events: (contig, start_kb, end_kb, copies)."""
    fine = {}
    for c in CONTIGS:
        pl = 2 if (c != "chrX" or sex == "female") else 1
        lam = np.full(NB, cov * pl / 2.0)
        for ec, s, e, k in events:
            if ec == c:
                lam[s:e] = cov * k / 2.0
        fine[c] = rng.poisson(lam).astype(np.float32)
    return fine


def _pon(rng, cov=30, n=8):
    genomes, coarse, sexes, fams = [], {}, {}, {}
    for i in range(n):
        g, sex = f"n{i}", ("male" if i % 2 else "female")
        fine = _genome(rng, cov, sex)
        cs = cnlib.to_coarse(fine)
        coarse[g] = (cs, float(np.median(np.concatenate([cs["chr1"]]))))
        genomes.append(g); sexes[g] = sex; fams[g] = g
    return cnlib.PanelOfNormals(genomes, coarse, sexes, fams, {})


def _segment(fine, sex, pon):
    old = cnlib.AUTOSOMES
    cnlib.AUTOSOMES = ["chr1"]
    try:
        return cnlib.segment_genome(fine, sex, {}, pon, contigs=CONTIGS)[0]
    finally:
        cnlib.AUTOSOMES = old


def test_read_bins_without_pandas(tmp_path, monkeypatch):
    import gzip
    p = tmp_path / "d.regions.bed.gz"
    with gzip.open(p, "wt") as fh:
        for i in range(25):
            fh.write(f"chr1\t{i*1000}\t{i*1000+1000}\t{i}.25\n")
        fh.write("chrUn\t0\t1000\t3\n")
    monkeypatch.setattr(cnlib, "pd", None)
    b = cnlib.read_bins(str(p))
    assert list(b) == ["chr1"] and len(b["chr1"]) == 25 and b["chr1"][3] == pytest.approx(3.25)


def test_null_genome_has_no_segments():
    rng = np.random.default_rng(1)
    pon = _pon(rng)
    assert _segment(_genome(rng, 30, "female"), "female", pon) == []


@pytest.mark.parametrize("sex,base", [("male", 1), ("female", 2)])
def test_dup_trp_segments_and_class(sex, base):
    """DUP (+1) - TRP (+2) - DUP (+1) on chrX (male 1-2-3-2, female 2-3-4-3)."""
    rng = np.random.default_rng(2)
    pon = _pon(rng)
    ev = [("chrX", 1000, 1300, base + 1), ("chrX", 1300, 1500, base + 2), ("chrX", 1500, 1800, base + 1)]
    segs = _segment(_genome(rng, 30, sex, ev), sex, pon)
    states = [s["state"] for s in segs]
    assert states == [base + 1, base + 2, base + 1]
    for s, (_, a, b, _) in zip(segs, ev):
        assert abs(s["start"] - a * 1000) <= 10_000 and abs(s["end"] - b * 1000) <= 10_000   # attach tol
    inv = dict(id="r1", chrom="chrX", pos=1_300_000, chr2="chrX", end=1_500_000, svtype="INV",
               svlen=200_000, supp="", qual=None, inverted=True, line="")
    cands = cnlink.candidates(segs, [inv])
    assert len(cands) == 1
    assert cands[0]["cls"] == "DUP-TRP/INV-DUP" and cands[0]["decision"] == "locus"


def test_single_dup_is_simple_cnv():
    seg = dict(chrom="chr1", start=1_000_000, end=1_300_000, state=3, ploidy=2)
    dup = dict(id="d", chrom="chr1", pos=1_000_500, chr2="chr1", end=1_299_500, svtype="DUP",
               svlen=299_000, supp="", qual=None, inverted=False, line="")
    c = cnlink.candidates([seg], [dup])
    assert c[0]["decision"] == "simple_cnv"


def test_cn_only_and_public_filter():
    a = dict(chrom="chr1", start=1_000_000, end=1_200_000, state=3, ploidy=2)
    b = dict(chrom="chr1", start=1_200_000, end=1_400_000, state=1, ploidy=2)
    assert cnlink.candidates([a, b], [])[0]["decision"] == "cn_only"
    bnd = dict(id="x", chrom="chr1", pos=1_200_100, chr2="chr1", end=1_399_900, svtype="BND",
               svlen=0, supp="", qual=None, inverted=False, line="")
    assert cnlink.candidates([a, b], [bnd])[0]["decision"] == "locus"
    other = [("g2", "f2", [dict(a), dict(b)])]
    c = cnlink.candidates([a, b], [bnd], others=other)[0]
    assert c["decision"] == "public_cn" and c["public_counts"] == [1, 1]


def test_inverted_bnd_alt_parsing(tmp_path):
    v = tmp_path / "m.vcf"
    v.write_text("##fileformat=VCFv4.2\n"
                 "chrX\t100\tb1\tN\tN]chrX:5000]\t.\tPASS\tSVTYPE=BND\n"
                 "chrX\t100\tb2\tN\tN[chrX:5000[\t.\tPASS\tSVTYPE=BND\n"
                 "chrX\t100\ti1\tN\t<INV>\t.\tPASS\tSVTYPE=INV;END=900\n")
    r = {x["id"]: x for x in cnlink.read_records(str(v))}
    assert r["b1"]["inverted"] and not r["b2"]["inverted"] and r["i1"]["inverted"]
    assert r["b1"]["end"] == 5000 and r["i1"]["end"] == 900


def test_nested_scan_finds_short_compressed_triplication():
    """A 120 kb TRP observed at half a copy above the DUP (as PC L01): the HMM keeps one
    gain segment, the nested scan splits it into DUP / TRP / DUP."""
    rng = np.random.default_rng(3)
    pon = _pon(rng)
    ev = [("chr1", 1000, 1500, 3), ("chr1", 1500, 1620, 3.5), ("chr1", 1620, 1700, 3)]
    segs = _segment(_genome(rng, 38, "female", ev), "female", pon)
    assert [s["state"] for s in segs] == [3, 4, 3]
    assert abs(segs[1]["start"] - 1_500_000) <= 5000 and abs(segs[1]["end"] - 1_620_000) <= 5000
    assert segs[1]["nested"]


def test_nested_scan_quiet_on_plain_gain():
    rng = np.random.default_rng(4)
    pon = _pon(rng)
    segs = _segment(_genome(rng, 30, "female", [("chr1", 1000, 1600, 3)]), "female", pon)
    assert [s["state"] for s in segs] == [3]


def test_attach_widens_over_masked_bins_and_inv_other_end_slack():
    # gain 1,000,000-1,300,000 whose start sits after 30 kb of masked bins
    g = dict(chrom="chr1", start=1_000_000, end=1_300_000, state=3, ploidy=2,
             start_mask_before=30_000, start_mask_after=0, end_mask_before=0, end_mask_after=0)
    t = dict(chrom="chr1", start=1_300_000, end=1_400_000, state=4, ploidy=2)
    inv = dict(id="i", chrom="chr1", pos=975_000, chr2="chr1", end=1_250_000, svtype="INV", svlen=275_000,
               supp="", qual=None, inverted=True, line="")
    c = cnlink.candidates([g, t], [inv])[0]
    assert len(c["records"]) == 1 and c["cls"] == "DUP-TRP/INV-DUP"
    g2 = dict(g, start_mask_before=0)
    assert cnlink.candidates([g2, t], [inv])[0]["records"] == []
    # one end at the TRP end, other end 60 kb past the candidate (outer DUP undetected)
    inv2 = dict(inv, pos=1_400_500, end=1_460_000)
    assert len(cnlink.candidates([g2, t], [inv2])[0]["records"]) == 1
