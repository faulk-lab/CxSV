"""External SV catalogue conversion (acquisition/catalogues): schema rules,
per-source adapters and breakpoint liftover bookkeeping. Pure Python: no
bcftools, liftOver or catalogue data needed."""
import io
import sys

import pytest

from conftest import REPO

sys.path.insert(0, str(REPO / "acquisition" / "catalogues"))

import adapters  # noqa: E402
import junction_sites as js  # noqa: E402
import liftover  # noqa: E402

# hs1 .fai order is by length: chr7 then chrX then chr9, chr8
ORDER = js.ContigOrder(["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
                        "chrX", "chr9", "chr8", "chr11", "chr10", "chr12", "chr13",
                        "chr14", "chr15", "chr16", "chr17", "chr18", "chr20", "chrY",
                        "chr19", "chr22", "chr21", "chrM"])


def write(sites, groups=()):
    fh = io.StringIO()
    counts = js.write_sites(fh, sites, "SRC", "v1", list(groups), ORDER)
    rows = [l.split("\t") for l in fh.getvalue().splitlines()]
    return rows, counts


def test_header_is_shared_schema():
    rows, _ = write([], groups=["AFR"])
    assert rows[0] == ["chrom1", "pos1", "chrom2", "pos2", "svtype", "svlen", "ac",
                       "an", "af", "source", "version", "af_AFR"]


def test_bnd_oriented_in_fai_order_not_lexical():
    # chr8 comes AFTER chr9 in hs1 .fai order, and chrX before chr9
    s = js.make_site("chr8", 100, "chr9", 500, "BND")
    rows, _ = write([s])
    assert rows[1][:4] == ["chr9", "500", "chr8", "100"]
    s = js.make_site("chr9", 5, "chrX", 7, "BND")
    rows, _ = write([s])
    assert rows[1][:4] == ["chrX", "7", "chr9", "5"]


def test_missing_values_and_af_from_ac_an():
    rows, _ = write([js.make_site("chr1", 10, "chr1", 100, "DEL", 90, ac=1, an=4)])
    assert rows[1][6:9] == ["1", "4", "0.25"]
    rows, _ = write([js.make_site("chr1", 10, "chr1", 10, "INS", 60)])
    assert rows[1][6:9] == [".", ".", "."]


def test_validate_rejects_bad_sites():
    with pytest.raises(ValueError):
        write([js.make_site("chr1", 10, "chr2", 100, "DEL", 90)])
    with pytest.raises(ValueError):
        write([js.make_site("chr1", 10, "chr1", 12, "INS", 60)])
    with pytest.raises(ValueError):
        write([js.make_site("chr1", 10, "chr1", 100, "CPX", 90)])


def test_contig_and_type_normalisation():
    assert js.normalise_contig("1") == "chr1"
    assert js.normalise_contig("MT") == "chrM"
    assert js.normalise_contig("CHM13#0#chr5") == "chr5"
    assert js.normalise_contig("chrUn_KI270302v1") is None
    assert js.normalise_svtype("<INS:ME:ALU>") == "INS"
    assert js.normalise_svtype("CTX") == "BND"
    assert js.normalise_svtype("TRA") == "BND"
    assert js.normalise_svtype("CNV") is None


def test_sequence_resolved_alleles():
    ref = "A" + "C" * 60
    assert js.classify_sequence_allele(100, ref, "A") == ("DEL", 100, 160, 60)
    assert js.classify_sequence_allele(100, "A", "A" + "G" * 75) == ("INS", 100, 100, 75)
    assert js.classify_sequence_allele(100, "A", "A" + "G" * 49) is None
    inv = js.classify_sequence_allele(100, "ACGT" * 20, "TGCA" * 20, inv_flag=True)
    assert inv == ("INV", 100, 179, 80)


def test_symbolic_records():
    s, why = adapters.vcf_sv_sites("chr1", 2708, "C", "<DEL>", "DEL", None, -62, ac=2, an=130)
    assert why is None and (s[0]["pos2"], s[0]["svlen"], s[0]["svtype"]) == (2770, 62, "DEL")
    s, _ = adapters.vcf_sv_sites("chr1", 2708, "C", "<INS>", "INS", None, 169)
    assert (s[0]["pos1"], s[0]["pos2"], s[0]["svlen"]) == (2708, 2708, 169)
    s, _ = adapters.vcf_sv_sites("chr2", 50, "N", "N[chr1:900[", "BND")
    assert (s[0]["chrom2"], s[0]["pos2"]) == ("chr1", 900)
    assert adapters.vcf_sv_sites("chr1", 1, "C", "<DEL>", "DEL", None, -60, ac=0)[1] == "ac_zero"
    assert adapters.vcf_sv_sites("chr1", 1, "C", "<DEL>", "DEL", None, -20)[1] == "below_min_len"


def test_gnomad_record_types():
    kw = dict(chr2=None, end2=None, cpx_intervals=None, ac="3", an="100", af="0.03", pop_af={"AFR": 0.1})
    s, why = adapters.gnomad_sites("chr1", 100, "<DEL>", "PASS", "DEL", "900", "800", **kw)
    assert why is None and s[0]["pop_af"] == {"AFR": 0.1} and s[0]["pos2"] == 900
    assert adapters.gnomad_sites("chr1", 100, "<DEL>", "HIGH_NCR", "DEL", "900", "800", **kw)[1] == "filter_HIGH_NCR"
    assert adapters.gnomad_sites("chr1", 100, "<CNV>", "PASS", "CNV", "900", "800", **kw)[1] == "svtype_CNV"
    kw.update(chr2="chr5", end2="7000")
    s, _ = adapters.gnomad_sites("chr1", 100, "<CTX>", "PASS", "CTX", "101", "-1", **kw)
    assert (s[0]["svtype"], s[0]["chrom2"], s[0]["pos2"], s[0]["svlen"]) == ("BND", "chr5", 7000, None)
    kw.update(chr2="chr1", end2=None, cpx_intervals="DUP_chr1:100-2000,INV_chr1:100-2000")
    s, _ = adapters.gnomad_sites("chr1", 100, "<CPX>", "PASS", "CPX", "2000", "1900", **kw)
    assert [x["svtype"] for x in s] == ["BND", "DUP", "INV"]


def test_hgsvc3_cpx_tiling():
    s, why, tiled = adapters.hgsvc3_cpx_sites("chr1", 1000, 1300, 300, "DEL(100):NML(100):INVDUP(100)", 2, 124)
    assert tiled and why is None
    assert [(x["svtype"], x["pos1"], x["pos2"]) for x in s] == [
        ("BND", 1000, 1300), ("DEL", 1000, 1100), ("INV", 1200, 1300), ("DUP", 1200, 1300)]
    s, _, tiled = adapters.hgsvc3_cpx_sites("chr1", 1000, 1300, 300, "DUP(282):DUP(179)", 1, 124)
    assert not tiled and [x["svtype"] for x in s] == ["BND"]


def test_apply_lift_rules():
    sites = [js.make_site("chr1", 10, "chr1", 100, "DEL", 90),
             js.make_site("chr1", 20, "chr1", 20, "INS", 60),
             js.make_site("chr1", 30, "chr1", 300, "DUP", 270),
             js.make_site("chr1", 40, "chr2", 400, "BND"),
             js.make_site("chr1", 50, "chr1", 500, "INV", 450)]
    mapped = {("chr1", 10): ("chr1", 1010), ("chr1", 100): ("chr1", 1100),
              ("chr1", 20): ("chr1", 1020),
              ("chr1", 30): ("chr1", 1030), ("chr1", 300): ("chr3", 5),
              ("chr1", 40): ("chr1", 1040), ("chr2", 400): ("chr7", 9),
              ("chr1", 50): ("chr1", 1050)}
    kept, stats = liftover.apply_lift(sites, mapped)
    assert [k["svtype"] for k in kept] == ["DEL", "INS", "BND"]
    assert kept[1]["pos2"] == 1020
    assert stats["DUP"]["lost_split_contig"] == 1
    assert stats["INV"]["lost_unmapped_end"] == 1
