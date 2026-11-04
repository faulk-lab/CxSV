"""Curated literature set (roadmap item 11): junction model, parsers, liftover
bookkeeping, the literature_match annotation and the add-on builder's cutting.
Pure Python (no liftOver, pysam or source files needed)."""
import io
import sys
import zipfile

import pytest

from conftest import REPO

sys.path.insert(0, str(REPO / "acquisition" / "literature"))
sys.path.insert(0, str(REPO / "acquisition" / "catalogues"))
sys.path.insert(0, str(REPO / "tests" / "literature"))
sys.path.insert(0, str(REPO / "tests" / "positive_control"))

import literature_set as L          # noqa: E402
import build_literature_set as B    # noqa: E402
import match_literature as M        # noqa: E402
import score_literature_set as S    # noqa: E402
import literature_addons as A       # noqa: E402
import cxsv_scoring as sc           # noqa: E402

NAMES = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chrX", "chr9", "chr8", "chr11", "chr10",
         "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr20", "chrY", "chr19", "chr22",
         "chr21", "chrM"]
ORDER = {c: i for i, c in enumerate(NAMES)}


# ─── junction model ──────────────────────────────────────────────────────────
@pytest.mark.parametrize("s1,s2,want", [("+", "-", "DEL"), ("-", "+", "DUP"), ("+", "+", "INV"), ("-", "-", "INV")])
def test_svtype_from_sides(s1, s2, want):
    assert L.svtype_of("chr1", 100, s1, "chr1", 900, s2, ORDER) == want


def test_sides_travel_with_ends_when_reordered():
    # given high end first: (900 '-', 100 '+') is the DEL-like join (100 +, 900 -)
    j = L.make_junction("E", 1, "chr1", 900, "-", "chr1", 100, "+", "bp", ORDER)
    assert (j["pos1"], j["side1"], j["pos2"], j["side2"], j["svtype"]) == (100, "+", 900, "-", "DEL")
    # BND ends follow hs1 .fai order: chrX sorts before chr9
    j = L.make_junction("E", 2, "chr9", 5, "+", "chrX", 7, "-", "bp", ORDER)
    assert (j["chrom1"], j["chrom2"], j["svtype"]) == ("chrX", "chr9", "BND")
    with pytest.raises(ValueError):
        L.make_junction("E", 3, "chr1", 1, "?", "chr1", 2, "+", "bp", ORDER)


def test_derivative_junctions_and_reference_adjacency():
    # der: pter..100 (+) | 301..400 (-) | 101..200 (+) | 401..qter (+)  = an inverted, displaced segment
    frags = [("A", "chr1", "pter", 100, "+"), ("B", "chr1", 301, 400, "-"), ("C", "chr1", 101, 200, "+"),
             ("D", "chr1", 401, "qter", "+")]
    js = L.derivative_junctions(frags)
    assert [(j[3], j[4], j[6], j[7]) for j in js] == [(100, "+", 400, "+"), (301, "-", 101, "-"),
                                                      (200, "+", 401, "-")]
    assert [L.svtype_of(*j[2:], ORDER) for j in js] == ["INV", "INV", "DEL"]
    # consecutive reference-adjacent fragments are not a junction
    assert L.derivative_junctions([("A", "chr1", 1, 100, "+"), ("B", "chr1", 101, 200, "+")]) == []
    assert L.derivative_junctions([("A", "chr1", 101, 200, "-"), ("B", "chr1", 1, 100, "-")]) == []


def test_fix_derivative_boundaries_reattaches_telomeric_fragment():
    ders = {"der(7)": [("Z", "chr7", "pter", 1000, "+"), ("AD", "chr2", 50, 90, "-")],
            "der(8)": [("AE", "chr7", 1001, "qter", "+"), ("AF", "chr8", "pter", 500, "+")]}
    fixed, moved = L.fix_derivative_boundaries(ders)
    assert [f[0] for f in fixed["der(7)"]] == ["Z", "AD", "AE"]
    assert [f[0] for f in fixed["der(8)"]] == ["AF"]
    assert moved == ["AE: der(8) -> der(7)"]
    ok = {"der(1)": [("A", "chr1", "pter", 10, "+"), ("B", "chr1", 20, "qter", "+")]}
    assert L.fix_derivative_boundaries(ok) == (ok, [])


def test_breakends_merge_within_distance():
    js = [dict(chrom1="chr1", pos1=100, chrom2="chr1", pos2=5000), dict(chrom1="chr1", pos1=250, chrom2="chr2", pos2=9)]
    assert L.breakends(js, 300) == 3


# ─── parsers ─────────────────────────────────────────────────────────────────
def make_xlsx(rows):
    """Minimal .xlsx with inline strings / numbers for read_xlsx."""
    def cell(r, c, v):
        ref = f"{chr(65 + c)}{r + 1}"
        if isinstance(v, (int, float)):
            return f'<c r="{ref}"><v>{v}</v></c>'
        return f'<c r="{ref}" t="inlineStr"><is><t>{v}</t></is></c>'
    sheet = ('<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheetData>'
             + "".join(f'<row r="{i + 1}">' + "".join(cell(i, j, v) for j, v in enumerate(r) if v != "") + "</row>"
                       for i, r in enumerate(rows)) + "</sheetData></worksheet>")
    wb = ('<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
          'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets>'
          '<sheet name="Sheet1" sheetId="1" r:id="rId1"/></sheets></workbook>')
    rels = ('<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            '<Relationship Id="rId1" Target="worksheets/sheet1.xml" '
            'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet"/></Relationships>')
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        z.writestr("xl/workbook.xml", wb)
        z.writestr("xl/_rels/workbook.xml.rels", rels)
        z.writestr("xl/worksheets/sheet1.xml", sheet)
    return buf


def test_read_xlsx_and_parse_nazaryan(tmp_path):
    rows = [["Name", "Category", "Chromosome", "Junction", "Side 1", "Side 2", "Side 1:\nOrientation"],
            ["Patient 72", "Deletions only", 7, 1, 70609299, 74047984],
            ["", "", "", 2, 72423000, 74049001],
            ["P9", "x", 21, 1, 100, ""]]
    p = tmp_path / "s2.xlsx"
    p.write_bytes(make_xlsx(rows).getvalue())
    sheet = L.read_xlsx(str(p))["Sheet1"]
    assert sheet[1][0] == "Patient 72" and sheet[2][3] == "2"
    types = {("P72", "1"): dict(svtype="INV", side1="+", side2="+", basis="hgvs"),
             ("P72", "2"): dict(svtype="INV", side1="-", side2="-", basis="hgvs"),
             ("P9", "1"): dict(svtype=".", side1=".", side2=".", basis="one_sided")}
    js, skipped = L.parse_nazaryan2018(sheet, types, ORDER)
    assert [(j["junction_id"], j["chrom1"], j["pos1"], j["pos2"], j["svtype"]) for j in js] == \
        [("NP18_P72.J1", "chr7", 70609299, 74047984, "INV"), ("NP18_P72.J2", "chr7", 72423000, 74049001, "INV")]
    assert skipped == [{"event_id": "NP18_P9", "junction": 1, "reason": "one_sided"}]
    bad = dict(types)
    bad[("P72", "1")] = dict(svtype="DEL", side1="+", side2="+", basis="x")
    with pytest.raises(ValueError):
        L.parse_nazaryan2018(sheet, bad, ORDER)


def test_grochowski_model_is_the_E01_structure():
    rows = [dict(case="BAB1", dup1_start=1000, dup1_end=2000, trp_start=2000, trp_end=2500, dup2_start=2500,
                 dup2_end=3000)]
    js, cn, flags = L.parse_grochowski2024(rows, ORDER)
    assert [(j["pos1"], j["side1"], j["pos2"], j["side2"]) for j in js] == [(2500, "+", 3000, "+"),
                                                                           (1000, "-", 2000, "-")]
    assert [c["delta"] for c in cn] == [1, 2, 1] and flags == []


# ─── liftover bookkeeping ────────────────────────────────────────────────────
def test_lifted_end_flips_side_in_reversed_block():
    mapped = {"GRCh37": {("chr1", 100): ("chr1", 5000), ("chr1", 101): ("chr1", 4999),
                         ("chr1", 200): ("chr1", 7000), ("chr1", 201): ("chr1", 7001)}}
    assert B.lifted_end("GRCh37", "chr1", 100, mapped) == ("chr1", 5000, True)
    assert B.lifted_end("GRCh37", "chr1", 200, mapped) == ("chr1", 7000, False)
    assert B.lifted_end("GRCh37", "chr1", 300, mapped) is None


# ─── scoring-test VCF: sides survive the round trip through the scorer ──────
@pytest.mark.parametrize("s1,s2", [("+", "-"), ("-", "+"), ("+", "+"), ("-", "-")])
def test_synthetic_vcf_alt_encodes_sides(tmp_path, s1, s2):
    j = dict(chrom1="chr1", pos1=1000, side1=s1, chrom2="chr1", pos2=9000, side2=s2,
             svtype=L.svtype_of("chr1", 1000, s1, "chr1", 9000, s2, ORDER))
    vcf = tmp_path / "m.vcf"
    S.write_vcf(str(vcf), [j], {"chr1": 100000}, "0/1")
    rec = [ln.split("\t") for ln in open(vcf) if not ln.startswith("#")][0]
    assert sc.record_orientation(rec[0], int(rec[1]), rec[4], j["svtype"]) == {(s1, s2)}
    assert f"SVTYPE={j['svtype']}" in rec[7] and "SUPP_VEC=111" in rec[7]


def test_synthetic_depth_follows_copy_changes(tmp_path):
    bed = tmp_path / "seg.bed"
    bed.write_text("chr1\t0\t1000\tL\tflank_left\t0\nchr1\t1000\t2000\tL\tsegment\t1\n"
                   "chrX\t0\t1000\tL\tsegment\t2\n")
    out = tmp_path / "depth.tsv"
    S.synth_depth(str(bed), str(out), [dict(chrom="chr1", start=1000, end=2000, delta=-1)], "male")
    sums = [int(ln.rstrip().split("\t")[-1]) for ln in open(out)]
    assert sums == [30000, 15000, 15000]       # 30x flank; one copy lost; male chrX one copy


# ─── literature_match annotation ─────────────────────────────────────────────
def test_match_literature_levels():
    lit_j = [dict(event_id="E1", junction_id="E1.J1", svtype="INV", chrom1="chr1", pos1="1000", chrom2="chr1",
                  pos2="9000", precision="bp", lift_status="lifted"),
             dict(event_id="E1", junction_id="E1.J2", svtype="INV", chrom1="chr1", pos1="5000", chrom2="chr1",
                  pos2="20000", precision="bp", lift_status="lifted")]
    lit_e = [dict(event_id="E1", phenotype="test phenotype", region_hs1="chr1:1000-20000"),
             dict(event_id="CV1", phenotype="clinvar region", region_hs1="chr2:800-2500"),
             dict(event_id="BIG", phenotype="large region", region_hs1="chr3:1-5000000")]
    loci = [dict(Chr="chr1", Start="900", End="20100"), dict(Chr="chr1", Start="900", End="9100"),
            dict(Chr="chr2", Start="1000", End="2000"), dict(Chr="chr3", Start="1", End="5")]
    js = [dict(svtype="INV", chrom1="chr1", pos1=1200, chrom2="chr1", pos2=8800, loci={"chr1:900-20100", "chr1:900-9100"}),
          dict(svtype="INV", chrom1="chr1", pos1=5100, chrom2="chr1", pos2=19900, loci={"chr1:900-20100"}),
          dict(svtype="DEL", chrom1="chr2", pos1=1000, chrom2="chr2", pos2=2000, loci={"chr2:1000-2000"})]
    tol = {"bp": 500, "approximate": 10000, "region": 50000}
    rows = {r["Locus_ID"]: r for r in M.annotate(loci, js, lit_j, lit_e, tol, 0.75)}
    assert rows["chr1:900-20100"]["Literature_Match"] == "event"
    assert rows["chr1:900-20100"]["Literature_Events"] == "E1(2/2)"
    assert rows["chr1:900-9100"]["Literature_Match"] == "junction"
    assert rows["chr2:1000-2000"]["Literature_Match"] == "region"
    assert rows["chr3:1-5"]["Literature_Match"] == "none"
    # a junction 600 bp off misses at 500 bp tolerance
    far = [dict(js[0], pos1=1600, loci={"chr1:900-9100"})]
    # ... and the locus covers only 43% of E1's region: no region match at the default 0.5
    assert M.annotate([loci[1]], far, lit_j, lit_e, tol, 0.75)[0]["Literature_Match"] == "none"
    assert M.annotate([loci[1]], far, lit_j, lit_e, tol, 0.75, region_min_overlap=0)[0]["Literature_Match"] == "region"
    # a 4 bp locus inside a 5 Mb published region is not a region match (reciprocal overlap)
    assert rows["chr3:1-5"]["Literature_Region_Events"] == 0


def test_reciprocal_overlap():
    assert M.reciprocal(0, 100, 50, 150, 0.5)
    assert not M.reciprocal(0, 100, 60, 1000, 0.5)
    assert M.reciprocal(0, 100, 60, 1000, 0)
    assert not M.reciprocal(0, 100, 200, 300, 0)


# ─── add-on builder: windowed cutting ────────────────────────────────────────
def test_derivative_contigs_cut_long_fragments():
    W = 100
    frags = [dict(chrom="chr1", start0=0, end0=10_000, strand="+"),      # long, first: tail only
             dict(chrom="chr2", start0=500, end0=600, strand="-"),       # short: whole
             dict(chrom="chr1", start0=20_000, end0=30_000, strand="+"),  # long: head closes, tail opens
             dict(chrom="chr3", start0=0, end0=5_000, strand="-")]       # long, last: head only
    cs = A.derivative_contigs(frags, W)
    assert cs == [[("chr1", 9_900, 10_000, "+"), ("chr2", 500, 600, "-"), ("chr1", 20_000, 20_100, "+")],
                  [("chr1", 29_900, 30_000, "+"), ("chr3", 4_900, 5_000, "-")]]


def test_addon_templates_reproduce_published_junctions():
    js = [dict(chrom1="chrX", pos1=1_000_001, side1="-", chrom2="chrX", pos2=1_300_001, side2="-"),
          dict(chrom1="chrX", pos1=1_400_000, side1="+", chrom2="chrX", pos2=1_450_000, side2="+")]
    segs, paths = A.tmpl_dup_trp_inv_dup(js, 1000)
    assert segs["B"] == ["chrX", 1_000_000, 1_300_000] and segs["C"] == ["chrX", 1_300_000, 1_400_000]
    assert segs["D"] == ["chrX", 1_400_000, 1_450_000]
    js = [dict(chrom1="chr17", pos1=100, side1="+", chrom2="chr17", pos2=500, side2="+"),
          dict(chrom1="chr17", pos1=300, side1="-", chrom2="chr17", pos2=700, side2="-")]
    segs, paths = A.tmpl_del_inv_del(js, 50)
    assert segs["D1"] == ["chr17", 100, 299] and segs["I"] == ["chr17", 299, 500] and segs["D2"] == ["chr17", 500, 699]
    assert paths == {"der": ["p+", "I-", "q+"]}
