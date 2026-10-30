"""ClinGen region list (roadmap item 9): parser and liftover bookkeeping in
acquisition/resources/build_resources.py. Pure Python: liftOver is replaced by
an in-memory point map."""
import sys

from conftest import REPO

sys.path.insert(0, str(REPO / "acquisition" / "resources"))
sys.path.insert(0, str(REPO / "acquisition" / "catalogues"))

import build_resources as br  # noqa: E402

HEADER = ("#ISCA ID\tISCA Region Name\tcytoBand\tGenomic Location\tHaploinsufficiency Score\t"
          "Haploinsufficiency Description\tTriplosensitivity Score\tTriplosensitivity Description\t"
          "Date Last Evaluated\n")
TEXT = ("#ClinGen Region Curation Results\n#28 Apr,2026\n#Genomic Locations are reported on GRCh38 (hg38)\n"
        + HEADER
        + "ISCA-1\tA region\t1p36\tchr1:1001-5000\t3\tSufficient\t0\tNo evidence\t2025-01-01\n"
          "ISCA-2\tB region\tXq28\tchrX:20001-30000\t0\tNo evidence\t3\tSufficient\t2025-01-02\n"
          "ISCA-3\tC region\t7q34\ttbd\t40\tUnlikely\t\t\t2025-01-03\n"
          "ISCA-4\tD region\t2p\tchr2:101-200\t\t\t1\tLittle\t2025-01-04\n"
          "ISCA-5\tE region\t3p\tchr3:1001-2000\t2\tEmerging\t0\tNo evidence\t2025-01-05\n"
          "ISCA-6\tF region\t4p\tchr4:1001-2000\t1\tLittle\t0\tNo evidence\t2025-01-06\n")


def parsed():
    return br.parse_clingen_regions(TEXT.splitlines(True))


def test_parser_reads_version_scores_and_coordinates():
    version, regions = parsed()
    assert version == "2026-04-28" and len(regions) == 6
    a = regions[0]
    assert (a["isca_id"], a["chrom"], a["start"], a["end"]) == ("ISCA-1", "chr1", 1000, 5000)   # 1-based -> BED
    assert a["hi_score"] == "3" and a["ts_score"] == "0"
    assert regions[1]["ts_score"] == "3" and regions[1]["chrom"] == "chrX"
    assert regions[2]["chrom"] is None and regions[2]["hg38_location"] == "tbd"
    assert regions[3]["hi_score"] == "" and regions[3]["ts_score"] == "1"


def test_parser_needs_header():
    try:
        br.parse_clingen_regions(["#x\n", "ISCA-1\tA\t1p\tchr1:1-2\n"])
    except ValueError as e:
        assert "header" in str(e)
    else:
        raise AssertionError("no error without the #ISCA ID header")


def test_lift_regions_statuses():
    _, regions = parsed()
    # end bases (1-based): ISCA-1 1001 / 5000 -> inverted block on chr1 (ends swap);
    # ISCA-2 -> same contig, span 2.5x (dropped); ISCA-4 one end unmapped;
    # ISCA-5 ends on two contigs; ISCA-6 span +20% (kept, flagged)
    m = {("chr1", 1001): ("chr1", 9000), ("chr1", 5000): ("chr1", 5001),
         ("chrX", 20001): ("chrX", 1), ("chrX", 30000): ("chrX", 25000),
         ("chr2", 101): ("chr2", 50),
         ("chr3", 1001): ("chr3", 1), ("chr3", 2000): ("chr5", 1000),
         ("chr4", 1001): ("chr4", 1001), ("chr4", 2000): ("chr4", 2200)}
    pieces = {"ISCA-1": [("chr1", 5000, 7000), ("chr1", 7100, 9000)], "ISCA-6": [("chr4", 1000, 2200)],
              "ISCA-5": [("chr3", 5000, 5970)],                 # 97% of the region lifts in one piece: rescued
              "ISCA-4": [("chr2", 40, 80)]}                     # 40% only: stays unlifted
    out, stats = br.lift_regions(regions, lambda pts: {p: m[p] for p in pts if p in m}, lambda ivs: pieces)
    by = {r["isca_id"]: r for r in out}
    a = by["ISCA-1"]
    assert a["lift_status"] == "lifted" and (a["chrom"], a["start"], a["end"]) == ("chr1", 5000, 9000)
    assert a["interval_pieces"] == 2 and a["interval_mapped_fraction"] == "0.975"
    assert by["ISCA-2"]["lift_status"] == "unlifted_span_changed" and by["ISCA-2"]["chrom"] is None
    assert by["ISCA-3"]["lift_status"] == "unlifted_no_coordinates"
    assert by["ISCA-4"]["lift_status"] == "unlifted_end_unmapped"
    e = by["ISCA-5"]
    assert e["lift_status"] == "lifted_interval" and (e["chrom"], e["start"], e["end"]) == ("chr3", 5000, 5970)
    f = by["ISCA-6"]
    assert f["lift_status"] == "lifted_span_changed" and (f["start"], f["end"]) == (1000, 2200)
    assert stats["lifted"] == 1 and stats["lifted_span_changed"] == 1 and stats["lifted_interval"] == 1 and sum(
        v for k, v in stats.items() if not k.startswith("lifted_split")) == 6
    assert stats["lifted_split_in_hs1 (QC: >1 chain piece)"] == 1


def test_region_table_round_trip(tmp_path):
    _, regions = parsed()
    ident = {(r["chrom"], p): (r["chrom"], p) for r in regions if r["chrom"] for p in (r["start"] + 1, r["end"])}
    out, _ = br.lift_regions(regions, lambda pts: {p: ident[p] for p in pts})
    tab, rep = tmp_path / "reg.tsv", tmp_path / "reg.liftover.tsv"
    br.write_region_tables(out, str(tab), str(rep), {"chain": "test.chain", "chain_md5": "0"})
    lines = tab.read_text().splitlines()
    assert lines[0].split("\t") == br.REGION_COLUMNS and len(lines) == 1 + 5      # tbd region not written
    assert rep.read_text().startswith("## chain: test.chain") and "unlifted_no_coordinates" in rep.read_text()
    # the scorer's loader reads the table
    import score_cxsv
    regs = score_cxsv.load_regions(str(tab))
    assert {r["isca_id"] for r in regs} == {"ISCA-1", "ISCA-2", "ISCA-4", "ISCA-5", "ISCA-6"}
    assert next(r for r in regs if r["isca_id"] == "ISCA-1")["hi"] == "3"
