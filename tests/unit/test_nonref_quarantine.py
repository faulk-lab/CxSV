"""Non-reference contig quarantine in the scorer (ROADMAP item 18; PI decision 2026-11-13) and the
locus-id check of bundle validation (check 10). BOC_5883 (panel-v1.3 test run) had a BND_cluster
locus chrN:62452837-62454007 built from two Sniffles2 BNDs (chr5:1 and chr17:1) whose ALT mate
contig is 'chrN' while CHR2=chrY."""
import json
import subprocess
import sys

import cxsv_junctions as cj
import cxsv_scoring as sc
from conftest import REPO

sys.path.insert(0, str(REPO / "tests" / "lab"))
import check_locus_ids as cli  # noqa: E402

REF = {"chr1", "chr5", "chr17", "chrY"}


def _row(chrom, s, e, cls="BND_cluster"):
    return {"Chr": chrom, "Start": str(s), "End": str(e), "Cmplx_Class": cls, "Panel_Status": "public"}


def _j(c1, p1, c2, p2, rid, loci):
    return {"chrom1": c1, "pos1": p1, "chrom2": c2, "pos2": p2, "record_ids": rid, "loci": loci}


def _boc5883_like():
    q = "chrN:62452837-62454007"
    recs = {"S.BND.1": {"chrom": "chr5", "pos": 1, "alt": "N[chrN:62452838[", "svtype": "TRA"},
            "S.BND.2": {"chrom": "chr17", "pos": 1, "alt": "N[chrN:62454006[", "svtype": "TRA"},
            "S.DEL.3": {"chrom": "chr1", "pos": 1000, "alt": "<DEL>", "svtype": "DEL"},
            "S.BND.4": {"chrom": "chr17", "pos": 5000, "alt": "N]chrY:2000]", "svtype": "BND"}}
    status = [_row("chrN", 62452837, 62454007), _row("chr1", 900, 2100, "Complex_other"),
              _row("chr17", 0, 5037000, "CN_linked_loss"), _row("chr17", 4000, 6000, "Nested_SV")]
    by_locus = {q: [_j("chr5", 1, "chrN", 62452838, "S.BND.1", [q]), _j("chr17", 1, "chrN", 62454006, "S.BND.2", [q])],
                "chr1:900-2100": [_j("chr1", 1000, "chr1", 2000, "S.DEL.3", ["chr1:900-2100"])],
                "chr17:0-5037000": [_j("chr17", 1, "chrN", 62454006, "S.BND.2", ["chr17:0-5037000"])],
                "chr17:4000-6000": [_j("chr17", 5000, "chrY", 2000, "S.BND.4", ["chr17:4000-6000"]),
                                    _j("chr17", 4500, "chrN", 7, "S.BND.2", ["chr17:4000-6000"])]}
    return status, by_locus, recs


def test_bnd_mate_contig():
    assert sc.bnd_mate_contig("N[chrN:62452838[") == "chrN"
    assert sc.bnd_mate_contig("]chr2:123]N") == "chr2"
    assert sc.bnd_mate_contig("<DEL>") is None and sc.bnd_mate_contig(None) is None


def test_chrN_bnd_cluster_quarantined_normal_loci_untouched():
    status, by_locus, recs = _boc5883_like()
    kept, quar, dropped = sc.quarantine_nonreference(status, by_locus, recs, REF)
    ids = [q["Locus_ID"] for q in quar]
    assert ids == ["chrN:62452837-62454007", "chr17:0-5037000"]
    assert quar[0]["reason"].startswith("locus contig chrN") and quar[0]["nonref_junctions"] == 2
    assert quar[1]["reason"].startswith("every junction")           # its only junction is the artefact
    assert [f"{r['Chr']}:{r['Start']}-{r['End']}" for r in kept] == ["chr1:900-2100", "chr17:4000-6000"]
    assert kept[0] is status[1]                                      # normal locus: the same row object
    assert list(dropped) == ["chr17:4000-6000"] and dropped["chr17:4000-6000"] == [{"chrN"}]
    line = sc.quarantine_line(quar, dropped)
    assert line.startswith("2 locus/loci removed from scoring and tiering; 1 junction(s) dropped from 1")
    assert "chrN" in line


def test_nothing_quarantined_on_reference_contigs():
    status, by_locus, recs = _boc5883_like()
    kept, quar, dropped = sc.quarantine_nonreference(status[1:2], by_locus, recs, REF)
    assert quar == [] and dropped == {} and kept == status[1:2]
    assert sc.quarantine_line([], {}).startswith("0 locus/loci removed")


def test_alt_mate_contig_alone_triggers(tmp_path):
    """chrom2 says chrY (as SURVIVOR's CHR2 can) but the record ALT names chrN."""
    recs = {"R": {"chrom": "chr5", "pos": 10, "alt": "N[chrN:5[", "svtype": "TRA"}}
    assert sc.junction_nonreference({"chrom1": "chr5", "chrom2": "chrY", "record_ids": "R"}, recs, REF) == {"chrN"}


def test_locus_id_check(tmp_path):
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("chr1\t1000\t6\t60\t61\nchr17\t5037000\t6\t60\t61\n")
    t = tmp_path / "s.tsv"
    t.write_text("Rank\tLocus_ID\n1\tchr1:0-2\n2\tchr17:0-5037000\n")
    assert cli.main(["--fai", str(fai), str(t)]) == 0
    for bad, why in (("chrN:62452837-62454007", "not in the reference"), ("chr1:500-500", "start >= end"),
                     ("chr1:10-1001", "past the end"), ("chr1-10-20", "not CONTIG")):
        t.write_text(f"Rank\tLocus_ID\n1\t{bad}\n")
        assert cli.main(["--fai", str(fai), str(t)]) == 1
        assert why in cli.problem(bad, {"chr1": 1000, "chr17": 5037000})


def test_validate_bundle_runs_the_check():
    v = (REPO / "tests" / "lab" / "validate_bundle.sh").read_text()
    assert "check_locus_ids.py" in v and "quarantined_loci.tsv" in v


def test_rules_wire_the_quarantine_table():
    ps = (REPO / "workflow" / "rules" / "patient_scoring.smk").read_text()
    cr = (REPO / "workflow" / "rules" / "clinician_report.smk").read_text()
    assert "--out-quarantine {output.quarantine}" in ps and "--quarantine {input.quarantine}" in cr


def test_scorer_files_are_not_fingerprinted():
    sys.path.insert(0, str(REPO / "workflow" / "lib"))
    import cxsv_fingerprint as fpm
    for f in ("scripts/score_cxsv.py", "lib/cxsv_scoring.py", "scripts/clinician_report.py",
              "rules/patient_scoring.smk", "rules/clinician_report.smk"):
        assert f not in fpm.CALLING_CODE
