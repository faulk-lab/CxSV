"""Sample discovery (workflow/lib/cxsv_inputs.py)."""
import pytest

import cxsv_inputs as ci


def touch(p):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("")
    return p


def test_directory_mode_bam_and_cram(tmp_path):
    touch(tmp_path / "bam" / "c1" / "A.bam")
    touch(tmp_path / "bam" / "c1" / "B.cram")
    table, notes = ci.build_sample_table({"c1": {"platform": "ONT"}}, tmp_path / "bam")
    assert sorted(table) == [("c1", "A"), ("c1", "B")]
    a, b = table[("c1", "A")], table[("c1", "B")]
    assert a.index.endswith("A.bam.bai") and b.index.endswith("B.cram.crai")
    assert (a.individual_id, a.platform, b.is_cram) == ("A", "ONT", True)
    assert "1 CRAM, 1 BAM" in notes[0]


def test_directory_mode_rejects_same_sample_twice(tmp_path):
    touch(tmp_path / "c1" / "A.bam")
    touch(tmp_path / "c1" / "A.cram")
    with pytest.raises(ci.SampleSheetError, match="both"):
        ci.samples_from_dir("c1", tmp_path / "c1")


def write_sheet(path, rows, header=("sample_id", "alignment", "individual_id", "platform", "coverage_x")):
    path.write_text("\t".join(header) + "\n" + "".join("\t".join(r) + "\n" for r in rows))
    return path


def test_sheet_min_coverage_individual_and_platform(tmp_path):
    for n in ("s1.bam", "s1.bam.bai", "s2.cram", "s2.cram.crai"):
        touch(tmp_path / n)
    sheet = write_sheet(tmp_path / "sheet.tsv", [
        ("s1", str(tmp_path / "s1.bam"), "HG001", "ONT", "30"),
        ("s2", str(tmp_path / "s2.cram"), "HG001", "PACBIO", "12"),
    ])
    recs, dropped = ci.samples_from_sheet("c", sheet, min_coverage=17)
    assert dropped == 1 and [r.sample_id for r in recs] == ["s1"]
    assert (recs[0].individual_id, recs[0].platform, recs[0].coverage_x) == ("HG001", "ONT", "30")


def test_sheet_reports_every_problem(tmp_path):
    sheet = write_sheet(tmp_path / "sheet.tsv", [
        ("s1", str(tmp_path / "missing.bam"), "", "", "20"),
        ("s1", str(tmp_path / "x.txt"), "", "", "20"),
        ("s3", str(tmp_path / "s3.bam"), "", "", "n/a"),
    ])
    with pytest.raises(ci.SampleSheetError) as e:
        ci.samples_from_sheet("c", sheet, min_coverage=17)
    msg = str(e.value)
    assert "alignment not found" in msg and "must end in .bam or .cram" in msg and "not a number" in msg


def test_min_coverage_requires_sheet(tmp_path):
    with pytest.raises(ci.SampleSheetError, match="needs a sample sheet"):
        ci.build_sample_table({"c": {"min_coverage": 17}}, tmp_path)
