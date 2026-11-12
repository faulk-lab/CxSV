"""Bundle validation check 7 (tests/lab/compare_patient_runs.py): the clinician report of a
bundle run and of a repository run differ by design in the run-identity row only ("Pipeline
bundle": bundle id / "not run from an offline bundle"; "Report date ...; code <commit12>").
The report text is compared with those values masked; any other text change is still flagged."""
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tests" / "lab"))
import compare_patient_runs as cpr  # noqa: E402

docx = pytest.importorskip("docx")
pytest.importorskip("reportlab")


def make_report(base, bundle, code, sex="female"):
    from reportlab.lib.pagesizes import letter
    from reportlab.platypus import Paragraph, SimpleDocTemplate, Table
    from reportlab.lib.styles import getSampleStyleSheet
    rows = [["Pipeline bundle", bundle], ["Report date", "2026-11-12" + (f"; code {code}" if code else "")],
            ["Sex-chromosome check", f"{sex} (read-depth sex call); sex_call {sex}"]]
    d = docx.Document()
    d.add_paragraph("CxSV research report")
    t = d.add_table(rows=len(rows), cols=2)
    for i, (k, v) in enumerate(rows):
        t.cell(i, 0).text, t.cell(i, 1).text = k, v
    d.save(f"{base}.clinician_report.docx")
    SimpleDocTemplate(f"{base}.clinician_report.pdf", pagesize=letter, invariant=1).build(
        [Paragraph("CxSV research report", getSampleStyleSheet()["Normal"]), Table(rows)])


def tree(root, bundle, code, sex="female"):
    rep = root / "samples" / "run" / "S" / "report"
    rep.mkdir(parents=True)
    make_report(rep / "S", bundle, code, sex)
    (root / "samples" / "run" / "S" / "x.tsv").write_text("a\tb\n")
    return root


def run(a, b):
    return subprocess.run([sys.executable, str(REPO / "tests/lab/compare_patient_runs.py"), str(a), str(b)],
                          capture_output=True, text=True)


def test_identity_row_is_an_expected_difference(tmp_path):
    a = tree(tmp_path / "bundle", "cxsv-bundle-20261112-1c9c80d0e2c6", "1c9c80d0e2c6")
    b = tree(tmp_path / "repo", "not run from an offline bundle", None)
    r = run(a, b)
    assert r.returncode == 0, r.stdout
    assert r.stdout.count("differs as expected (run identity row); report text identical") == 2
    assert "0 unexpected differences" in r.stdout


def test_other_text_change_is_flagged(tmp_path):
    a = tree(tmp_path / "bundle", "cxsv-bundle-20261112-1c9c80d0e2c6", "1c9c80d0e2c6", sex="female")
    b = tree(tmp_path / "repo", "not run from an offline bundle", None, sex="male")
    r = run(a, b)
    assert r.returncode == 1
    assert r.stdout.count("\tDIFFERS\t") == 2


def test_text_extraction_not_empty(tmp_path):
    a = tree(tmp_path / "t", "cxsv-bundle-20261112-1c9c80d0e2c6", "1c9c80d0e2c6")
    p = a / "samples" / "run" / "S" / "report" / "S.clinician_report"
    for ext in ("docx", "pdf"):
        t = cpr.report_text(f"{p}.{ext}")
        assert t and "<RUN-IDENTITY>" in t and "sex_call female" in t and "code " not in t
