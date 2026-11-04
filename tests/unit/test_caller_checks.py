"""SVIM output check (workflow/scripts/check_svim_output.py, rules/caller_checks.smk)."""
import pytest

import check_svim_output as cso
import cxsv_fingerprint as fpm

CONTIGS = [f"chr{i}" for i in (1, 2, 3, 4)] + ["chrM", "chrX", "chrY"]
DONE_LOG = "2026-10-30 14:23:16,090 [INFO   ]  Draw plots..\n2026-10-30 14:23:33,056 [INFO   ]  Done.\n"


def vcf(path, counts):
    with open(path, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        for c in CONTIGS:
            fh.write(f"##contig=<ID={c},length=1000000>\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for c in CONTIGS:
            for i in range(counts.get(c, 0)):
                fh.write(f"{c}\t{100 + i}\tsv{c}_{i}\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL\n")
    return str(path)


@pytest.fixture
def sample(tmp_path):
    def make(svim, sniffles, log=DONE_LOG):
        (tmp_path / "svim.log").write_text(log)
        return (vcf(tmp_path / "svim.vcf", svim), str(tmp_path / "svim.log"), vcf(tmp_path / "sniffles.vcf", sniffles),
                str(tmp_path / "out.tsv"))
    return make


FULL_SNIFFLES = {"chr1": 50, "chr2": 40, "chr3": 30, "chr4": 25, "chrM": 2, "chrX": 30, "chrY": 3}


def test_complete_output_passes_and_low_count_contigs_are_not_judged(sample):
    s, log, ref, out = sample({"chr1": 900, "chr2": 800, "chr3": 600, "chr4": 500, "chrX": 400}, FULL_SNIFFLES)
    assert cso.main(["--svim-vcf", s, "--svim-log", log, "--reference-vcf", ref, "--out", out]) == 0
    rows = [ln.split("\t") for ln in open(out) if not ln.startswith(("#", "contig"))]
    st = {r[0]: r[3].strip() for r in rows}
    assert st["chr4"] == "ok" and st["chrY"] == "not_judged" and st["chrM"] == "not_judged"


def test_output_cut_after_a_contig_fails_and_names_the_point(sample, capsys):
    """The 2026-10-31 failure: missing chr7 cache file -> no SVIM calls from there on."""
    s, log, ref, out = sample({"chr1": 900, "chr2": 800}, FULL_SNIFFLES)
    assert cso.main(["--svim-vcf", s, "--svim-log", log, "--reference-vcf", ref, "--out", out]) == 1
    err = capsys.readouterr().err
    assert "output ends at chr2" in err and "chr3,chr4,chrX" in err
    assert not (pytest.importorskip("pathlib").Path(out).exists())


def test_single_missing_contig_fails(sample):
    s, log, ref, out = sample({"chr1": 900, "chr3": 600, "chr4": 500, "chrX": 400}, FULL_SNIFFLES)
    probs, summary, _ = cso.check(s, log, ref, 20)
    assert summary["contigs_missing"] == "chr2" and summary["status"] == "FAIL"
    assert "output ends at" not in probs[0]          # calls resume after chr2: not a cut-off


@pytest.mark.parametrize("log,expect", [
    ("[E::cram_decode_slice] Failure to decode slice\n" + DONE_LOG, "decode/truncation/reference error"),
    ("OSError: truncated file\n" + DONE_LOG, "decode/truncation/reference error"),
    ("2026-10-30 [INFO   ]  Draw plots..\n", "does not end with SVIM's 'Done.'"),
])
def test_log_errors_fail(sample, log, expect):
    s, lg, ref, out = sample({"chr1": 900, "chr2": 800, "chr3": 600, "chr4": 500, "chrX": 400}, FULL_SNIFFLES, log)
    probs, _, _ = cso.check(s, lg, ref, 20)
    assert any(expect in p for p in probs)


def test_sliced_alignment_passes(sample):
    """Regression slices: both callers see the same few contigs; the rest are not judged."""
    s, log, ref, out = sample({"chr3": 40}, {"chr3": 5, "chr4": 1})
    assert cso.check(s, log, ref, 20)[0] == []


def test_both_workflows_and_the_panel_build_run_the_check():
    from conftest import REPO
    for wf in ("panel", "patient"):
        assert 'include: "../rules/caller_checks.smk"' in (REPO / f"workflow/{wf}/Snakefile").read_text()
    assert '"qc/svim_check.tsv"' in (REPO / "workflow/panel_build/Snakefile").read_text()
    assert 'svim_ok   = f"{S}/qc/svim_check.tsv"' in (REPO / "workflow/rules/patient_scoring.smk").read_text()


def test_check_is_outside_the_calling_fingerprint():
    """Adding the check must not change the calling fingerprint (no panel re-calling)."""
    assert "rules/caller_checks.smk" not in fpm.CALLING_CODE
    assert "scripts/check_svim_output.py" not in fpm.CALLING_CODE
    assert "svim_check" not in fpm.CALLING_PARAM_KEYS
