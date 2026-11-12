"""Patient intake helpers (workflow/lib/cxsv_intake.py, roadmap item 16k)."""
import io
import struct

import pytest

import cxsv_intake as ci
from conftest import REPO

FLAGSTAT = """1210 + 0 in total (QC-passed reads + QC-failed reads)
1000 + 0 primary
0 + 0 secondary
210 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
1198 + 0 mapped (99.01% : N/A)
988 + 0 primary mapped (98.80% : N/A)
"""
GOOD_PG = ("@PG\tID:minimap2\tPN:minimap2\tVN:2.24-r1122\tCL:minimap2 -t 26 -y -ax map-ont --eqx -k 17 -K 4g "
           "-R @RG\\tID:P1\\tSM:P1 /refs/hs1.map-ont.k17.mmi -")
HEADER = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:10\n@RG\tID:P1\tSM:P1\tPL:ONT\n" + GOOD_PG + "\n"


def sheet(tmp_path, rows, header="sample_id\tinput\tindividual_id\tplatform\tinput_reference\tdescription"):
    p = tmp_path / "sheet.tsv"
    p.write_text(header + "\n" + "\n".join(rows) + "\n")
    return p


def touch(p):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("x")
    return str(p)


def test_input_format():
    assert [ci.input_format(p) for p in ("a.mod.bam", "a.CRAM", "r.fastq.gz", "r.fq", "r.FASTQ")] == \
        ["bam", "cram", "fastq", "fastq", "fastq"]
    with pytest.raises(ci.IntakeSheetError):
        ci.input_format("reads.pod5")


def test_sheet_multi_input_defaults_and_formats(tmp_path):
    a, b = touch(tmp_path / "fc1.bam"), touch(tmp_path / "fc2.bam")
    fq = touch(tmp_path / "r.fq.gz")
    recs = ci.read_intake_sheet(sheet(tmp_path, [f"P1\t{a},{b}\tFAM1\t\t\tdorado sup", f"P2\t{fq}\t\tont\t\t"]))
    assert recs[0].inputs == (a, b) and recs[0].individual_id == "FAM1" and recs[0].platform == "ONT"
    assert recs[0].format_summary == "bam" and recs[0].description == "dorado sup"
    assert recs[1].individual_id == "P2" and recs[1].platform == "ONT" and recs[1].formats == ("fastq",)


def test_sheet_minimal_columns(tmp_path):
    a = touch(tmp_path / "x.cram")
    recs = ci.read_intake_sheet(sheet(tmp_path, [f"P1\t{a}"], header="sample_id\tinput"))
    assert recs[0].sample_id == "P1" and recs[0].input_reference == ""


def test_sheet_lists_every_problem(tmp_path):
    a = touch(tmp_path / "x.bam")
    rows = [f"P1\t{a}\t\tONT\t\t", f"P1\t{a}\t\tONT\t\t",           # duplicate
            f"bad id\t{a}\t\tONT\t\t",                            # bad characters
            f"P3\t{a}\t\tPACBIO\t\t",                             # map-ont recipe only
            f"P4\t{tmp_path}/missing.bam\t\tONT\t\t",             # missing input
            f"P5\t{a},{a}\t\tONT\t\t",                            # same input twice
            f"P6\t{a}\t\tONT\t{tmp_path}/nope.fa\t",              # missing input_reference
            f"P7\treads.pod5\t\tONT\t\t"]                         # unsupported type
    with pytest.raises(ci.IntakeSheetError) as e:
        ci.read_intake_sheet(sheet(tmp_path, rows))
    msg = str(e.value)
    for frag in ("duplicate sample_id 'P1'", "sample_id 'bad id'", "platform 'PACBIO' not supported",
                 "input not found", "same input is listed twice", "input_reference not found", "unsupported input type"):
        assert frag in msg, frag
    assert "7 problem(s)" in msg      # the first P1 row is valid


def test_sheet_requires_columns(tmp_path):
    with pytest.raises(ci.IntakeSheetError, match="lacks required column 'input'"):
        ci.read_intake_sheet(sheet(tmp_path, ["P1\tx"], header="sample_id\talignment"))


def test_keep_tags_filters_non_sam_tokens_and_other_tags():
    c = "runid=abc qs:f:12.5 MM:Z:C+m?,0,3; ML:B:C,200,13 RG:Z:run1_dorado ch=7 MN:i:100 xx:Z:keepnot"
    assert ci.keep_tags(c) == "qs:f:12.5\tMM:Z:C+m?,0,3;\tML:B:C,200,13\tMN:i:100"
    assert ci.keep_tags("") == "" and ci.keep_tags("MM=Z:bad") == ""


def test_normalise_fastq(tmp_path):
    fq = tmp_path / "r.fq"
    fq.write_text("@r1 runid=1 MM:Z:C+m?,0; ML:B:C,9\nACGT\n+\nIIII\n@r2\tqs:f:9.5\tch=3\nAC\n+\nII\n@r3\nA\n+\nI\n")
    out = io.StringIO()
    assert ci.normalise_fastq([str(fq)], out) == 3
    assert out.getvalue().splitlines()[::4] == ["@r1\tMM:Z:C+m?,0;\tML:B:C,9", "@r2\tqs:f:9.5", "@r3"]


@pytest.mark.parametrize("text,err", [("@r1\nACGT\n+\nIII\n", "lengths differ"), ("@r1\nACGT\n", "truncated")])
def test_normalise_fastq_rejects_bad_records(tmp_path, text, err):
    fq = tmp_path / "r.fq"
    fq.write_text(text)
    with pytest.raises(ci.IntakeCheckError, match=err):
        ci.normalise_fastq([str(fq)], io.StringIO())


def write_mmi(path, w=10, k=17, n=25, flag=0):
    path.write_bytes(b"MMI\x02" + struct.pack("<5I", w, k, 14, n, flag) + b"\0" * 16)
    return path


def test_mmi_check(tmp_path):
    assert ci.check_mmi(write_mmi(tmp_path / "ok.mmi"), 25)["k"] == 17
    with pytest.raises(ci.IntakeCheckError, match="recipe needs k=17"):
        ci.check_mmi(write_mmi(tmp_path / "k15.mmi", k=15))
    with pytest.raises(ci.IntakeCheckError, match="sequences"):
        ci.check_mmi(write_mmi(tmp_path / "n.mmi", n=24), 25)
    with pytest.raises(ci.IntakeCheckError, match="homopolymer"):
        ci.check_mmi(write_mmi(tmp_path / "hpc.mmi", flag=1))
    (tmp_path / "x.mmi").write_bytes(b"BAM\x01" + b"\0" * 30)
    with pytest.raises(ci.IntakeCheckError, match="not a minimap2 index"):
        ci.check_mmi(tmp_path / "x.mmi")


def test_parse_flagstat_and_fastq_count():
    fs = ci.parse_flagstat(FLAGSTAT)
    assert (fs["primary"], fs["primary_mapped"], fs["unmapped"], fs["supplementary"]) == (1000, 988, 12, 210)
    assert ci.parse_samtools_fastq_count("[M::bam2fq_mainloop] discarded 0 singletons\n"
                                         "[M::bam2fq_mainloop] processed 1000 reads\n") == 1000


def test_verify_pass_and_read_count_failure():
    res = ci.verify("P1", 1000, FLAGSTAT, HEADER, mod_records_in=2000, mod_records_out=1990)
    assert res["verdict"] == "PASS" and res["out_unmapped"] == 12
    with pytest.raises(ci.IntakeCheckError, match="input primary 1001 != output primary mapped 988"):
        ci.verify("P1", 1001, FLAGSTAT, HEADER)
    with pytest.raises(ci.IntakeCheckError, match="output sample has none"):
        ci.verify("P1", 1000, FLAGSTAT, HEADER, mod_records_in=50, mod_records_out=0)


@pytest.mark.parametrize("mutate,err", [
    (lambda h: h.replace("SM:P1\tPL", "SM:pos_control_3559\tPL"), "@RG SM"),       # the BOC_3559 case
    (lambda h: h + "@RG\tID:x\tSM:P1\n", "2 @RG lines"),
    (lambda h: h.replace("VN:2.24-r1122", "VN:2.28-r1209"), "minimap2 2.28-r1209"),
    (lambda h: h.replace(" --eqx", ""), "lacks --eqx"),
    (lambda h: h.replace(" -y", ""), "lacks -y"),
    (lambda h: h.replace("-k 17", "-k 15"), "lacks -k 17"),
    (lambda h: h.replace("-ax map-ont", "-ax map-hifi"), "lacks -x map-ont"),
    (lambda h: h + GOOD_PG.replace("ID:minimap2", "ID:minimap2.1") + "\n", "2 minimap2 @PG"),
])
def test_check_header_catches_recipe_drift(mutate, err):
    assert ci.check_header(HEADER, "P1") == []
    probs = ci.check_header(mutate(HEADER), "P1")
    assert any(err in p for p in probs), probs


def test_check_header_accepts_split_preset():
    h = HEADER.replace("-ax map-ont", "-y -x map-ont -a").replace(" -y -y", " -y")
    assert ci.check_header(h, "P1") == []


def test_patient_sheet_rows_load_as_patient_run(tmp_path):
    import cxsv_inputs
    recs = [ci.IntakeRecord("P1", ("a.bam",), "FAM1"), ci.IntakeRecord("P2", ("b.fq",), "P2")]
    for r in recs:
        touch(tmp_path / r.sample_id / f"{r.sample_id}.hs1.cram")
        touch(tmp_path / r.sample_id / f"{r.sample_id}.hs1.cram.crai")
    rows = ci.patient_sheet_rows(recs, lambda s: tmp_path / s / f"{s}.hs1.cram")
    ci.write_tsv(rows, tmp_path / "ps.tsv", ["sample_id", "alignment", "index", "individual_id", "platform"])
    table, _ = cxsv_inputs.build_sample_table({"run": {"samples": str(tmp_path / "ps.tsv")}}, "bamfiles")
    assert table[("run", "P1")].individual_id == "FAM1" and table[("run", "P2")].is_cram


def test_intake_is_outside_the_calling_fingerprint():
    import cxsv_fingerprint as fpm
    assert not any("intake" in rel for rel in fpm.CALLING_CODE)
    assert "intake" not in fpm.CALLING_PARAM_KEYS


def test_recipe_is_the_same_in_the_script_and_the_library():
    sh = (REPO / "workflow/scripts/intake_align.sh").read_text()
    assert 'RECIPE_MM="2.24-r1122"' in sh and ci.RECIPE["minimap2_version"] == "2.24-r1122"
    assert '-y -ax map-ont --eqx -k 17 -K "$MM_BATCH"' in sh
    assert "samtools fastq -T MM,ML,MN,qs" in sh and ",".join(ci.TAGS) == "MM,ML,MN,qs"


def test_container_env_pins_minimap2_2_24():
    explicit = (REPO / "envs/cxsv-pipeline.explicit.txt").read_text()
    assert "/minimap2-2.24-" in explicit
