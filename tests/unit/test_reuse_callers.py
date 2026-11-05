"""tests/lab/lib/reuse_callers.sh and its use in tests/lab/pos_control_rerun_downstream.sh.

The behaviour (snakemake --touch on a toy workflow) is tested in CI by
tests/synthetic/check_rerun_touch.sh; these checks need only bash.
"""
import re
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
LIB = REPO / "tests/lab/lib/reuse_callers.sh"
RERUN = REPO / "tests/lab/pos_control_rerun_downstream.sh"


def _bash(script: str) -> subprocess.CompletedProcess:
    return subprocess.run(["bash", "-c", script], capture_output=True, text=True)


def _upstream_rules() -> list[str]:
    r = _bash(f'source "{LIB}"; printf "%s\\n" "${{CXSV_UPSTREAM_RULES[@]}}"')
    assert r.returncode == 0, r.stderr
    return r.stdout.split()


def test_scripts_parse():
    for p in (LIB, RERUN, REPO / "tests/synthetic/check_rerun_touch.sh"):
        r = subprocess.run(["bash", "-n", str(p)], capture_output=True, text=True)
        assert r.returncode == 0, f"{p}: {r.stderr}"


def test_upstream_rules_exist_in_workflow():
    defined = set()
    for smk in (REPO / "workflow/rules").glob("*.smk"):
        defined |= set(re.findall(r"^(?:local)?rule (\w+):", smk.read_text(), re.M))
    rules = _upstream_rules()
    assert rules and set(rules) <= defined, set(rules) - defined


def test_upstream_rules_are_exactly_the_guarded_ones_plus_the_reference_check():
    """touch_upstream may only cover what the step-2 guard forbids re-running, plus check_cram_reference."""
    guard = re.search(r"\(run_sniffles\|run_cutesv\|run_svim\|filter_\|merge_sv_callers\)", RERUN.read_text())
    assert guard, "step-2 caller guard not found"
    rules = set(_upstream_rules())
    assert "check_cram_reference" in rules
    for r in rules - {"check_cram_reference"}:
        assert re.match(r"run_(sniffles|cutesv|svim)$|filter_|merge_sv_callers$", r), r
    assert not any(r.startswith("sample_") for r in rules)


def test_same_reference_manifest(tmp_path):
    a, b, c = tmp_path / "a.tsv", tmp_path / "b.tsv", tmp_path / "c.tsv"
    a.write_text("name\tlength\tmd5\nchr1\t1\tx\n"); b.write_text(a.read_text()); c.write_text("other\n")
    ok = lambda x, y: _bash(f'source "{LIB}"; same_reference_manifest "{x}" "{y}"').returncode == 0
    assert ok(a, b)
    assert not ok(a, c)
    assert not ok(a, tmp_path / "missing.tsv")


def test_rerun_script_touches_only_with_identical_manifest_and_without_forcerun():
    s = RERUN.read_text()
    i_check, i_touch = s.index("same_reference_manifest"), s.index("touch_upstream \"/work/$S_REL/merged_sv.vcf\" \"${X[@]}\" \"${SM_BASE[@]}\"")
    assert i_check < i_touch < s.index('hdr "2. dry run"')
    base = re.search(r"SM_BASE=\((.*?)\)\n", s, re.S).group(1)
    assert "--forcerun" not in base and "--rerun-triggers mtime" in base


def _fake_snakemake(tmp_path, rules):
    """A stand-in for snakemake: prints a plan with RULES and logs each call's arguments."""
    log = tmp_path / "calls.log"
    fake = tmp_path / "fake_snakemake.sh"
    plan = "".join(f"rule {r}:\\n" for r in rules)
    fake.write_text(f'#!/usr/bin/env bash\necho "$*" >> "{log}"\nprintf "{plan}"\n')
    fake.chmod(0o755)
    return fake, log


def test_touch_upstream_touches_after_listing_upstream_only_plan(tmp_path):
    fake, log = _fake_snakemake(tmp_path, ["check_cram_reference", "run_svim", "merge_sv_callers"])
    r = _bash(f'source "{LIB}"; touch_upstream T/merged_sv.vcf "{fake}" --cores 1')
    assert r.returncode == 0, r.stdout + r.stderr
    calls = log.read_text().splitlines()
    assert calls == ["--cores 1 --touch -n T/merged_sv.vcf", "--cores 1 --touch T/merged_sv.vcf"]


def test_touch_upstream_refuses_plan_with_downstream_rule(tmp_path):
    fake, log = _fake_snakemake(tmp_path, ["run_svim", "sample_vcf_to_bed"])
    r = _bash(f'source "{LIB}"; touch_upstream T/merged_sv.vcf "{fake}"')
    assert r.returncode == 1
    assert "sample_vcf_to_bed outside CXSV_UPSTREAM_RULES" in r.stdout
    assert log.read_text().splitlines() == ["--touch -n T/merged_sv.vcf"]   # no touching call
