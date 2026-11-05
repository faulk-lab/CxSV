"""Downstream-only rebuild (scripts/panel_build/link_upstream.py): which files of
an old build are reused as hard links, and the post-build check that the old
tree did not change."""
import os
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts" / "panel_build"))
import link_upstream as lu  # noqa: E402


def write(p, text="x\n"):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(text)
    return p


@pytest.fixture
def trees(tmp_path):
    old, new = tmp_path / "old" / "results" / "panel", tmp_path / "new" / "results" / "panel"
    s = old / "samples" / "C" / "S1"
    for rel in ("merged_sv.vcf", "callers/sniffles.snf", "callers/svim/variants.vcf",
                "callability/callable.bed.gz", "cxsv/clusters.bed", "logs/svim.log", "logs/vcf_to_bed.log",
                "benchmarks/run_svim.tsv"):
        write(s / rel)
    write(old / "panel" / "sniffles_population" / "sniffles_population.vcf")
    write(old / "panel" / "plots" / "a.pdf")
    write(old / "release" / "release.json")
    dry = tmp_path / "dryrun.txt"
    dry.write_text(
        "rule sample_cluster_sv:\n"
        f"    input: {s}/cxsv/merged_sv_sorted.bed\n"
        f"    output: {s}/cxsv/clusters.bed\n"
        f"    log: {s}/logs/cluster_sv.log\n"
        "    jobid: 1\n"
        "localrule plot_panel_sharing:\n"
        f"    output: {old}/panel/plots\n"
        f"    benchmark: {s}/benchmarks/x.tsv\n")
    return old, new, dry, tmp_path


def run(old, new, dry, inv, *include_logs):
    a = lu.argparse.Namespace(old=str(old), new=str(new), dryrun=str(dry), inventory=str(inv), apply=True,
                              include_log=list(include_logs), verify=None)
    lu.build(a)


def test_planned_paths(trees):
    old, _, dry, _ = trees
    planned, rules = lu.planned_paths(dry)
    assert rules == {"sample_cluster_sv": 1, "plot_panel_sharing": 1}
    s = old / "samples" / "C" / "S1"
    assert {str(s / "cxsv/clusters.bed"), str(s / "logs/cluster_sv.log"), str(old / "panel/plots"),
            str(s / "benchmarks/x.tsv")} == planned


def test_links_upstream_only(trees):
    old, new, dry, tmp = trees
    inv = tmp / "new" / "provenance" / "reused.tsv"
    run(old, new, dry, inv, "svim.log")
    s = new / "samples" / "C" / "S1"
    for rel in ("merged_sv.vcf", "callers/sniffles.snf", "callability/callable.bed.gz", "logs/svim.log"):
        assert os.stat(s / rel).st_ino == os.stat(old / "samples" / "C" / "S1" / rel).st_ino
    assert not (s / "cxsv" / "clusters.bed").exists()          # planned output: re-made, never linked
    assert not (s / "logs" / "vcf_to_bed.log").exists()        # logs are not linked unless named
    assert not (s / "benchmarks").exists()
    assert not (new / "panel" / "plots").exists()               # under a planned directory output
    assert not (new / "release").exists()                       # the old release is never reused
    assert (new / "panel" / "sniffles_population" / "sniffles_population.vcf").exists()
    rows = [l for l in inv.read_text().splitlines() if not l.startswith("#")]
    assert len(rows) == 1 + 6
    assert lu.verify(str(inv))


def test_verify_detects_change(trees):
    old, new, dry, tmp = trees
    inv = tmp / "inv.tsv"
    run(old, new, dry, inv)
    assert lu.verify(str(inv))
    p = old / "samples" / "C" / "S1" / "merged_sv.vcf"
    p.write_text("changed in place\n")
    assert not lu.verify(str(inv))


def test_refuses_foreign_file(trees):
    old, new, dry, tmp = trees
    write(new / "samples" / "C" / "S1" / "merged_sv.vcf", "other\n")
    with pytest.raises(SystemExit):
        run(old, new, dry, tmp / "inv.tsv")
