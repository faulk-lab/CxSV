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


def run(old, new, dry, inv, *include_logs, exclude=()):
    a = lu.argparse.Namespace(old=str(old), new=str(new), dryrun=str(dry), inventory=str(inv), apply=True,
                              include_log=list(include_logs), verify=None, exclude=list(exclude))
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


@pytest.fixture
def cn_trees(tmp_path):
    """panel-v1.3 case: OLD already has the CN track; the new code re-runs cn_depth and
    cn_linkage of S1 (planned) but not of S2 (held), and cn_params_json is params-only."""
    old, new = tmp_path / "old" / "results" / "panel", tmp_path / "new" / "results" / "panel"
    write(old / "cn_params.json", "{}\n")
    write(old / "tool_versions.tsv")
    for smp in ("S1", "S2"):
        s = old / "samples" / "C" / smp
        for rel in ("merged_sv.vcf", "cn/depth.regions.bed.gz", "cn/depth.regions.bed.gz.csi",
                    "cn/depth.mosdepth.summary.txt", "cn/segments.tsv", "cn/stats.json",
                    "cn/work/linkage_candidates.tsv", "cn/work/sub/x.tmp"):
            write(s / rel)
    s1 = old / "samples" / "C" / "S1"
    dry = tmp_path / "dryrun.txt"
    dry.write_text(
        "rule cn_depth:\n"
        f"    output: {s1}/cn/depth.regions.bed.gz\n"
        "    jobid: 1\n"
        "rule sample_cn_linkage:\n"
        f"    input: {old}/cn_params.json, {s1}/cn/segments.tsv\n"
        f"    output: {s1}/cn/cn_loci.tsv\n"
        f"    log: {s1}/logs/cn_linkage.log\n"
        "    jobid: 2\n")
    return old, new, dry, tmp_path


def test_params_only_and_cn_side_files_not_linked(cn_trees, capsys):
    old, new, dry, tmp = cn_trees
    inv = tmp / "new" / "provenance" / "reused.tsv"
    run(old, new, dry, inv)
    out = capsys.readouterr().out
    s1, s2 = new / "samples" / "C" / "S1", new / "samples" / "C" / "S2"
    assert not (new / "cn_params.json").exists()               # params-only: re-made from the new config
    assert (new / "tool_versions.tsv").exists()
    assert (s1 / "merged_sv.vcf").exists()                      # upstream of S1: still reused
    for rel in ("cn/depth.regions.bed.gz", "cn/depth.regions.bed.gz.csi", "cn/depth.mosdepth.summary.txt",
                "cn/segments.tsv", "cn/stats.json", "cn/work/linkage_candidates.tsv", "cn/work/sub/x.tmp"):
        assert not (s1 / rel).exists(), rel                     # planned cn dir: nothing in it is linked
        assert (s2 / rel).exists(), rel                         # held cn dir: reused as before
    assert "'params-only': 1" in out and "'cn side files': 6" in out and "'planned': 1" in out
    rows = [l for l in inv.read_text().splitlines() if not l.startswith("#")]
    assert len(rows) == 1 + 1 + 2 + 7                           # header, tool_versions, 2 merged, S2 cn
    assert lu.verify(str(inv))


def test_exclude_regex(cn_trees):
    old, new, dry, tmp = cn_trees
    run(old, new, dry, tmp / "inv.tsv", exclude=[r"^tool_versions\.tsv$"])
    assert not (new / "tool_versions.tsv").exists()
    assert (new / "samples" / "C" / "S2" / "cn" / "stats.json").exists()


def test_side_file_reason():
    dirs = {"samples/C/S1/cn"}
    assert lu.side_file_reason("cn_params.json", dirs) == "params-only"
    assert lu.side_file_reason("samples/C/S1/cn/work/a.tsv", dirs) == "cn side files"
    assert lu.side_file_reason("samples/C/S1/cn/depth.regions.bed.gz.csi", dirs) == "cn side files"
    assert lu.side_file_reason("samples/C/S2/cn/work/a.tsv", dirs) is None
    assert lu.side_file_reason("samples/C/S1/callers/svim/variants.vcf", dirs) is None
    assert lu.side_file_reason("panel/x/cn_params.json", dirs) is None


def test_linked_side_files_check(cn_trees):
    """The rebuild driver's post-link check: a hard-linked side file in NEW stops the run."""
    old, new, dry, tmp = cn_trees
    s1o, s1n = old / "samples" / "C" / "S1", new / "samples" / "C" / "S1"
    (s1n / "cn" / "work").mkdir(parents=True)
    newdry = tmp / "dry_new.txt"
    newdry.write_text(dry.read_text().replace(str(old), str(new)))
    assert lu.linked_side_files(str(newdry), str(new)) == []
    os.link(s1o / "cn" / "work" / "linkage_candidates.tsv", s1n / "cn" / "work" / "linkage_candidates.tsv")
    os.link(old / "cn_params.json", new / "cn_params.json")
    assert lu.linked_side_files(str(newdry), str(new)) == ["cn_params.json",
                                                           "samples/C/S1/cn/work/linkage_candidates.tsv"]
