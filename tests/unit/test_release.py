"""Panel release identity, manifest verification (tamper detection) and callability (roadmap item 2)."""
import gzip
import subprocess
import sys

import pytest

import cxsv_callability as cc
import cxsv_release as cr
from conftest import REPO
from test_panel_scripts import classify, panel_release, read

SCRIPTS = REPO / "workflow" / "scripts"
FP = "3f2a9c1d0b4e" + "0" * 52


# ─── release id and source ──────────────────────────────────────────────────
def test_release_id_format_and_validation():
    assert cr.release_id("cxsv-panel", "20261029", FP) == "cxsv-panel-20261029-3f2a9c1d0b4e"
    for bad in [("bad name", "20261029", FP), ("x", "2026-10-29", FP), ("x", "20261029", "abc")]:
        with pytest.raises(ValueError):
            cr.release_id(*bad)


def test_source_from_git_archive_commit_file(tmp_path):
    (tmp_path / "workflow").mkdir()
    (tmp_path / "workflow" / "COMMIT").write_text("a" * 40 + "\n")
    info = cr.source_info(tmp_path)
    assert info["commit"] == "a" * 40 and info["dirty"] is False
    (tmp_path / "workflow" / "COMMIT").write_text("$Format:%H$\n")      # not exported by git archive
    assert cr.source_info(tmp_path)["commit"] == "unknown"


# ─── manifest: verify and refuse ────────────────────────────────────────────
def make_release(tmp_path):
    rel = tmp_path / "release"
    (rel / "callability").mkdir(parents=True)
    (rel / "release.json").write_text('{"release_id": "x-20261029-3f2a9c1d0b4e"}\n')
    (rel / "panel_junction_sites.tsv").write_text("chrom1\tpos1\nchr1\t100\n")
    (rel / "callability" / "samples.tsv").write_text("mask\n")
    cr.write_manifest(rel)
    return rel


def test_fresh_release_verifies(tmp_path):
    rel = make_release(tmp_path)
    res = cr.verify(rel)
    assert res["ok"] and res["n_files"] == 3, res
    assert [r["file"] for r in cr.read_manifest(rel)] == ["callability/samples.tsv", "panel_junction_sites.tsv",
                                                         "release.json"]


@pytest.mark.parametrize("tamper, expect", [
    (lambda r: (r / "panel_junction_sites.tsv").write_text("chrom1\tpos1\nchr1\t101\n"), "sha256 differs"),
    (lambda r: (r / "panel_junction_sites.tsv").write_text("chrom1\tpos1\nchr1\t1000\n"), "size differs"),
    (lambda r: (r / "callability" / "samples.tsv").unlink(), "missing: callability/samples.tsv"),
    (lambda r: (r / "extra.tsv").write_text("x\n"), "not in manifest: extra.tsv"),
])
def test_tampered_release_is_refused(tmp_path, tamper, expect):
    rel = make_release(tmp_path)
    tamper(rel)
    res = cr.verify(rel)
    assert not res["ok"] and any(p.startswith(expect) for p in res["problems"]), res["problems"]
    run = subprocess.run([sys.executable, SCRIPTS / "verify_panel_release.py", rel, "--out", tmp_path / "ok.txt"],
                         capture_output=True, text=True)
    assert run.returncode == 1 and "status: FAILED" in run.stdout and not (tmp_path / "ok.txt").exists()


def test_release_without_manifest_is_refused(tmp_path):
    (tmp_path / "r").mkdir()
    assert not cr.verify(tmp_path / "r")["ok"]


# ─── callability rules ──────────────────────────────────────────────────────
def test_breakpoint_callable_needs_one_full_flank():
    m = cc.Mask({"chr1": [(1000, 2000), (2000, 2050), (5000, 6000)]})     # touching intervals merge
    # 1100 both sides inside; 1000 right side only; 2050 left side only (1000-2050 after merging);
    # 2051 and 2100 neither; 4950 straddles the gap; 5000 right side; 6050 past the end
    pos = [1100, 1000, 2050, 2051, 2100, 4950, 5000, 6050]
    got = m.breakpoints_callable("chr1", pos, 100).tolist()
    assert got == [True, True, True, False, False, False, True, False]
    assert not m.breakpoints_callable("chr2", [1500], 100).any()                    # contig without mask
    assert m.callable_bp() == 2050


def test_junction_needs_both_ends():
    m = cc.Mask({"chr1": [(0, 10_000)], "chr2": [(0, 500)]})
    js = [{"chrom1": "chr1", "pos1": 500, "chrom2": "chr1", "pos2": 9000},
          {"chrom1": "chr1", "pos1": 500, "chrom2": "chr1", "pos2": 20_000},
          {"chrom1": "chr1", "pos1": 500, "chrom2": "chr2", "pos2": 300}]
    assert cc.junctions_callable(m, js, 100).tolist() == [True, False, True]


def write_masks(d, masks):
    """masks: {(cohort, sample, individual, platform): [(chrom, s, e)]}."""
    d.mkdir(parents=True, exist_ok=True)
    with open(d / cc.SAMPLES_TSV, "w") as fh:
        fh.write("\t".join(cc.SAMPLES_COLUMNS) + "\n")
        for (cohort, sample, ind, plat), iv in masks.items():
            name = cc.mask_name(cohort, sample)
            with gzip.open(d / name, "wt") as g:
                for c, s, e in iv:
                    g.write(f"{c}\t{s}\t{e}\n")
            fh.write(f"{name}\t{cohort}\t{sample}\t{ind}\t{plat}\t{sum(e - s for _, s, e in iv)}\n")
    return d


def test_patient_loci_report_callability_next_to_carriers(tmp_path):
    common = [(10_000, 20_000, "INV"), (12_000, 25_000, "DEL"), (21_000, 30_000, "INV")]
    out = panel_release(tmp_path, [
        ("A", "A", "ONT", [("chr1", 10_000, 30_000)], [(*j, "chr1:10000-30000") for j in common]),
        ("B", "B", "ONT", [], []),
        ("C", "C", "PACBIO", [], []),
    ])
    whole = [("chr1", 0, 1_000_000)]
    masks = write_masks(tmp_path / "callability", {
        ("c", "A", "A", "ONT"): whole,
        ("c", "B", "B", "ONT"): [("chr1", 0, 22_000)],             # misses the 25-30 kb breakpoints
        ("c", "C", "C", "PACBIO"): whole,
    })
    loci = [("chr1", 10_000, 30_000), ("chr1", 60_000, 70_000)]
    jrows = [(*j, "chr1:10000-30000") for j in common] + [(60_000, 70_000, "INV", "chr1:60000-70000")]
    rows, res = classify(tmp_path, out, loci, jrows, callability_dir=masks, flank=100)
    got = [(r["Panel_Status"], r["Panel_Carriers"], r["Panel_Callable"], r["Panel_Individuals"]) for r in rows]
    assert got == [("public", "1", "2", "3"), ("private", "0", "2", "3")]
    js = read(res["junctions"])
    assert [j["panel_callable"] for j in js] == ["3", "2", "2", "2"]
    # site callability on the panel side: the three common sites
    subprocess.run([sys.executable, SCRIPTS / "panel_site_callability.py", "--sites", out["sites"],
                    "--callability-dir", masks, "--panel-individuals", out["inds"], "--flank", "100",
                    "--out", tmp_path / "sc.tsv", "--summary", tmp_path / "scs.tsv"], check=True)
    sc = read(tmp_path / "sc.tsv")
    assert [(r["n_callable"], r["callable_ONT"], r["callable_PACBIO"], r["n_individuals"]) for r in sc] == \
        [("3", "2", "1", "3"), ("2", "1", "1", "3"), ("2", "1", "1", "3")]
    summ = {(r["level"], r["name"]): r["value"] for r in read(tmp_path / "scs.tsv")}
    assert summ[("sites", "callable_in_all_individuals")] == "1"
