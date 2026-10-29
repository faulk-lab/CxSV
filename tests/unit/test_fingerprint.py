"""Calling fingerprint (workflow/lib/cxsv_fingerprint.py)."""
import copy
import shutil

import cxsv_fingerprint as fpm
from conftest import REPO

PARAMS = {"sv_params": {"min_sv_len": 50, "min_support": 3}, "filter_params": {"min_qual": 10},
          "survivor_params": {"max_dist": 300}, "cxsv_params": {"cluster_distance": 50000},
          "junction_params": {"collapse_distance": 300}}


def make_fp(tmp_path, params=PARAMS, tools=None, contigs=None):
    m = tmp_path / "md5_manifest.tsv"
    m.write_text("SN\tLN\tM5\n" + "".join(f"{n}\t{l}\t{h}\n" for n, l, h in
                 (contigs or [("chr1", 100, "a" * 32), ("chr2", 50, "b" * 32)])))
    t = tmp_path / "tools.tsv"
    t.write_text("".join(f"{k}\t{v}\n" for k, v in (tools or {"sniffles": "2.8.1", "cuteSV": "2.1.4"}).items()))
    return fpm.build(dict(params, threads={"sv_callers": 8}), REPO / "workflow", m, t)


def test_identical_inputs_are_compatible(tmp_path):
    a, b = make_fp(tmp_path), make_fp(tmp_path)
    assert a["sha256"] == b["sha256"] and fpm.compare(a, b) == []


def test_threads_are_not_fingerprinted(tmp_path):
    a = make_fp(tmp_path)
    assert "threads" not in a["params"]


def test_parameter_tool_and_reference_differences_are_reported(tmp_path):
    a = make_fp(tmp_path)
    p = copy.deepcopy(PARAMS)
    p["sv_params"]["min_support"] = 2
    b = make_fp(tmp_path, params=p, tools={"sniffles": "2.8.2", "cuteSV": "2.1.4"},
                contigs=[("chr1", 100, "c" * 32), ("chr2", 50, "b" * 32)])
    diffs = fpm.compare(b, a)
    assert any(d.startswith("params.sv_params.min_support") for d in diffs)
    assert any(d.startswith("tools.sniffles") for d in diffs)
    assert any(d.startswith("reference.chr1") for d in diffs)
    assert a["sha256"] != b["sha256"]


def test_code_hash_ignores_comments_but_not_code(tmp_path):
    wd = tmp_path / "workflow"
    shutil.copytree(REPO / "workflow", wd)
    base = fpm.code_digest(wd)["sha256"]
    target = wd / "scripts" / "classify_cxsv.py"
    target.write_text(target.read_text() + "\n# a comment\n\n")
    assert fpm.code_digest(wd)["sha256"] == base
    target.write_text(target.read_text() + "X = 1\n")
    assert fpm.code_digest(wd)["sha256"] != base


def test_every_fingerprinted_code_file_exists():
    for rel in fpm.CALLING_CODE:
        assert (REPO / "workflow" / rel).is_file(), rel
