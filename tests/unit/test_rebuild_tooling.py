"""Downstream-only rebuild tooling for builds with the copy-number track
(panel-v1.2, roadmap item 15): make_aggregate_yaml.py copy-number keys,
check_aggregate_config.py, rebuild_downstream.sh --plan-old and the early key
check, tree_snapshot.py."""
import json
import os
import stat
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
PB = REPO / "scripts" / "panel_build"
sys.path.insert(0, str(PB))
import check_aggregate_config as cac  # noqa: E402
import make_aggregate_yaml as may  # noqa: E402
import tree_snapshot as ts  # noqa: E402


@pytest.fixture
def builds(tmp_path):
    old, new = tmp_path / "old", tmp_path / "new"
    old.mkdir(); new.mkdir()
    (old / "aggregate.yaml").write_text(json.dumps({
        "panel_results_root": str(old / "results" / "panel"),
        "release": {"name": "cxsv-panel", "date": "20261105", "description": "old"}}))
    (tmp_path / "desc.txt").write_text("panel-v1.2\n  test   build\n")
    depth = tmp_path / "cn_depth_v1"; depth.mkdir()
    strata = tmp_path / "strata.tsv"; strata.write_text("cohort\tsample_id\tcn_stratum\n")
    return old, new, tmp_path / "desc.txt", depth, strata


def test_make_aggregate_yaml_writes_cn_keys(builds):
    old, new, desc, depth, strata = builds
    may.main([str(old), str(new), str(desc), "--cn-depth-precomputed", str(depth), "--cn-strata-table", str(strata)])
    cfg = json.load(open(new / "aggregate.yaml"))
    assert cfg["panel_results_root"] == str(new.resolve() / "results" / "panel")
    assert cfg["release"]["description"] == "panel-v1.2 test build" and cfg["release"]["date"] is None
    assert cfg["cn_depth_precomputed"] == str(depth) and cfg["cn_strata_table"] == str(strata)
    assert cac.problems(str(new / "aggregate.yaml")) == []


def test_make_aggregate_yaml_without_cn_keys_and_missing_path(builds):
    old, new, desc, depth, strata = builds
    may.main([str(old), str(new), str(desc)])
    cfg = json.load(open(new / "aggregate.yaml"))
    assert "cn_depth_precomputed" not in cfg
    with pytest.raises(SystemExit):
        may.main([str(old), str(new), str(desc), "--cn-depth-precomputed", str(depth / "nope")])


def test_check_aggregate_config(tmp_path, builds):
    _, _, _, depth, strata = builds
    p = tmp_path / "a.yaml"
    p.write_text(json.dumps({"cn_depth_precomputed": "none"}))
    msg = cac.problems(str(p))
    assert any("cn_depth_precomputed unset" in m for m in msg) and any("cn_strata_table unset" in m for m in msg)
    p.write_text(json.dumps({"cn_depth_precomputed": str(depth / "x"), "cn_strata_table": str(strata)}))
    assert [m for m in cac.problems(str(p))] == [f"{p}: cn_depth_precomputed = {depth / 'x'} does not exist"]
    p.write_text(json.dumps({"cn_depth_precomputed": str(strata), "cn_strata_table": str(strata)}))
    assert len(cac.problems(str(p))) == 1          # a file is not a depth directory
    p.write_text("not json: [")
    assert "not readable as JSON" in cac.problems(str(p))[0]
    assert cac.main([str(p)]) == 1


def _new_code(new):
    pb = new / "code" / "scripts" / "panel_build"
    pb.mkdir(parents=True)
    (pb / "lib.sh").write_text((PB / "lib.sh").read_text())
    (new / "code" / "workflow").mkdir()
    (new / "code" / "COMMIT").write_text("abc\n")
    (new / "code" / "workflow" / "COMMIT").write_text("abc\n")


def _fake_snakemake(tmp_path):
    b = tmp_path / "bin"; b.mkdir()
    f = b / "snakemake"
    f.write_text("""#!/usr/bin/env bash
cfg=""; for ((i=1;i<=$#;i++)); do [ "${!i}" = --configfile ] && { j=$((i+1)); cfg="${!j}"; }; done
case " $* " in
  *" --list-rules "*) printf 'sample_vcf_to_bed\\ncn_depth\\nrun_svim\\ncalling_fingerprint\\n' ;;
  *" -n "*) echo "config $cfg"; echo "results_root $(python3 -c "import json,sys; print(json.load(open(sys.argv[1]))['panel_results_root'])" "$cfg")"
            printf 'rule sample_vcf_to_bed:\\n    output: /x/samples/C/S/cxsv/merged.bed\\n\\nJob stats:\\njob count\\nsample_vcf_to_bed 1\\ntotal 1\\n\\n' ;;
esac
""")
    f.chmod(f.stat().st_mode | stat.S_IEXEC)
    return b


def _run(args, tmp_path, binpath=None):
    env = dict(os.environ, CXSV_SIF="none", HOME=str(tmp_path), CXSV_REBUILD_WORKDIR=str(tmp_path / "wd"))
    if binpath:
        env["PATH"] = f"{binpath}:{env['PATH']}"
    return subprocess.run(["bash", str(PB / "rebuild_downstream.sh"), *args], capture_output=True, text=True, env=env)


def test_rebuild_downstream_stops_without_cn_keys(builds, tmp_path):
    old, new, desc, depth, strata = builds
    may.main([str(old), str(new), str(desc)])
    r = _run([str(old), str(new), "2", "--dry-run"], tmp_path)
    assert r.returncode == 1 and "copy-number keys missing" in r.stdout and "cn_strata_table unset" in r.stderr


def test_rebuild_downstream_plan_old(builds, tmp_path):
    old, new, desc, depth, strata = builds
    may.main([str(old), str(new), str(desc), "--cn-depth-precomputed", str(depth), "--cn-strata-table", str(strata)])
    _new_code(new)
    r = _run([str(old), str(new), "2", "--plan-old"], tmp_path, _fake_snakemake(tmp_path))
    assert r.returncode == 0, r.stdout + r.stderr
    out = (new / "provenance" / "dryrun_oldtree.txt").read_text()
    assert f"results_root {old / 'results' / 'panel'}" in out           # the OLD tree is planned against
    cfg = json.load(open(tmp_path / "wd" / "plan_old" / "aggregate_oldtree.yaml"))
    assert cfg["cn_depth_precomputed"] == str(depth)
    assert "old-tree plan:" in r.stdout


def test_tree_snapshot_compare(tmp_path):
    root = tmp_path / "b"; (root / "results" / "panel" / "release").mkdir(parents=True)
    (root / "results" / "panel" / "release" / "r.json").write_text("{}")
    (root / "x.txt").write_text("1")
    a, b = tmp_path / "a.tsv", tmp_path / "b.tsv"
    assert ts.main([str(root), "--out", str(a), "--sha256", "results/panel/release"]) == 0
    rows = ts.read(str(a))
    assert rows["results/panel/release/r.json"].split("\t")[5] != "" and rows["x.txt"].split("\t")[5] == ""
    assert ts.main([str(root), "--out", str(b), "--sha256", "results/panel/release"]) == 0
    assert ts.main(["--compare", str(a), str(b)]) == 0
    os.utime(root / "x.txt", ns=(1, 1))
    (root / "new.txt").write_text("")
    assert ts.main([str(root), "--out", str(b), "--sha256", "results/panel/release"]) == 0
    d = ts.compare(str(a), str(b))
    assert any(x.startswith("changed: x.txt") for x in d) and "only after: new.txt" in d
    assert ts.main(["--compare", str(a), str(b)]) == 1
