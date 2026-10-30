"""Panel build planning (scripts/panel_build/): shard partition, balancing,
config handling. The runner itself is exercised end-to-end by
tests/synthetic/run_panel_build_synthetic.sh."""
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "scripts" / "panel_build"))
import make_shards as ms  # noqa: E402
import shardlib  # noqa: E402

CM = {g: {"group": g, "test_coverage_x": 30.0, "cpu_h_per_x": 0.5, "wall_h_per_x": 0.3,
          "busy_threads": 2.0, "peak_rss_gb": 10.0, "kept_bytes_per_x": 1e8, "stage_mb_s": 500.0,
          "ref_speed": 1.0}
      for g in ("ONT_1019", "LRSC_R10", "LRSC_R9", "PacBio_HiFi")}


def node(threads=32, ram=120, speed=1.0, nic=500, **kw):
    d = dict(threads=threads, ram_gb=ram, reserve_threads=0, reserve_ram_gb=0, slots="auto",
             cores_per_sample=6, sv_threads=4, speed=speed, nic_mb_s=nic, stage_max_gb=2000)
    d.update(kw)
    return d


def samples(n=40):
    out = []
    for i in range(n):
        g = ("ONT_1019", "LRSC_R10", "LRSC_R9", "PacBio_HiFi")[i % 4]
        cov = 10 + (i * 7) % 40
        out.append(ms.estimate_sample({"cohort": f"c{i % 4}", "sample_id": f"S{i:03d}", "individual_id": f"I{i % 30}",
                                       "platform": "ONT", "group": g, "coverage_x": float(cov),
                                       "alignment": f"/x/S{i}.bam", "bytes": int(cov * 3e9)}, CM))
    return out


def test_platform_groups():
    assert ms.platform_group("1kGP_ONT_1019_CHM13", "") == "ONT_1019"
    assert ms.platform_group("1kGP_LRSC_ONT_CHM13", "R9") == "LRSC_R9"
    assert ms.platform_group("1kGP_LRSC_ONT_CHM13", "R10") == "LRSC_R10"
    assert ms.platform_group("1kGP_PacBio_CHM13", "") == "PacBio_HiFi"
    assert ms.platform_group("other", "", "ONT") == "ONT"


def test_estimates_scale_with_coverage_and_never_shrink_memory():
    s = ms.estimate_sample({"group": "LRSC_R10", "coverage_x": 60.0}, CM)
    assert s["est_cpu_h"] == pytest.approx(30.0) and s["est_wall_h"] == pytest.approx(18.0)
    assert s["est_peak_rss_gb"] == pytest.approx(20.0)
    low = ms.estimate_sample({"group": "LRSC_R10", "coverage_x": 5.0}, CM)
    assert low["est_peak_rss_gb"] == pytest.approx(10.0)


def test_slots_auto_is_bounded_by_threads_and_memory():
    assert ms.node_capacity("a", node(threads=32, ram=200), CM, 1.0)["slots"] == 16   # 32 / 2 busy
    assert ms.node_capacity("b", node(threads=32, ram=60), CM, 1.0)["slots"] == 6     # 60 / 10 GB
    assert ms.node_capacity("c", node(slots=3), CM, 1.0)["slots"] == 3


def test_assignment_is_a_partition_and_balanced_by_cost():
    ss = samples(80)
    caps = {n: ms.node_capacity(n, nd, CM, 1.0) for n, nd in
            {"fast": node(), "slow": node(speed=2.0), "thin": node(threads=8, ram=40)}.items()}
    where = ms.assign(ss, caps)
    ms.check_partition(ss, where)
    rows, totals = ms.summarise(ss, where, caps, CM, 600)
    fin = {r["node"]: r["finish_h"] for r in rows}
    assert max(fin.values()) / min(fin.values()) < 1.35, fin
    n = {r["node"]: r["samples"] for r in rows}
    assert n["fast"] > n["slow"] and n["fast"] > n["thin"], n
    assert totals["samples"] == 80 and totals["individuals"] == 30


def test_partition_check_catches_gaps():
    ss = samples(4)
    where = {(s["cohort"], s["sample_id"]): "a" for s in ss[:3]}
    with pytest.raises(ValueError, match="not a partition"):
        ms.check_partition(ss, where)


def test_io_bound_node_gets_fewer_bytes():
    ss = samples(60)
    caps = {"lan10g": ms.node_capacity("lan10g", node(nic=600), CM, 1.0),
            "lan1g": ms.node_capacity("lan1g", node(nic=110), CM, 1.0)}
    where = ms.assign(ss, caps)
    b = {k: sum(s["bytes"] for s in ss if where[(s["cohort"], s["sample_id"])] == k) for k in caps}
    assert b["lan10g"] > b["lan1g"]


def test_node_too_small_for_a_sample_is_skipped():
    ss = [ms.estimate_sample({"cohort": "c", "sample_id": "BIG", "individual_id": "I", "platform": "ONT",
                              "group": "PacBio_HiFi", "coverage_x": 30.0, "alignment": "/x", "bytes": int(400e9)}, CM)]
    caps = {"small": ms.node_capacity("small", node(stage_max_gb=300), CM, 1.0),
            "big": ms.node_capacity("big", node(), CM, 1.0)}
    assert ms.assign(ss, caps)[("c", "BIG")] == "big"


def test_site_config_has_no_coverage_floor_and_lists_three_cohorts():
    site = shardlib.load_yaml(str(REPO / "config/sites/faulk_lab_panel.yaml"))
    assert sorted(site["panel_cohorts"]) == ["1kGP_LRSC_ONT_CHM13", "1kGP_ONT_1019_CHM13", "1kGP_PacBio_CHM13"]
    assert all("min_coverage" not in spec for spec in site["panel_cohorts"].values())


def test_mini_yaml_reads_the_site_and_fleet_configs_like_pyyaml():
    yaml = pytest.importorskip("yaml")
    for f in ("config/sites/faulk_lab_panel.yaml", "config/sites/faulk_lab_fleet.yaml"):
        text = (REPO / f).read_text()
        assert shardlib._mini_yaml(text) == yaml.safe_load(text), f


def test_node_overlay_replaces_cohorts_and_never_carries_calling_params():
    site = {"reference": "r.fa", "ref_cache": "c", "annotation": {"genes_bed": "g"},
            "panel_cohorts": {"A": {"samples": "a.tsv"}, "B": {"samples": "b.tsv"}},
            "sv_params": {"min_support": 99}, "release": {"name": "x"}}
    ov = shardlib.node_overlay(site, {"A": {"samples": "shard_a.tsv"}}, "/res", 4)
    assert ov["panel_cohorts"] == {"A": {"samples": "shard_a.tsv"}}
    assert "sv_params" not in ov and ov["release"] == {"name": "x"}
    assert ov["threads"]["sv_callers"] == 4 and ov["panel_results_root"] == "/res"


def test_staged_row_points_at_local_copies():
    r = shardlib.staged_row({"alignment": "/nas/c/X.cram", "index": ""}, "/local/stage", "coh")
    assert r["alignment"] == "/local/stage/coh/X.cram" and r["index"] == "/local/stage/coh/X.cram.crai"


def test_fleet_config_nodes_and_local_disks():
    fleet = shardlib.load_yaml(str(REPO / "config/sites/faulk_lab_fleet.yaml"))
    assert set(fleet["nodes"]) == {"freya", "shadowfax", "athena", "betsy", "princess", "keurig"}
    for n, nd in fleet["nodes"].items():
        # tmpfs /tmp on betsy, princess, keurig: staging and temp files must be on disk
        assert not nd["stage_dir"].startswith("/tmp") and not nd["tmp_dir"].startswith("/tmp"), n


def test_runner_memory_gate_counts_running_samples(monkeypatch):
    """MemAvailable cannot see memory a just-started sample will still take."""
    import run_node
    r = run_node.Runner.__new__(run_node.Runner)
    r.cfg = {"reserve_ram_gb": 6, "threads": 32, "reserve_threads": 2}
    r.plan = [{"cohort": "c", "sample_id": f"s{i}", "est_peak_rss_gb": "20"} for i in range(3)]
    r.plan.append({"cohort": "c", "sample_id": "big", "est_peak_rss_gb": "50"})
    r.running = {}
    monkeypatch.setattr(run_node.Runner, "mem_available_gb", staticmethod(lambda: 58.0))
    monkeypatch.setattr(run_node.Runner, "mem_total_gb", staticmethod(lambda: 61.0))
    monkeypatch.setattr(run_node.os, "getloadavg", lambda: (1.0, 1.0, 1.0))
    assert r.can_start(r.plan[0])
    r.running = {("c", "s0"): {}, ("c", "s1"): {}}
    assert not r.can_start(r.plan[2])        # 20 + 20 + 20 > 61 - 6, though MemAvailable looks fine
    r.cfg["rss_commit_factor"] = 0.5
    assert r.can_start(r.plan[2])            # (20 + 20 + 20) x 0.5 <= 55
    r.running, r.cfg["rss_commit_factor"] = {}, 1.0
    assert r.can_start(r.plan[3])            # alone on the node: only MemAvailable decides
