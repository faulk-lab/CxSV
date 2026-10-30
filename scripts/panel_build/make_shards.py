#!/usr/bin/env python3
"""Plan the full panel build: assign every panel sample to one fleet node.

    make_shards.py --site config/sites/faulk_lab_panel.yaml \
        --fleet config/sites/faulk_lab_fleet.yaml \
        --registry /mnt/sargasso1/CxSV/datasets/registry/sample_registry.tsv \
        --cost-model <cost_model.tsv> [--build DIR] [--write] [--force]
        [--exclude NODE ...] [--set NODE.KEY=VALUE ...] [--json OUT]
        [--cohort-order C1,C2,...] [--cohorts C1,C2,...]

Without --write it only prints the estimate (what-if mode: try --exclude and
--set, e.g. --set freya.reserve_threads=16). With --write it creates

  <build>/plan/shard_plan.tsv      one row per sample: node, order, estimates
  <build>/plan/node_summary.tsv    per node: samples, CPU-h, bytes, hours
  <build>/plan/plan.json           inputs, fleet, cost model, totals
  <build>/nodes/<node>/sheets/<cohort>.tsv   the node's rows of each sheet
  <build>/nodes/<node>/node.yaml   snakemake overlay (sheets + shared results)
  <build>/nodes/<node>/runner.json settings for run_node.py
  <build>/aggregate.yaml           overlay for the aggregate step (all sheets)

Every sample of the site's cohort sheets is assigned to exactly one node (the
plan is checked to be a partition). Balancing is by estimated cost, not by
sample count: a sample's CPU-hours, critical-path wall hours and peak RSS
come from the measured cost model (per platform group: a coverage term scaled
by the sheet's coverage_x plus an alignment-size term), its NAS bytes from the
alignment's size. A node's finish time is the largest of CPU-bound (CPU-h /
usable threads), slot-bound (sum of sample wall-hours / concurrent samples)
and I/O-bound (bytes / the node's NAS read rate, at most an equal share of
the NAS's measured aggregate rate, since all nodes stage at once) time; samples are placed longest first on the node whose finish
time grows least (LPT). Within a node samples run longest first, or, with
--cohort-order, cohort by cohort in that order (a build staged by cohort: the
first cohort finishes on every node before the next starts). --cohorts limits
the estimate to some cohorts (what-if only; a written plan always covers every
cohort of the site). An existing plan is never rewritten without --force.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import shardlib  # noqa: E402

PLAN_COLUMNS = ["cohort", "sample_id", "individual_id", "platform", "group", "coverage_x",
                "alignment", "bytes", "est_cpu_h", "est_wall_h", "est_peak_rss_gb",
                "est_kept_bytes", "node", "order"]


# ─── inputs ──────────────────────────────────────────────────────────────────
def platform_group(cohort: str, chemistry: str, platform: str = "") -> str:
    """Cost-model group of a sample: what its per-x cost was measured on.
    Cohorts other than the lab's three fall back to their platform label."""
    if cohort == "1kGP_ONT_1019_CHM13":
        return "ONT_1019"
    if cohort == "1kGP_PacBio_CHM13":
        return "PacBio_HiFi"
    if cohort == "1kGP_LRSC_ONT_CHM13":
        return f"LRSC_{chemistry or 'R10'}"
    return platform or "default"


def load_samples(site: dict, registry_rows: list[dict], stat_sizes: bool = False) -> list[dict]:
    reg = {r["file_id"]: r for r in registry_rows}
    out = []
    for cohort, spec in sorted(site["panel_cohorts"].items()):
        _, rows = shardlib.read_sheet(spec["samples"])
        for r in rows:
            fid = f"{cohort}/{os.path.basename(r['alignment'])}"
            g = reg.get(fid, {})
            size = int(g.get("size_bytes") or 0)
            if stat_sizes or not size:
                size = os.path.getsize(r["alignment"])
            out.append({"cohort": cohort, "sample_id": r["sample_id"],
                        "individual_id": r.get("individual_id") or r["sample_id"],
                        "platform": r.get("platform") or spec.get("platform", ""),
                        "group": platform_group(cohort, g.get("chemistry", ""),
                                                r.get("platform") or spec.get("platform", "")),
                        "coverage_x": float(r.get("coverage_x") or 0), "alignment": r["alignment"],
                        "bytes": size})
    return out


def read_cost_model(path: str) -> dict:
    _, rows = shardlib.read_sheet(path)
    return {r["group"]: {k: (float(v) if k not in ("group", "test_samples") else v)
                         for k, v in r.items()} for r in rows}


def estimate_sample(s: dict, cm: dict) -> dict:
    m = cm[s["group"]]
    cov = s["coverage_x"]
    ratio = cov / m["test_coverage_x"] if m["test_coverage_x"] else 1.0
    s = dict(s)
    gb = s.get("bytes", 0) / 1e9
    s["est_cpu_h"] = m["cpu_h_per_x"] * cov + m.get("cpu_h_per_gb", 0.0) * gb
    s["est_wall_h"] = m["wall_h_per_x"] * cov + m.get("wall_h_per_gb", 0.0) * gb
    # memory grows with depth in dense regions; never assume less than measured
    s["est_peak_rss_gb"] = m["peak_rss_gb"] * max(1.0, ratio)
    s["est_kept_bytes"] = m["kept_bytes_per_x"] * cov
    return s


# ─── nodes ───────────────────────────────────────────────────────────────────
def node_capacity(name: str, n: dict, cm: dict, ref_speed: float, nas_share_mb_s: float = 1e9) -> dict:
    usable = max(1, int(n["threads"]) - int(n.get("reserve_threads", 0)))
    busy = max(m["busy_threads"] for m in cm.values())
    rss = max(m["peak_rss_gb"] for m in cm.values())
    mem = float(n["ram_gb"]) - float(n.get("reserve_ram_gb", 0))
    if str(n.get("slots", "auto")) == "auto":
        slots = max(1, min(math.floor(usable / busy), math.floor(mem / rss)))
    else:
        slots = int(n["slots"])
    return {"node": name, "usable_threads": usable, "slots": slots,
            "factor": float(n.get("speed", ref_speed)) / ref_speed,
            # all nodes stage at once: each gets about an equal share of the NAS
            "io_mb_s": min(float(n.get("nic_mb_s", 100)), nas_share_mb_s), "mem_gb": mem,
            "stage_max_gb": float(n.get("stage_max_gb", 1e9))}


def finish_h(load: dict, cap: dict) -> dict:
    cpu = load["cpu_h"] * cap["factor"] / cap["usable_threads"]
    slot = load["wall_h"] * cap["factor"] / cap["slots"]
    io = load["bytes"] / (cap["io_mb_s"] * 1e6) / 3600
    return {"cpu_bound_h": cpu, "slot_bound_h": slot, "io_bound_h": io,
            "finish_h": max(cpu, slot, io, load["max_wall_h"] * cap["factor"])}


def assign(samples: list[dict], caps: dict) -> dict:
    """LPT: longest (CPU-h) first onto the node whose finish time grows least."""
    loads = {n: {"cpu_h": 0.0, "wall_h": 0.0, "bytes": 0, "max_wall_h": 0.0, "n": 0} for n in caps}
    where = {}
    for s in sorted(samples, key=lambda x: (-x["est_cpu_h"], x["cohort"], x["sample_id"])):
        best = None
        for n, cap in caps.items():
            if s["bytes"] / 1e9 > cap["stage_max_gb"] or s["est_peak_rss_gb"] > cap["mem_gb"]:
                continue
            ld = dict(loads[n])
            ld["cpu_h"] += s["est_cpu_h"]
            ld["wall_h"] += s["est_wall_h"]
            ld["bytes"] += s["bytes"]
            ld["max_wall_h"] = max(ld["max_wall_h"], s["est_wall_h"])
            f = finish_h(ld, cap)["finish_h"]
            if best is None or f < best[0] - 1e-9:
                best = (f, n, ld)
        if best is None:
            raise ValueError(f"no enabled node can take {s['cohort']}/{s['sample_id']} "
                             f"({s['bytes'] / 1e9:.0f} GB, {s['est_peak_rss_gb']:.1f} GB RSS)")
        _, n, ld = best
        ld["n"] += 1
        loads[n] = ld
        where[(s["cohort"], s["sample_id"])] = n
    return where


def check_partition(samples: list[dict], where: dict) -> None:
    keys = [(s["cohort"], s["sample_id"]) for s in samples]
    if len(set(keys)) != len(keys):
        raise ValueError("duplicate (cohort, sample_id) in the cohort sheets")
    missing = [k for k in keys if k not in where]
    extra = [k for k in where if k not in set(keys)]
    if missing or extra:
        raise ValueError(f"plan is not a partition: {len(missing)} unassigned, {len(extra)} unknown")


def summarise(samples, where, caps, cm, nas_mb_s) -> tuple[list[dict], dict]:
    rows = []
    for n, cap in caps.items():
        mine = [s for s in samples if where[(s["cohort"], s["sample_id"])] == n]
        ld = {"cpu_h": sum(s["est_cpu_h"] for s in mine), "wall_h": sum(s["est_wall_h"] for s in mine),
              "bytes": sum(s["bytes"] for s in mine),
              "max_wall_h": max([s["est_wall_h"] for s in mine] or [0])}
        f = finish_h(ld, cap)
        groups = defaultdict(int)
        for s in mine:
            groups[s["group"]] += 1
        rows.append({"node": n, "samples": len(mine), **{f"n_{g}": groups.get(g, 0) for g in sorted(cm)},
                     "slots": cap["slots"], "usable_threads": cap["usable_threads"],
                     "speed_factor": round(cap["factor"], 3),
                     "cpu_h_node": round(ld["cpu_h"] * cap["factor"], 1),
                     "nas_tb": round(ld["bytes"] / 1e12, 2),
                     "kept_gb": round(sum(s["est_kept_bytes"] for s in mine) / 1e9, 1),
                     **{k: round(v, 1) for k, v in f.items()}})
    tot_bytes = sum(s["bytes"] for s in samples)
    totals = {"samples": len(samples), "individuals": len({s["individual_id"] for s in samples}),
              "cpu_h_ref": round(sum(s["est_cpu_h"] for s in samples), 1),
              "cpu_h_fleet": round(sum(r["cpu_h_node"] for r in rows), 1),
              "nas_tb": round(tot_bytes / 1e12, 2),
              "nas_floor_h": round(tot_bytes / (nas_mb_s * 1e6) / 3600, 1),
              "kept_tb": round(sum(s["est_kept_bytes"] for s in samples) / 1e12, 2),
              "makespan_h": round(max([r["finish_h"] for r in rows] + [tot_bytes / (nas_mb_s * 1e6) / 3600]), 1)}
    return rows, totals


# ─── writing ─────────────────────────────────────────────────────────────────
def write_plan(build, site, fleet, cm_path, samples, where, caps, rows, totals, force=False,
               cohort_order=()):
    plan_dir = f"{build}/plan"
    if os.path.exists(f"{plan_dir}/shard_plan.tsv") and not force:
        sys.exit(f"{plan_dir}/shard_plan.tsv exists; a plan is never rewritten without --force "
                 f"(re-sharding a running build: docs/PANEL_BUILD.md)")
    os.makedirs(plan_dir, exist_ok=True)
    results = f"{build}/results/panel"
    order = defaultdict(int)
    by_node = defaultdict(list)
    rank = {c: i for i, c in enumerate(cohort_order)}
    for s in sorted(samples, key=lambda x: (rank.get(x["cohort"], len(rank)), -x["est_wall_h"],
                                            x["cohort"], x["sample_id"])):
        n = where[(s["cohort"], s["sample_id"])]
        order[n] += 1
        s["node"], s["order"] = n, order[n]
        by_node[n].append(s)
    plan_rows = [{k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in s.items()}
                 for n in caps for s in by_node[n]]
    shardlib.write_sheet(f"{plan_dir}/shard_plan.tsv", PLAN_COLUMNS, plan_rows)
    shardlib.write_sheet(f"{plan_dir}/node_summary.tsv", list(rows[0].keys()), rows)
    fl_nodes = fleet["nodes"]
    for n, mine in by_node.items():
        nd = f"{build}/nodes/{n}"
        os.makedirs(f"{nd}/sheets", exist_ok=True)
        cohorts = {}
        for cohort, spec in site["panel_cohorts"].items():
            header, rows_ = shardlib.read_sheet(spec["samples"])
            ids = {s["sample_id"] for s in mine if s["cohort"] == cohort}
            keep = [r for r in rows_ if r["sample_id"] in ids]
            if keep:
                shardlib.write_sheet(f"{nd}/sheets/{cohort}.tsv", header, keep)
                cohorts[cohort] = {"samples": f"{nd}/sheets/{cohort}.tsv",
                                   "platform": spec.get("platform", "unknown")}
        nf = fl_nodes[n]
        shardlib.dump_json(f"{nd}/node.yaml",
                           shardlib.node_overlay(site, cohorts, results, int(nf["sv_threads"])))
        shardlib.dump_json(f"{nd}/runner.json", {
            "node": n, "build": build, "code": f"{build}/code", "sif": fleet["sif"],
            "results_root": results, "slots": caps[n]["slots"],
            "cores_per_sample": int(nf["cores_per_sample"]), "sv_threads": int(nf["sv_threads"]),
            "threads": int(nf["threads"]), "reserve_threads": int(nf.get("reserve_threads", 0)),
            "reserve_ram_gb": float(nf.get("reserve_ram_gb", 0)),
            "stage_dir": nf["stage_dir"], "tmp_dir": nf["tmp_dir"],
            "stage_max_gb": float(nf.get("stage_max_gb", 1000)),
            "max_concurrent_stage": int(fleet.get("max_concurrent_stage", 1)),
            "stage_bwlimit_mb_s": int(nf.get("stage_bwlimit_mb_s", 0) or 0),
            "max_attempts": int(fleet.get("max_attempts", 2)),
            "active_hours": nf.get("active_hours") or "",
            "peak_rss_gb": {g: m["peak_rss_gb"] for g, m in read_cost_model(cm_path).items()},
            "site_cohort_platform": {c: sp.get("platform", "") for c, sp in site["panel_cohorts"].items()},
        })
    agg = shardlib.node_overlay(site, {c: {"samples": sp["samples"], "platform": sp.get("platform", "unknown")}
                                       for c, sp in site["panel_cohorts"].items()}, results, 16)
    if fleet.get("panel_individuals") and not agg.get("panel_individuals"):
        agg["panel_individuals"] = fleet["panel_individuals"]
    shardlib.dump_json(f"{build}/aggregate.yaml", agg)
    shardlib.dump_json(f"{plan_dir}/plan.json", {"fleet": fleet, "cost_model": read_cost_model(cm_path),
                                                 "capacity": caps, "totals": totals,
                                                 "site_cohorts": site["panel_cohorts"]})


def apply_sets(fleet: dict, sets: list[str]) -> None:
    for s in sets or []:
        k, _, v = s.partition("=")
        node, _, key = k.partition(".")
        fleet["nodes"][node][key] = shardlib._scalar(v)


def plan(site, fleet, registry_rows, cm, exclude=(), stat_sizes=False, cohorts=None):
    ref_speed = float(next(iter(cm.values()))["ref_speed"])
    active = [n for n, nd in fleet["nodes"].items() if nd.get("enabled", True) and n not in exclude]
    share = float(fleet.get("nas_aggregate_mb_s", 600)) / max(1, len(active))
    caps = {n: node_capacity(n, fleet["nodes"][n], cm, ref_speed, share) for n in active}
    samples = [estimate_sample(s, cm) for s in load_samples(site, registry_rows, stat_sizes)
               if not cohorts or s["cohort"] in cohorts]
    where = assign(samples, caps)
    check_partition(samples, where)
    rows, totals = summarise(samples, where, caps, cm, float(fleet.get("nas_aggregate_mb_s", 600)))
    return samples, where, caps, rows, totals


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--site", required=True)
    ap.add_argument("--fleet", required=True)
    ap.add_argument("--registry", default=None, help="registry TSV (chemistry, sizes); optional")
    ap.add_argument("--cost-model", required=True)
    ap.add_argument("--build", default=None, help="default: fleet build_root")
    ap.add_argument("--exclude", nargs="*", default=[])
    ap.add_argument("--set", nargs="*", default=[], dest="sets")
    ap.add_argument("--write", action="store_true")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--json", default=None, help="write the estimate (rows + totals) here")
    ap.add_argument("--cohort-order", default="", help="comma list: run cohorts in this order on every node")
    ap.add_argument("--cohorts", default="", help="comma list: estimate these cohorts only (no --write)")
    a = ap.parse_args()
    site = shardlib.load_yaml(a.site)
    fleet = shardlib.load_yaml(a.fleet)
    apply_sets(fleet, a.sets)
    reg = shardlib.read_sheet(a.registry)[1] if a.registry else []
    cm = read_cost_model(a.cost_model)
    only = [c for c in a.cohorts.split(",") if c]
    if only and a.write:
        sys.exit("--cohorts is for estimates only; a written plan covers every cohort (use --cohort-order)")
    samples, where, caps, rows, totals = plan(site, fleet, [r for r in reg if r.get("status") == "verified"],
                                              cm, a.exclude, cohorts=only)
    cols = list(rows[0].keys())
    print("\t".join(cols))
    for r in rows:
        print("\t".join(str(r[c]) for c in cols))
    print(json.dumps(totals))
    if a.json:
        shardlib.dump_json(a.json, {"rows": rows, "totals": totals, "exclude": a.exclude, "set": a.sets})
    if a.write:
        build = a.build or fleet["build_root"]
        write_plan(build, site, fleet, a.cost_model, samples, where, caps, rows, totals, a.force,
                   [c for c in a.cohort_order.split(",") if c])
        print(f"plan written under {build}/")


if __name__ == "__main__":
    main()
