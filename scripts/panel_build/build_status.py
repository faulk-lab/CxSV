#!/usr/bin/env python3
"""Progress and cost summary of a panel build, from the plan, the per-node
ledgers and the per-sample cost records. Read-only on everything it reads.

    build_status.py BUILD [--write] [--require-complete] [--costs]

Prints one line per node: planned / done / running / staging / failed /
pending samples, CPU-hours spent (from cost records), estimated CPU-hours
left, hours since the node's last ledger event, and an ETA from the node's
own completion rate. --write also writes (freya only; temp + rename)
  BUILD/status/status.tsv       the same table
  BUILD/status/samples.tsv      one row per sample: node, state, attempts, seconds
  BUILD/status/node_costs.tsv   per node and rule: CPU-h, wall-h, max RSS (--costs)
--require-complete exits 1 unless every planned sample is done (aggregate.sh
uses it before the aggregate step).
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import sys
import time
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import shardlib  # noqa: E402


def read_ledger(path):
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def sample_states(plan, ledgers, results_root):
    """state per (cohort, sample): done | running | staging | failed | pending (+ attempts, seconds)."""
    last = {}
    attempts = defaultdict(int)
    seconds = {}
    for rows in ledgers.values():
        for r in rows:
            if not r["sample_id"]:
                continue
            k = (r["cohort"], r["sample_id"])
            last[k] = r
            if r["event"] == "started":
                attempts[k] += 1
            if r["event"] == "done":
                seconds[k] = r["seconds"]
    out = {}
    for p in plan:
        k = (p["cohort"], p["sample_id"])
        done = os.path.exists(f"{results_root}/samples/{k[0]}/{k[1]}/panel_build/cost.tsv")
        ev = last.get(k, {}).get("event", "")
        if done:
            st = "done"
        elif ev in ("started", "adopted"):
            st = "running"
        elif ev in ("staging", "staged"):
            st = "staging" if ev == "staging" else "staged"
        elif ev == "gave_up":
            st = "failed"
        else:
            st = "pending"
        out[k] = {"state": st, "attempts": attempts.get(k, 0), "seconds": seconds.get(k, ""),
                  "last_event": ev, "last_time": last.get(k, {}).get("time", "")}
    return out


def read_costs(plan, results_root):
    rows = []
    for p in plan:
        f = f"{results_root}/samples/{p['cohort']}/{p['sample_id']}/panel_build/cost.tsv"
        if os.path.exists(f):
            with open(f) as fh:
                for r in csv.DictReader(fh, delimiter="\t"):
                    r["node"] = p["node"]
                    rows.append(r)
    return rows


def fnum(x, default=0.0):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("build")
    ap.add_argument("--write", action="store_true")
    ap.add_argument("--require-complete", action="store_true")
    ap.add_argument("--costs", action="store_true")
    a = ap.parse_args()
    _, plan = shardlib.read_sheet(f"{a.build}/plan/shard_plan.tsv")
    results = f"{a.build}/results/panel"
    nodes = sorted({p["node"] for p in plan})
    ledgers = {n: read_ledger(f"{a.build}/nodes/{n}/ledger.tsv") for n in nodes}
    states = sample_states(plan, ledgers, results)
    costs = read_costs(plan, results) if (a.costs or a.write) else []
    cpu_done = defaultdict(float)
    for r in costs:
        cpu_done[r["node"]] += fnum(r["cpu_time"]) / 3600
    now = time.time()
    table = []
    for n in nodes:
        mine = [p for p in plan if p["node"] == n]
        cnt = defaultdict(int)
        for p in mine:
            cnt[states[(p["cohort"], p["sample_id"])]["state"]] += 1
        left = sum(fnum(p["est_cpu_h"]) for p in mine if states[(p["cohort"], p["sample_id"])]["state"] != "done")
        evs = [r for r in ledgers[n] if r.get("epoch")]
        last = (now - fnum(evs[-1]["epoch"])) / 3600 if evs else float("nan")
        first_start = next((fnum(r["epoch"]) for r in evs if r["event"] == "start_runner"), None)
        eta = ""
        if first_start and cnt["done"]:
            rate = cnt["done"] / max(1e-9, (now - first_start) / 3600)
            eta = f"{(len(mine) - cnt['done']) / rate:.1f}"
        runner_alive = bool(evs) and evs[-1]["event"] != "exit_runner"
        table.append({"node": n, "planned": len(mine), "done": cnt["done"], "running": cnt["running"],
                      "staging": cnt["staging"] + cnt["staged"], "failed": cnt["failed"],
                      "pending": cnt["pending"], "cpu_h_spent": f"{cpu_done[n]:.1f}",
                      "est_cpu_h_left": f"{left:.1f}", "h_since_last_event": f"{last:.2f}",
                      "runner": "up" if runner_alive else "down", "eta_h": eta,
                      "paused": "yes" if os.path.exists(f"{a.build}/nodes/{n}/PAUSE") else ""})
    tot = {k: sum(int(t[k]) for t in table) for k in ("planned", "done", "running", "staging", "failed", "pending")}
    cols = list(table[0].keys())
    print("\t".join(cols))
    for t in table:
        print("\t".join(str(t[c]) for c in cols))
    print("TOTAL\t" + "\t".join(str(tot[k]) for k in ("planned", "done", "running", "staging", "failed", "pending")))
    stale = [t["node"] for t in table if t["running"] and fnum(t["h_since_last_event"]) > 24]
    if stale:
        print(f"WARNING: no ledger event for >24 h on {', '.join(stale)} while samples run")
    if a.write:
        os.makedirs(f"{a.build}/status", exist_ok=True)
        shardlib.write_sheet(f"{a.build}/status/status.tsv", cols, table)
        srows = [{"cohort": p["cohort"], "sample_id": p["sample_id"], "node": p["node"],
                  **states[(p["cohort"], p["sample_id"])]} for p in plan]
        shardlib.write_sheet(f"{a.build}/status/samples.tsv", list(srows[0].keys()), srows)
    if a.costs or a.write:
        agg = defaultdict(lambda: {"samples": set(), "cpu_h": 0.0, "wall_h": 0.0, "max_rss_gb": 0.0})
        for r in costs:
            for key in ((r["node"], r["rule"]), (r["node"], "ALL"), ("ALL", r["rule"]), ("ALL", "ALL")):
                g = agg[key]
                g["samples"].add((r["cohort"], r["sample_id"]))
                g["cpu_h"] += fnum(r["cpu_time"]) / 3600
                g["wall_h"] += fnum(r["s"]) / 3600
                g["max_rss_gb"] = max(g["max_rss_gb"], fnum(r["max_rss"]) / 1024)
        crow = [{"node": k[0], "rule": k[1], "samples": len(v["samples"]), "cpu_h": f"{v['cpu_h']:.2f}",
                 "job_wall_h": f"{v['wall_h']:.2f}", "max_rss_gb": f"{v['max_rss_gb']:.2f}"}
                for k, v in sorted(agg.items())]
        if crow and a.write:
            shardlib.write_sheet(f"{a.build}/status/node_costs.tsv", list(crow[0].keys()), crow)
        if a.costs:
            for r in crow:
                if r["rule"] == "ALL":
                    print(f"cost\t{r['node']}\t{r['samples']} samples\t{r['cpu_h']} CPU-h\t"
                          f"{r['job_wall_h']} job-wall-h\tmax RSS {r['max_rss_gb']} GB")
    if a.require_complete and tot["done"] != tot["planned"]:
        print(f"not complete: {tot['done']}/{tot['planned']} done", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
