#!/usr/bin/env python3
"""Write a node-style config overlay that lists a subset of the registry samples.

    cost_overlay.py --site config/sites/faulk_lab_panel.yaml --root ROOT \
        --threads 4 [--stage-dir DIR] [--tag NAME] COHORT/SAMPLE [...] > overlay.yaml

Used for the per-sample cost measurement (scripts/panel_build/measure_cost.sh):
the overlay points panel_cohorts at one-cohort sheets under ROOT/sheets/ that
hold only the named rows of the site's sheets, and panel_results_root at
ROOT/results/panel. Same format as the node overlays of make_shards.py.
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import shardlib  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--site", required=True)
    ap.add_argument("--root", required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--stage-dir", default=None)
    ap.add_argument("--tag", default="all")
    ap.add_argument("keys", nargs="+")
    a = ap.parse_args()
    site = shardlib.load_yaml(a.site)
    want = {}
    for k in a.keys:
        c, s = k.split("/", 1)
        want.setdefault(c, set()).add(s)
    os.makedirs(f"{a.root}/sheets", exist_ok=True)
    cohorts = {}
    for c, ids in want.items():
        spec = site["panel_cohorts"][c]
        header, rows = shardlib.read_sheet(spec["samples"])
        keep = [r for r in rows if r["sample_id"] in ids]
        missing = ids - {r["sample_id"] for r in keep}
        if missing:
            sys.exit(f"{c}: not in {spec['samples']}: {sorted(missing)}")
        if a.stage_dir:
            keep = [shardlib.staged_row(r, a.stage_dir, c) for r in keep]
        path = f"{a.root}/sheets/{a.tag}/{c}.tsv"
        os.makedirs(os.path.dirname(path), exist_ok=True)
        shardlib.write_sheet(path, header, keep)
        cohorts[c] = {"samples": path, "platform": spec.get("platform", "unknown")}
    ov = shardlib.node_overlay(site, cohorts, f"{a.root}/results/panel", a.threads)
    json.dump(ov, sys.stdout, indent=2)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
