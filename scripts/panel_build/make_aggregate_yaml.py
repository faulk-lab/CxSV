#!/usr/bin/env python3
"""NEW_BUILD/aggregate.yaml for a downstream-only rebuild: OLD's config with the
new results root and release description (panel-v1.1 rebuild), and, for a build
with the copy-number track (panel-v1.2, roadmap item 15), the run-configuration
keys rebuild_downstream.sh requires.

    make_aggregate_yaml.py OLD_BUILD NEW_BUILD DESCRIPTION_FILE \
        [--cn-depth-precomputed DIR] [--cn-strata-table TSV]

--cn-depth-precomputed: root of a full-panel depth pass (<DIR>/<cohort>/<sample>/
depth.regions.bed.gz, depth.ok, depth.sha256; calibration/cn_linkage/cn_depth_build.sh).
--cn-strata-table: TSV cohort, sample_id, cn_stratum (calibration/cn_linkage/panel_strata_v1.tsv).
Both must exist; they are written as absolute paths. The file is JSON (a YAML subset).
"""
import argparse
import json
import os
import sys


def build(old, new, desc, cn_depth=None, cn_strata=None):
    cfg = json.load(open(os.path.join(old, "aggregate.yaml")))
    cfg["panel_results_root"] = os.path.join(os.path.abspath(new), "results", "panel")
    cfg["release"]["description"] = " ".join(open(desc).read().split())
    cfg["release"]["date"] = None
    for key, path in (("cn_depth_precomputed", cn_depth), ("cn_strata_table", cn_strata)):
        if path is None:
            continue
        if not os.path.exists(path):
            raise SystemExit(f"--{key.replace('_', '-')}: {path} does not exist")
        cfg[key] = os.path.abspath(path)
    tmp = os.path.join(new, "aggregate.yaml.tmp")
    with open(tmp, "w") as fh:
        json.dump(cfg, fh, indent=2)
        fh.write("\n")
    os.replace(tmp, os.path.join(new, "aggregate.yaml"))
    return cfg


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("old_build")
    ap.add_argument("new_build")
    ap.add_argument("description_file")
    ap.add_argument("--cn-depth-precomputed", default=None)
    ap.add_argument("--cn-strata-table", default=None)
    a = ap.parse_args(argv)
    cfg = build(a.old_build, a.new_build, a.description_file, a.cn_depth_precomputed, a.cn_strata_table)
    keys = ("panel_results_root", "release", "cn_depth_precomputed", "cn_strata_table")
    print(json.dumps({k: cfg[k] for k in keys if k in cfg}, indent=1))


if __name__ == "__main__":
    sys.exit(main())
