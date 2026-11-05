#!/usr/bin/env python3
"""NEW_BUILD/aggregate.yaml for a downstream-only rebuild: OLD's config with the
new results root and release description (panel-v1.1 rebuild).

    make_aggregate_yaml.py OLD_BUILD NEW_BUILD DESCRIPTION_FILE
"""
import json
import os
import sys

old, new, desc = sys.argv[1:4]
cfg = json.load(open(os.path.join(old, "aggregate.yaml")))
cfg["panel_results_root"] = os.path.join(os.path.abspath(new), "results", "panel")
cfg["release"]["description"] = " ".join(open(desc).read().split())
cfg["release"]["date"] = None
tmp = os.path.join(new, "aggregate.yaml.tmp")
with open(tmp, "w") as fh:
    json.dump(cfg, fh, indent=2)
    fh.write("\n")
os.replace(tmp, os.path.join(new, "aggregate.yaml"))
print(json.dumps({k: cfg[k] for k in ("panel_results_root", "release")}, indent=1))
