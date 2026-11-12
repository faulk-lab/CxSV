#!/usr/bin/env python3
"""Calling fingerprint of a code tree with the panel-v1.2-r2 reference manifest and tool
versions (workflow/lib/cxsv_fingerprint.py of that tree, as rule calling_fingerprint).
  fingerprint_of.py TREE [TREE ...]"""
import importlib.util, json, sys
import yaml
REL = "/mnt/sargasso1/CxSV/build/panel-v1.2-r2/results/panel/release"
MANIFEST = "/mnt/sargasso1/CxSV/bundles/cxsv-bundle-20261109-84d380e08cc4/reference/ref_cache/md5_manifest.tsv"
rel = json.load(open(f"{REL}/calling_fingerprint.json"))
for tree in sys.argv[1:]:
    spec = importlib.util.spec_from_file_location("fp", f"{tree}/workflow/lib/cxsv_fingerprint.py")
    fp = importlib.util.module_from_spec(spec); spec.loader.exec_module(fp)
    cfg = {}
    for f in ("calling.yaml",):
        cfg.update(yaml.safe_load(open(f"{tree}/config/{f}")) or {})
    d = fp.build(cfg, f"{tree}/workflow", MANIFEST, f"{REL}/tool_versions.tsv")
    print(tree, fp.digest(d)[:12], "full", fp.digest(d), "| diffs vs release:", fp.compare(d, rel)[:6])
print("release", rel.get("digest") or rel.get("fingerprint") or "", sorted(rel)[:8])
