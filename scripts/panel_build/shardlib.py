"""Helpers shared by the panel build scripts (scripts/panel_build/*.py).

Plain Python 3 standard library only, so the scripts run on any lab node's
system python as well as inside the pipeline image. YAML is read with PyYAML
when it is importable and otherwise with a small reader that covers the
subset used by config/sites/*.yaml (nested mappings, scalars, comments).
Overlays are written as JSON, which is valid YAML.
"""
from __future__ import annotations

import csv
import json
import os


def _scalar(v: str):
    v = v.strip()
    if not v or v in ("null", "~"):
        return None
    if v[0] in "\"'" and v[-1] == v[0]:
        return v[1:-1]
    if v in ("true", "false"):
        return v == "true"
    for t in (int, float):
        try:
            return t(v)
        except ValueError:
            pass
    return v


def _mini_yaml(text: str) -> dict:
    root: dict = {}
    stack = [(-1, root)]
    for raw in text.splitlines():
        line = raw.split(" #")[0].rstrip() if not raw.lstrip().startswith("#") else ""
        if not line.strip():
            continue
        ind = len(line) - len(line.lstrip())
        key, _, val = line.strip().partition(":")
        while stack[-1][0] >= ind:
            stack.pop()
        parent = stack[-1][1]
        if val.strip() in ("", "{}"):
            parent[key] = {}
            if val.strip() == "":
                stack.append((ind, parent[key]))
        else:
            parent[key] = _scalar(val)
    return root


def load_yaml(path: str) -> dict:
    text = open(path).read()
    try:
        import yaml  # type: ignore
        return yaml.safe_load(text) or {}
    except ImportError:
        return _mini_yaml(text)


def read_sheet(path: str):
    with open(path, newline="") as fh:
        rd = csv.DictReader((ln.replace("\r", "") for ln in fh), delimiter="\t")
        rows = list(rd)
        return list(rd.fieldnames or []), rows


def write_sheet(path: str, header, rows) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in header})
    os.replace(tmp, path)


# Site keys copied into every node / aggregate overlay unchanged. Fingerprinted
# calling sections are deliberately absent: they come from config/calling.yaml
# only, identically on every node (tests/unit/test_config_guard.py).
PASSTHROUGH = ("reference", "ref_cache", "annotation", "panel_individuals", "release",
               "panel_sites", "sniffles_population")


def node_overlay(site: dict, cohorts: dict, results_root: str, threads: int,
                 callability_threads: int = 4) -> dict:
    """A self-contained overlay: the site's paths plus this node's cohort sheets.

    The site overlay itself is NOT passed to the node's snakemake: its
    panel_cohorts would merge with the node's and bring every sample back."""
    ov = {k: site[k] for k in PASSTHROUGH if site.get(k) is not None}
    ov["panel_results_root"] = results_root
    ov["panel_cohorts"] = cohorts
    ov["threads"] = {"sv_callers": int(threads), "samtools": int(threads)}
    ov["panel_callability"] = dict(site.get("panel_callability") or {}, threads=int(callability_threads))
    return ov


def dump_json(path: str, obj) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, indent=2)
        fh.write("\n")
    os.replace(tmp, path)


def staged_path(path: str, stage_dir: str, cohort: str) -> str:
    return f"{stage_dir}/{cohort}/{os.path.basename(path)}"


def staged_row(row: dict, stage_dir: str, cohort: str) -> dict:
    """Sheet row with alignment and index pointing at node-local copies."""
    r = dict(row)
    aln = r["alignment"]
    idx = r.get("index") or aln + (".crai" if aln.endswith(".cram") else ".bai")
    r["alignment"] = staged_path(aln, stage_dir, cohort)
    r["index"] = staged_path(idx, stage_dir, cohort)
    return r
