#!/usr/bin/env python3
"""Write release.json: what a panel release is, where it came from, how it was made.

Read by people, by the patient workflow (callability definition, release id) and
by the bundle builder. Every value is derived from the release's own files, the
calling fingerprint, the panel config and the runtime; nothing is typed in.

Runtime: when the workflow runs inside the CxSV container, the image sets
CXSV_IMAGE (tag) and CXSV_ENV_LOCK_SHA256 (sha256 of the explicit conda lock it
was built from); the bundle's run script adds CXSV_IMAGE_ID (Docker image id)
and CXSV_SIF_SHA256. Outside a container these are null and CONDA_PREFIX is
recorded instead.
"""
import argparse
import csv
import datetime
import json
import os
import platform as pyplatform
import subprocess
import sys
from collections import Counter

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import cxsv_release as cr  # noqa: E402
import cxsv_fingerprint as fpm  # noqa: E402

SCHEMA = "cxsv-panel-release/1"


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def version_of(cmd):
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        text = (out.stdout or out.stderr).strip().splitlines()
        return text[0] if text else None
    except (OSError, subprocess.SubprocessError):
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--release-id-file", required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--test-release", action="store_true")
    ap.add_argument("--description", default="")
    ap.add_argument("--fingerprint", required=True)
    ap.add_argument("--sample-manifest", required=True)
    ap.add_argument("--panel-individuals", required=True)
    ap.add_argument("--sites", required=True)
    ap.add_argument("--loci", required=True)
    ap.add_argument("--callability-summary", required=True)
    ap.add_argument("--min-depth", type=int, required=True)
    ap.add_argument("--min-mapq", type=int, required=True)
    ap.add_argument("--flank", type=int, required=True)
    ap.add_argument("--site-tolerance", type=int, required=True)
    ap.add_argument("--site-ins-len-ratio", type=float, required=True)
    ap.add_argument("--tool-versions", required=True)
    ap.add_argument("--project-dir", required=True)
    ap.add_argument("--files", nargs="+", required=True, help="release files (relative names) besides release.json")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    rid = open(a.release_id_file).read().strip()
    with open(a.fingerprint) as fh:
        fp = json.load(fh)
    samples = read_tsv(a.sample_manifest)
    inds = read_tsv(a.panel_individuals)
    n_sites = sum(1 for _ in open(a.sites)) - 1
    n_loci = sum(1 for _ in open(a.loci)) - 1
    summ = {}
    per_sample_bp = {}
    for r in read_tsv(a.callability_summary):
        if r["level"] == "sample_callable_bp":
            per_sample_bp[r["name"]] = int(r["value"])
        elif r["level"] == "sites":
            summ[r["name"]] = int(r["value"])
    ind_platform = Counter(p for i in inds for p in i["platforms"].split(","))

    meta = {
        "schema": SCHEMA,
        "release_id": rid,
        "name": a.name,
        "date": rid[len(a.name) + 1: len(a.name) + 9],
        "test_release": bool(a.test_release),
        "description": a.description,
        "created_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "source": cr.source_info(a.project_dir),
        "calling_fingerprint": {"sha256": fp["sha256"], "file": "calling_fingerprint.json",
                                "code_sha256": fp["code"]["sha256"], "params": fp["params"]},
        "reference": {"n_contigs": len(fp["reference"]),
                      "total_length": sum(c["length"] for c in fp["reference"]),
                      "contig_md5": {c["name"]: c["md5"] for c in fp["reference"]}},
        "panel": {
            "n_samples": len(samples),
            "n_individuals": len(inds),
            "samples_by_cohort": dict(sorted(Counter(s["cohort"] for s in samples).items())),
            "samples_by_platform": dict(sorted(Counter(s["platform"] for s in samples).items())),
            "individuals_by_platform": dict(sorted(ind_platform.items())),
            "individuals_by_superpopulation": dict(sorted(Counter(i["superpopulation"] for i in inds).items())),
            "individuals_by_sex": dict(sorted(Counter(i["sex"] for i in inds).items())),
        },
        "contents": {"n_junction_sites": n_sites, "n_cxsv_loci": n_loci},
        "junction_sites": {"tolerance_bp": a.site_tolerance, "ins_len_ratio": a.site_ins_len_ratio,
                           "frequencies": "ac/an/af overall, af_<platform>, af_<superpopulation> "
                                          "(panel_junction_sites.tsv); carrier individuals per site "
                                          "(panel_junction_carriers.tsv)"},
        "callability": {
            "definition": ("a junction is callable in a panel sample when each of its two breakpoints has "
                           f"depth >= {a.min_depth} (reads of MAPQ >= {a.min_mapq}) over a {a.flank}-bp "
                           "window on at least one side; an individual is callable when any of its samples is"),
            "min_depth": a.min_depth, "min_mapq": a.min_mapq, "flank_bp": a.flank,
            "masks": "callability/samples.tsv + callability/<cohort>__<sample>.callable.bed.gz",
            "per_site": "panel_junction_callability.tsv",
            "sites": summ,
            "callable_bp_per_sample": per_sample_bp,
        },
        "tools": {**fpm.read_tool_versions(a.tool_versions),
                  "mosdepth": version_of(["mosdepth", "--version"]),
                  "snakemake": version_of(["snakemake", "--version"]),
                  "python": pyplatform.python_version()},
        "runtime": {
            "container_image": os.environ.get("CXSV_IMAGE"),
            "image_id": os.environ.get("CXSV_IMAGE_ID"),
            "sif_sha256": os.environ.get("CXSV_SIF_SHA256"),
            "env_lock_sha256": os.environ.get("CXSV_ENV_LOCK_SHA256"),
            "conda_prefix": None if os.environ.get("CXSV_IMAGE") else os.environ.get("CONDA_PREFIX"),
        },
        "files": sorted(set(a.files) | {cr.METADATA}),
        "integrity": f"{cr.MANIFEST}: file, bytes, sha256 of every file above",
    }
    with open(a.out, "w") as fh:
        json.dump(meta, fh, indent=2, sort_keys=False)
        fh.write("\n")


if __name__ == "__main__":
    main()
