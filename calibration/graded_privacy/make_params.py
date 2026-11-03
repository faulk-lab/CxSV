#!/usr/bin/env python3
"""scoring_params.json for one panel sample, as rule scoring_params writes it.

Weights, tiers, copy-number, clustering, fold-back settings and resources come
from config/scoring.yaml; matching from config/patient.yaml; collapse distance
and composite fraction from config/calling.yaml (all in this checkout). Patient
information: the sample's registry sex only (no HPO terms, no parents).
--weights-override (YAML) replaces top-level keys of `weights` (used to score
with a candidate privacy rule).
"""
import argparse
import json
import os
import sys

import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gp_common import REPO  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--release-json", required=True)
    ap.add_argument("--sex", default="")
    ap.add_argument("--weights-override", default=None)
    ap.add_argument("--config-dir", default=os.path.join(REPO, "config"),
                    help="config/ of the checkout whose scoring rules to use (default: this one)")
    ap.add_argument("--patient-from", default=None,
                    help="take the patient block from this scoring_params.json instead of --sex")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    cfg = {}
    for f in ("calling.yaml", "patient.yaml", "scoring.yaml"):
        with open(os.path.join(a.config_dir, f)) as fh:
            cfg.update(yaml.safe_load(fh) or {})
    W = cfg["weights"]
    if a.weights_override:
        with open(a.weights_override) as fh:
            W.update(yaml.safe_load(fh) or {})
    rel = json.load(open(a.release_json))
    m = cfg.get("panel_match", {})
    body = {
        "weights": W, "tiers": cfg["tiers"], "copy_number": cfg["copy_number"],
        "clustering": cfg["clustering"], "foldback": cfg["foldback"],
        "scoring_resources": {k: v for k, v in (cfg.get("scoring_resources") or {}).items() if v},
        "panel_match": {"tolerance": m.get("tolerance", 500), "ins_len_ratio": m.get("ins_len_ratio", 0.5),
                        "min_junction_fraction": m.get("min_junction_fraction", 0.75)},
        "collapse_distance": cfg["junction_params"]["collapse_distance"],
        "composite_fraction": cfg["cxsv_params"]["composite_fraction"],
        "panel_release_id": rel.get("release_id"),
        "patient": {"hpo_terms": [], "parents_affected": False,
                    **({"sex": a.sex} if a.sex in ("male", "female") else {})},
    }
    if a.patient_from:
        body["patient"] = json.load(open(a.patient_from))["patient"]
    tmp = a.out + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(body, fh, indent=2, sort_keys=True)
    os.replace(tmp, a.out)


if __name__ == "__main__":
    main()
