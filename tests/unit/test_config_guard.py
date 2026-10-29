"""Guards for the shared-calling-module guarantee (roadmap item 1).

* Both entry points include the same shared calling rule files.
* config/calling.yaml defines every fingerprinted section.
* No committed workflow or site config overrides a fingerprinted section, so a
  panel and a patient run started from this repository cannot silently differ.
"""
import re

import yaml

import cxsv_fingerprint as fpm
from conftest import REPO

SHARED = [
    "common.smk", "indexing.smk", "calling.smk", "per_sample.smk", "junctions.smk",
    "sample_report.smk", "fingerprint.smk",
]


def includes(snakefile):
    return re.findall(r'^include:\s*"\.\./rules/([\w.]+)"', snakefile.read_text(), re.M)


def test_both_workflows_share_the_calling_module():
    panel = includes(REPO / "workflow/panel/Snakefile")
    patient = includes(REPO / "workflow/patient/Snakefile")
    for rule_file in SHARED:
        assert rule_file in panel and rule_file in patient, rule_file
    assert panel[: len(SHARED)] == patient[: len(SHARED)] == SHARED


def test_calling_yaml_defines_all_fingerprinted_sections():
    cfg = yaml.safe_load((REPO / "config/calling.yaml").read_text())
    for key in fpm.CALLING_PARAM_KEYS:
        assert isinstance(cfg.get(key), dict) and cfg[key], key


def test_no_committed_config_overrides_fingerprinted_sections():
    paths = [REPO / "config/panel.yaml", REPO / "config/patient.yaml", *sorted((REPO / "config/sites").glob("*.yaml"))]
    assert paths
    for p in paths:
        cfg = yaml.safe_load(p.read_text()) or {}
        clash = [k for k in fpm.CALLING_PARAM_KEYS if k in cfg]
        assert not clash, f"{p.relative_to(REPO)} overrides fingerprinted section(s) {clash}"
