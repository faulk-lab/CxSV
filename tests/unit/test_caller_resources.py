"""Caller memory requests and retries (rules/caller_resources.smk, roadmap item 16m)."""
import re

import yaml

import cxsv_fingerprint as fpm
from conftest import REPO


def test_caller_resources_outside_the_fingerprint():
    assert "rules/caller_resources.smk" not in fpm.CALLING_CODE
    assert "caller_resources" not in fpm.CALLING_PARAM_KEYS
    calling = (REPO / "workflow/rules/calling.smk").read_text()
    assert "mem_mb" not in calling and "retries" not in calling


def test_config_requests_memory_for_every_caller_and_one_svim_retry():
    cfg = yaml.safe_load((REPO / "config/calling.yaml").read_text())["caller_resources"]
    assert set(cfg) == {"svim", "sniffles", "cutesv"}
    assert all(int(v["mem_mb"]) > 0 for v in cfg.values())
    assert cfg["svim"]["retries"] == 1


def test_both_workflows_include_it_after_the_caller_rules():
    for wf in ("panel", "patient"):
        inc = re.findall(r'^include:\s*"\.\./rules/([\w.]+)"', (REPO / f"workflow/{wf}/Snakefile").read_text(), re.M)
        assert "caller_resources.smk" in inc and inc.index("caller_resources.smk") > inc.index("calling.smk"), wf


def test_rules_file_names_existing_caller_rules():
    smk = (REPO / "workflow/rules/caller_resources.smk").read_text()
    calling = (REPO / "workflow/rules/calling.smk").read_text()
    for rule in re.findall(r'"(run_\w+)"', smk):
        assert f"rule {rule}:" in calling, rule


def test_bundle_runner_passes_a_memory_budget():
    sh = (REPO / "container/run_cxsv.sh").read_text()
    assert '--resources mem_mb="$MEM_MB"' in sh and "--mem-mb)" in sh
    assert "intake|patient|panel)" in sh
