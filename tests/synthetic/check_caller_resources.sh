#!/usr/bin/env bash
# =============================================================================
# Caller memory requests and retries (rules/caller_resources.smk, item 16m).
#
#   bash tests/synthetic/check_caller_resources.sh SYNTH_WORKDIR
#
# Uses the site.yaml written by run_synthetic.sh. A forced dry run of the
# patient workflow must show the config's mem_mb on run_svim / run_sniffles /
# run_cutesv and report retries=1 for SVIM, and the calling fingerprint's code
# digest must not depend on rules/caller_resources.smk or the caller_resources
# config (both are outside the fingerprint).
# =============================================================================
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
W="$(cd "${1:?usage: check_caller_resources.sh SYNTH_WORKDIR}" && pwd)"
cd "$REPO"
snakemake -s workflow/patient/Snakefile --configfile "$W/site.yaml" -n --forceall --cores 4 \
    --resources mem_mb=100000 > "$W/caller_resources_dryrun.log" 2>&1
python - "$W/caller_resources_dryrun.log" <<'PY'
import re, sys, yaml
log = open(sys.argv[1]).read()
cfg = yaml.safe_load(open("config/calling.yaml"))["caller_resources"]
fail = 0
for caller, rule in (("svim", "run_svim"), ("sniffles", "run_sniffles"), ("cutesv", "run_cutesv")):
    blocks = re.findall(rf"(?:local)?rule {rule}:\n((?:    [^\n]*\n)+)", log)   # the indented job block
    mems = {m for b in blocks for m in re.findall(r"mem_mb=(\d+)", b)}
    ok = blocks and mems == {str(cfg[caller]["mem_mb"])}
    print(("PASS" if ok else "FAIL"), f"{rule}: mem_mb {sorted(mems)} == {cfg[caller]['mem_mb']} in {len(blocks)} job(s)")
    fail += not ok
ok = re.search(r"\[caller_resources\] run_svim: mem_mb=\d+ x attempt, retries=1", log) is not None
print(("PASS" if ok else "FAIL"), "run_svim gets retries=1")
fail += not ok
sys.exit(1 if fail else 0)
PY
python - <<'PY'
import sys
sys.path.insert(0, "workflow/lib")
import cxsv_fingerprint as fpm
assert "rules/caller_resources.smk" not in fpm.CALLING_CODE and "caller_resources" not in fpm.CALLING_PARAM_KEYS
calling = open("workflow/rules/calling.smk").read()
assert "mem_mb" not in calling and "retries" not in calling, "caller resources belong in caller_resources.smk"
print("PASS caller resources are outside the calling fingerprint")
PY
