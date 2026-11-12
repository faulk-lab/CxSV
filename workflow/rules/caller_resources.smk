# ═══════════════════════════════════════════════════════════════════════════════
# CALLER MEMORY AND RETRIES  (shared; NOT part of the calling fingerprint)
#
# Gives the three fingerprinted caller rules of rules/calling.smk an explicit
# memory request and a retry count, from config caller_resources (calling.yaml,
# outside the fingerprinted sections), without editing calling.smk: resources
# and retries decide when a job runs, never what it calls, so a panel built
# before this file keeps its fingerprint (roadmap item 16m).
#
# Why: SVIM 2.0.0 aborted with std::bad_alloc on athena (94 GB) on the 39x
# BOC_3559 genome on 2026-11-11 (patient run, 24 cores); the rerun on freya
# passed. Run alone on the same CRAM on 2026-11-12 SVIM peaked at 5.4 GB and
# passed on athena and on betsy, so the abort was not SVIM's own footprint and
# did not reproduce. A memory request keeps the callers' sum within the node,
# and one retry covers a transient failure (docs/RELEASE_AND_BUNDLE.md
# "Node memory").
#
#   mem_mb    reserved per job. Effective only with a memory budget
#             (snakemake --resources mem_mb=N; run_cxsv.sh passes 90% of the
#             host's RAM, override with --mem-mb). Jobs then never run together
#             beyond N; a request above N is capped at N by snakemake.
#   retries   re-runs after a failure; attempt k reserves mem_mb x k, so a retry
#             runs with less beside it.
# Must be included after rules/calling.smk.
# ═══════════════════════════════════════════════════════════════════════════════

CALLER_RESOURCES = config.get("caller_resources") or {}
_CALLER_RULES = {"svim": "run_svim", "sniffles": "run_sniffles", "cutesv": "run_cutesv"}


def _mem_by_attempt(base_mb):
    return lambda wildcards, attempt: int(base_mb) * attempt


for _caller, _spec in CALLER_RESOURCES.items():
    if _caller not in _CALLER_RULES:
        raise ValueError(f"caller_resources: unknown caller '{_caller}' (known: {sorted(_CALLER_RULES)})")
    _rule = getattr(rules, _CALLER_RULES[_caller]).rule
    if _spec.get("mem_mb"):
        _rule.resources["mem_mb"] = _mem_by_attempt(_spec["mem_mb"])
    if _spec.get("retries") is not None:
        _rule.restart_times = int(_spec["retries"])
    print(f"[caller_resources] {_rule.name}: mem_mb={_spec.get('mem_mb')} x attempt, "
          f"retries={_spec.get('retries', 0)}", file=sys.stderr)
