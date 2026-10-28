# ═══════════════════════════════════════════════════════════════════════════════
# COMMON — shared by the panel and patient workflows
#
# The including Snakefile must define, before `include:`-ing this file:
#   WORKFLOW_NAME  "panel" or "patient"
#   COHORT_SPECS   {name: spec}  — background cohorts (panel) or patient runs
#   RESULTS_ROOT   output root for this workflow, e.g. results/panel
#
# Everything per sample lives under {SAMPLES_ROOT}/{cohort}/{sample}/ so the
# same person can appear in two cohorts without output collisions.
# Sample discovery itself is in workflow/lib/cxsv_inputs.py (unit-tested).
# ═══════════════════════════════════════════════════════════════════════════════

import os
import re
import sys

sys.path.insert(0, os.path.join(workflow.basedir, "..", "lib"))
import cxsv_inputs            # noqa: E402
import cxsv_fingerprint       # noqa: E402

WORKFLOW_DIR  = os.path.abspath(os.path.join(workflow.basedir, ".."))
SCRIPTS       = config.get("scripts_dir", "workflow/scripts")
REFERENCE     = config["reference"]
BAMFILES_ROOT = config.get("bamfiles_root", "bamfiles")
SAMPLES_ROOT  = f"{RESULTS_ROOT}/samples"
S             = f"{SAMPLES_ROOT}/{{cohort}}/{{sample}}"      # per-sample output prefix

if not COHORT_SPECS:
    raise ValueError(f"[{WORKFLOW_NAME}] no cohorts/runs configured -- see config/{WORKFLOW_NAME}.yaml")

SAMPLE_TABLE, _notes = cxsv_inputs.build_sample_table(COHORT_SPECS, BAMFILES_ROOT)
for _n in _notes:
    print(f"[{WORKFLOW_NAME}] {_n}", file=sys.stderr)
if not SAMPLE_TABLE:
    raise ValueError(f"[{WORKFLOW_NAME}] no samples found in any configured cohort/run")

COHORTS = sorted(COHORT_SPECS)
KEYS    = sorted(SAMPLE_TABLE)                     # [(cohort, sample), ...]
COHORT_OF_KEY = [k[0] for k in KEYS]
SAMPLE_OF_KEY = [k[1] for k in KEYS]


def per_sample(pattern, cohort=None):
    """Expand a per-sample path pattern over every (cohort, sample) pair,
    or over one cohort's samples."""
    keys = [k for k in KEYS if cohort is None or k[0] == cohort]
    return [pattern.format(cohort=c, sample=s) for c, s in keys]


# Pin wildcards to known values so sample IDs containing dots can never be
# split ambiguously against other patterns.
wildcard_constraints:
    cohort="|".join(re.escape(c) for c in COHORTS),
    sample="|".join(re.escape(s) for s in sorted(set(SAMPLE_OF_KEY))),

# ─── Reference cache (shared by every workflow and cohort) ───────────────────
# A CRAM stores reads as differences from the reference; htslib finds the
# sequence by each @SQ line's M5 tag via REF_PATH. The cache is built once
# (rule build_ref_cache) and every rule that reads alignments exports HTS_ENV.
# With REF_PATH set to the local cache only, htslib never downloads anything.
REF_CACHE        = os.path.abspath(config.get("ref_cache", "resources/ref_cache"))
REF_MD5_MANIFEST = f"{REF_CACHE}/md5_manifest.tsv"
HTS_ENV = f"export REF_PATH='{REF_CACHE}/%2s/%2s/%s' REF_CACHE='{REF_CACHE}/%2s/%2s/%s'; "


# ─── Input functions ─────────────────────────────────────────────────────────
def record(wildcards):
    return SAMPLE_TABLE[(wildcards.cohort, wildcards.sample)]


def alignment(wildcards):
    return record(wildcards).alignment


def caller_inputs(wildcards):
    """aln + idx for every sample; CRAMs also wait for their reference check."""
    r = record(wildcards)
    inputs = {"aln": r.alignment, "idx": r.index}
    if r.is_cram:
        inputs["refcheck"] = f"{SAMPLES_ROOT}/{r.cohort}/{r.sample_id}/refcheck.ok"
    return inputs
