# ═══════════════════════════════════════════════════════════════════════════════
# CALLING FINGERPRINT  (shared)
#
# Records what must match between a panel and a patient run: calling and
# discovery parameters, reference contigs (name/length/MD5), tool versions and
# the code that defines CxSV loci. See workflow/lib/cxsv_fingerprint.py.
# The panel release carries this file; the patient workflow compares its own
# against it before it labels anything private (rules/patient_private.smk).
# ═══════════════════════════════════════════════════════════════════════════════

import json

rule tool_versions:
    output: f"{RESULTS_ROOT}/tool_versions.tsv"
    shell:  "bash {SCRIPTS}/tool_versions.sh > {output}"


rule calling_fingerprint:
    """The fingerprinted code files and config/calling.yaml are inputs and the
    calling parameters are params, so the fingerprint is rebuilt whenever they
    change, also in an existing work directory or build tree (roadmap item 13:
    without them an old calling_fingerprint.json was kept after a code change,
    so a release or a patient compatibility check could carry a stale value)."""
    input:
        manifest = REF_MD5_MANIFEST,
        tools    = f"{RESULTS_ROOT}/tool_versions.tsv",
        code     = [os.path.join(WORKFLOW_DIR, f) for f in cxsv_fingerprint.CALLING_CODE],
        calling  = os.path.join(PROJECT_DIR, "config", "calling.yaml")
    output:
        f"{RESULTS_ROOT}/calling_fingerprint.json"
    params:
        calling = lambda wc: json.dumps(cxsv_fingerprint.calling_params(config), sort_keys=True)
    run:
        fp = cxsv_fingerprint.build(config, WORKFLOW_DIR, input.manifest, input.tools)
        empty = [t for t, v in fp["tools"].items() if not v]
        if empty:
            raise ValueError(f"could not read a version for: {', '.join(empty)} "
                             f"(see {input.tools}); is the cxsv-pipeline env active?")
        with open(output[0], "w") as fh:
            json.dump(fp, fh, indent=2, sort_keys=True)
            fh.write("\n")
