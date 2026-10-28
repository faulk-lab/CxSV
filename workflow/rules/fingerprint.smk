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
    input:
        manifest = REF_MD5_MANIFEST,
        tools    = f"{RESULTS_ROOT}/tool_versions.tsv"
    output:
        f"{RESULTS_ROOT}/calling_fingerprint.json"
    run:
        fp = cxsv_fingerprint.build(config, WORKFLOW_DIR, input.manifest, input.tools)
        empty = [t for t, v in fp["tools"].items() if not v]
        if empty:
            raise ValueError(f"could not read a version for: {', '.join(empty)} "
                             f"(see {input.tools}); is the cxsv-pipeline env active?")
        with open(output[0], "w") as fh:
            json.dump(fp, fh, indent=2, sort_keys=True)
            fh.write("\n")
