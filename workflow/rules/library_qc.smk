# ═══════════════════════════════════════════════════════════════════════════════
# LIBRARY QC  (per sample; NOT part of the calling fingerprint)
#
# Chimeric-library metrics that do not depend on CxSV discovery or scoring
# (scripts/library_qc.py, docs/LIBRARY_QC.md):
#   read level   supplementary-alignment fractions of primary MAPQ >= min_mapq reads
#                over 20 fixed 1-Mb windows of unique sequence
#                (config/library_qc_windows.hs1.bed)
#   call level   SVIM-only / single-caller BND/TRA and any-type records per Gb, from
#                merged_sv.vcf
#   flag         coverage-adjusted robust z of the SVIM-only rates against the panel
#                calibration of the sample's group (config/library_qc.yaml); FLAG when
#                any z > k
#
#   library_qc_params   {S}/qc/library_qc_params.json  (the config block used)
#   library_qc          {S}/qc/library_qc.tsv
#
# The group is the sample-sheet column qc_group, else patients.<sample>.qc_group,
# else library_qc.default_group[<platform>]. The flag changes no locus; with
# library_qc.conservative_mode true, score_cxsv.py applies the conservative Tier-1
# gate to a FLAGged library (rules/patient_scoring.smk).
# ═══════════════════════════════════════════════════════════════════════════════

import csv as _lqc_csv
import json as _lqc_json

LIBQC = config.get("library_qc", {}) or {}


def _lqc_path(p):
    return p if os.path.isabs(p) else os.path.join(PROJECT_DIR, p)


def library_qc_group(wildcards):
    """Calibration group: sample sheet qc_group, then patients.<sample>.qc_group, then
    library_qc.default_group by platform; NA when none applies."""
    spec = COHORT_SPECS.get(wildcards.cohort) or {}
    sheet = spec.get("samples") if isinstance(spec, dict) else None
    if sheet and os.path.exists(sheet):
        with open(sheet, newline="") as fh:
            for r in _lqc_csv.DictReader(fh, delimiter="\t"):
                if (r.get("sample_id") or "").strip() == wildcards.sample and (r.get("qc_group") or "").strip():
                    return r["qc_group"].strip()
    g = ((config.get("patients") or {}).get(wildcards.sample) or {}).get("qc_group")
    if g:
        return str(g)
    plat = (record(wildcards).platform or "").upper()
    return (LIBQC.get("default_group") or {}).get(plat, "NA")


rule library_qc_params:
    """The library QC configuration this sample is checked with."""
    output: f"{S}/qc/library_qc_params.json"
    params: body = lambda wc: _lqc_json.dumps({**LIBQC, "group": library_qc_group(wc)}, sort_keys=True)
    run:
        with open(output[0], "w") as fh:
            fh.write(params.body + "\n")


rule library_qc:
    """Chimeric-library metrics and flag of one sample (scripts/library_qc.py)."""
    input:
        unpack(caller_inputs),
        vcf     = f"{S}/merged_sv.vcf",
        params  = f"{S}/qc/library_qc_params.json",
        windows = _lqc_path(LIBQC.get("windows", "config/library_qc_windows.hs1.bed")),
        script  = f"{SCRIPTS}/library_qc.py"
    output: f"{S}/qc/library_qc.tsv"
    params:
        mapq  = int(LIBQC.get("min_mapq", 20)),
        group = library_qc_group
    log: f"{S}/logs/library_qc.log"
    shell:
        HTS_ENV + """
        python {input.script} --aln {input.aln} --ref {REFERENCE} --windows {input.windows} \
            --merged-vcf {input.vcf} --sample {wildcards.sample} --min-mapq {params.mapq} \
            --thresholds {input.params} --group {params.group} --out {output} 2> {log}
        """
