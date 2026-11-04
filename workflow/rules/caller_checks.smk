# ═══════════════════════════════════════════════════════════════════════════════
# CALLER OUTPUT CHECKS  (shared; NOT part of the calling fingerprint)
#
# Checks that run after the fingerprinted caller rules (rules/calling.smk) and
# fail a sample whose caller output is incomplete, without changing any caller
# rule or output (so no panel sample needs re-calling).
#
#   check_svim_output   {S}/qc/svim_check.tsv  (scripts/check_svim_output.py):
#                       SVIM log free of CRAM decode / truncation errors and
#                       SVIM calls on every contig where Sniffles2 has
#                       >= svim_check.min_reference_calls records
#
# Wired in as a per-sample target of the panel, panel_build and patient
# workflows; in the patient workflow score_cxsv also waits for it, so no score
# report is written from an incomplete SVIM call set.
# ═══════════════════════════════════════════════════════════════════════════════

SVIM_CHECK = config.get("svim_check", {}) or {}


rule check_svim_output:
    """Fail the sample when SVIM stopped early (decode error, contigs without calls)."""
    input:
        svim     = f"{S}/callers/svim/variants.vcf",
        sniffles = f"{S}/callers/sniffles.vcf",
        script   = f"{SCRIPTS}/check_svim_output.py"
    output: f"{S}/qc/svim_check.tsv"
    params:
        svim_log  = f"{S}/logs/svim.log",
        min_calls = int(SVIM_CHECK.get("min_reference_calls", 20))
    log: f"{S}/logs/check_svim_output.log"
    shell:
        "python {input.script} --svim-vcf {input.svim} --svim-log {params.svim_log} "
        "--reference-vcf {input.sniffles} --min-reference-calls {params.min_calls} "
        "--out {output} 2> {log}"
