# ═══════════════════════════════════════════════════════════════════════════════
# PHASE 1 — REFERENCE CACHE, CRAM REFERENCE CHECK, AND INDEXING
# Inputs may be BAM or CRAM (see rules/common.smk for sample discovery).
# ═══════════════════════════════════════════════════════════════════════════════

rule build_ref_cache:
    """
    MD5-keyed htslib cache of config['reference'], built once and shared by all
    cohorts. Lets every tool decode CRAMs on any node, regardless of the
    reference path recorded in the CRAM header. ~3.1 GB for hs1.
    """
    input:
        ref = REFERENCE
    output:
        manifest = REF_MD5_MANIFEST
    log:
        f"{REF_CACHE}.build.log"
    shell:
        "python {SCRIPTS}/build_ref_cache.py {input.ref} {REF_CACHE} > {log} 2>&1"


rule check_cram_reference:
    """
    Verify a CRAM's @SQ M5 tags against the reference cache before any caller
    reads it. Seconds per sample; stops a wrong-assembly CRAM up front.
    """
    input:
        cram     = alignment,
        manifest = REF_MD5_MANIFEST
    output:
        ok = f"{RESULTS}/refcheck/{{sample}}.ok"
    log:
        f"{RESULTS}/logs/refcheck/{{sample}}.log"
    shell:
        "python {SCRIPTS}/check_cram_reference.py {input.manifest} {input.cram} "
        "> {output.ok} 2> {log}"


# Index creation runs only in directory mode (bamfiles/<cohort>/): a sample
# sheet's indexes must already exist, so its data directories are never written.
rule index_bam:
    """samtools index for a BAM in bamfiles/<cohort>/ that has no .bai yet."""
    input:
        bam = f"{DATA}/{{sample}}.bam"
    output:
        bai = f"{DATA}/{{sample}}.bam.bai"
    threads: config["threads"]["samtools"]
    log:
        f"{RESULTS}/logs/index/{{sample}}.log"
    shell:
        "samtools index -@ {threads} {input.bam} > {log} 2>&1"


rule index_cram:
    """samtools index for a CRAM in bamfiles/<cohort>/ that has no .crai yet."""
    input:
        cram     = f"{DATA}/{{sample}}.cram",
        manifest = REF_MD5_MANIFEST
    output:
        crai = f"{DATA}/{{sample}}.cram.crai"
    threads: config["threads"]["samtools"]
    log:
        f"{RESULTS}/logs/index/{{sample}}.log"
    shell:
        HTS_ENV + "samtools index -@ {threads} {input.cram} > {log} 2>&1"
