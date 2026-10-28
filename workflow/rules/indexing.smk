# ═══════════════════════════════════════════════════════════════════════════════
# REFERENCE CACHE, CRAM REFERENCE CHECK, AND INDEXING  (shared)
# Inputs may be BAM or CRAM; see workflow/lib/cxsv_inputs.py for discovery.
# ═══════════════════════════════════════════════════════════════════════════════

rule build_ref_cache:
    """
    MD5-keyed htslib cache of config['reference'], built once and shared by all
    workflows and cohorts. Lets every tool decode CRAMs on any node, whatever
    reference path is recorded in the CRAM header. ~3.1 GB for hs1.
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
        ok = f"{S}/refcheck.ok"
    log:
        f"{S}/logs/refcheck.log"
    shell:
        "python {SCRIPTS}/check_cram_reference.py {input.manifest} {input.cram} "
        "> {output.ok} 2> {log}"


# Index creation runs only in directory mode (bamfiles/<cohort>/): a sample
# sheet's indexes must already exist, so its data directories are never written.
rule index_bam:
    """samtools index for a BAM in bamfiles/<cohort>/ that has no .bai yet."""
    input:
        bam = f"{BAMFILES_ROOT}/{{cohort}}/{{sample}}.bam"
    output:
        bai = f"{BAMFILES_ROOT}/{{cohort}}/{{sample}}.bam.bai"
    threads: config["threads"]["samtools"]
    log:
        f"{BAMFILES_ROOT}/{{cohort}}/.{{sample}}.index.log"
    shell:
        "samtools index -@ {threads} {input.bam} > {log} 2>&1"


rule index_cram:
    """samtools index for a CRAM in bamfiles/<cohort>/ that has no .crai yet."""
    input:
        cram     = f"{BAMFILES_ROOT}/{{cohort}}/{{sample}}.cram",
        manifest = REF_MD5_MANIFEST
    output:
        crai = f"{BAMFILES_ROOT}/{{cohort}}/{{sample}}.cram.crai"
    threads: config["threads"]["samtools"]
    log:
        f"{BAMFILES_ROOT}/{{cohort}}/.{{sample}}.index.log"
    shell:
        HTS_ENV + "samtools index -@ {threads} {input.cram} > {log} 2>&1"
