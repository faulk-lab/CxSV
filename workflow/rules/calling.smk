# ═══════════════════════════════════════════════════════════════════════════════
# PHASE 2 — ENSEMBLE SV CALLING
# Three callers run independently on each sample, from a local BAM or CRAM plus
# its index. For CRAMs, HTS_ENV (rules/common.smk) points htslib at the MD5
# reference cache and the sample must first pass check_cram_reference.
# Sniffles2 2.8.1, cuteSV 2.1.4 and SVIM 2.0.0 gave identical call sets from a
# matched BAM/CRAM pair (chr20: 511 / 2,064 / 10,577 calls, 0 discordant).
# ═══════════════════════════════════════════════════════════════════════════════

rule run_sniffles:
    """
    Sniffles2: long-read SV caller using split-read signatures.
    Outputs a VCF and a .snf (used for optional joint calling).
    """
    input:
        unpack(caller_inputs)
    output:
        vcf = f"{RESULTS}/callers/{{sample}}_sniffles.vcf",
        snf = f"{RESULTS}/callers/{{sample}}_sniffles.snf"
    params:
        min_sv_len  = config["sv_params"]["min_sv_len"],
        min_support = config["sv_params"]["min_support"]
    threads: config["threads"]["sv_callers"]
    log:
        f"{RESULTS}/logs/sniffles/{{sample}}.log"
    shell:
        HTS_ENV + """
        sniffles \
            --input      {input.aln} \
            --reference  {REFERENCE} \
            --threads    {threads} \
            --vcf        {output.vcf} \
            --snf        {output.snf} \
            --minsvlen   {params.min_sv_len} \
            --minsupport {params.min_support} \
            > {log} 2>&1
        """


rule run_cutesv:
    """
    CuteSV: signature-based SV caller.
    Uses a per-sample working directory for temporary split-read files.
    """
    input:
        unpack(caller_inputs)
    output:
        vcf = f"{RESULTS}/callers/{{sample}}_cutesv.vcf"
    params:
        workdir     = f"{RESULTS}/callers/{{sample}}_cutesv_work",
        min_sv_len  = config["sv_params"]["min_sv_len"],
        min_support = config["sv_params"]["min_support"]
    threads: config["threads"]["sv_callers"]
    log:
        f"{RESULTS}/logs/cutesv/{{sample}}.log"
    shell:
        HTS_ENV + """
        mkdir -p {params.workdir}
        cuteSV \
            {input.aln} \
            {REFERENCE} \
            {output.vcf} \
            {params.workdir} \
            -t {threads} \
            -s {params.min_support} \
            -l {params.min_sv_len} \
            --genotype \
            > {log} 2>&1
        """


rule run_svim:
    """
    SVIM: split-read and coverage-based SV caller.
    Outputs to a per-sample directory; variants.vcf is the main output.
    """
    input:
        unpack(caller_inputs)
    output:
        vcf = f"{RESULTS}/callers/{{sample}}_svim/variants.vcf"
    params:
        outdir           = f"{RESULTS}/callers/{{sample}}_svim",
        min_sv_len       = config["sv_params"]["min_sv_len"],
        min_mapq         = config["sv_params"]["min_mapq"],
        cluster_max_dist = config["sv_params"]["cluster_max_distance"]
    log:
        f"{RESULTS}/logs/svim/{{sample}}.log"
    shell:
        HTS_ENV + """
        mkdir -p {params.outdir}
        svim alignment \
            {params.outdir} \
            {input.aln} \
            {REFERENCE} \
            --min_sv_size {params.min_sv_len} \
            --min_mapq {params.min_mapq} \
            --cluster_max_distance {params.cluster_max_dist} \
            > {log} 2>&1
        """
