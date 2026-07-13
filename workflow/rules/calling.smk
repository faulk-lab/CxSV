# ═══════════════════════════════════════════════════════════════════════════════
# PHASE 2 — ENSEMBLE SV CALLING
# Three callers run independently on each sample.
# All take the local BAM + BAI as input — no streaming needed.
# ═══════════════════════════════════════════════════════════════════════════════

rule run_sniffles:
    """
    Sniffles2: long-read SV caller using split-read signatures.
    Outputs a VCF and a .snf (used for optional joint calling).
    """
    input:
        bam = f"{DATA}/{{sample}}.bam",
        bai = f"{DATA}/{{sample}}.bam.bai"
    output:
        vcf = f"{RESULTS}/callers/{{sample}}_sniffles.vcf",
        snf = f"{RESULTS}/callers/{{sample}}_sniffles.snf"
    params:
        min_sv_len  = config["sv_params"]["min_sv_len"],
        min_support = config["sv_params"]["min_support"]
    threads: config["threads"]["sv_callers"]
    log:
        f"logs/sniffles/{{sample}}.log"
    shell:
        """
        sniffles \
            --input      {input.bam} \
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
    Requires a local BAM + BAI file — works correctly here.
    Uses a per-sample working directory for temporary split-read files.
    """
    input:
        bam = f"{DATA}/{{sample}}.bam",
        bai = f"{DATA}/{{sample}}.bam.bai"
    output:
        vcf = f"{RESULTS}/callers/{{sample}}_cutesv.vcf"
    params:
        workdir     = f"{RESULTS}/callers/{{sample}}_cutesv_work",
        min_sv_len  = config["sv_params"]["min_sv_len"],
        min_support = config["sv_params"]["min_support"]
    threads: config["threads"]["sv_callers"]
    log:
        f"logs/cutesv/{{sample}}.log"
    shell:
        """
        mkdir -p {params.workdir}
        cuteSV \
            {input.bam} \
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
    Requires a local BAM + BAI file — works correctly here.
    Outputs to a per-sample directory; variants.vcf is the main output.
    """
    input:
        bam = f"{DATA}/{{sample}}.bam",
        bai = f"{DATA}/{{sample}}.bam.bai"
    output:
        vcf = f"{RESULTS}/callers/{{sample}}_svim/variants.vcf"
    params:
        outdir           = f"{RESULTS}/callers/{{sample}}_svim",
        min_sv_len       = config["sv_params"]["min_sv_len"],
        min_mapq         = config["sv_params"]["min_mapq"],
        cluster_max_dist = config["sv_params"]["cluster_max_distance"]
    log:
        f"logs/svim/{{sample}}.log"
    shell:
        """
        mkdir -p {params.outdir}
        svim alignment \
            {params.outdir} \
            {input.bam} \
            {REFERENCE} \
            --min_sv_size {params.min_sv_len} \
            --min_mapq {params.min_mapq} \
            --cluster_max_distance {params.cluster_max_dist} \
            > {log} 2>&1
        """
