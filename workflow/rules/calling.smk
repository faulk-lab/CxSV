# ═══════════════════════════════════════════════════════════════════════════════
# ENSEMBLE SV CALLING  (shared — part of the calling fingerprint)
#
# Three callers run independently on each sample, from a local BAM or CRAM plus
# its index. For CRAMs, HTS_ENV points htslib at the MD5 reference cache and the
# sample must first pass check_cram_reference.
# Sniffles2 2.8.1, cuteSV 2.1.4 and SVIM 2.0.0 gave identical call sets from a
# matched BAM/CRAM pair (chr20, 0 discordant calls).
#
# Any change to this file changes the calling fingerprint, so a patient can no
# longer be compared with panels built before the change.
# ═══════════════════════════════════════════════════════════════════════════════

rule run_sniffles:
    """Sniffles2: split-read SV caller. Writes a VCF and a .snf for joint calling."""
    input:
        unpack(caller_inputs)
    output:
        vcf = f"{S}/callers/sniffles.vcf",
        snf = f"{S}/callers/sniffles.snf"
    params:
        min_sv_len  = config["sv_params"]["min_sv_len"],
        min_support = config["sv_params"]["min_support"]
    threads: config["threads"]["sv_callers"]
    log:
        f"{S}/logs/sniffles.log"
    shell:
        HTS_ENV + """
        sniffles \
            --input      {input.aln} \
            --reference  {REFERENCE} \
            --threads    {threads} \
            --sample-id  {wildcards.sample} \
            --vcf        {output.vcf} \
            --snf        {output.snf} \
            --minsvlen   {params.min_sv_len} \
            --minsupport {params.min_support} \
            > {log} 2>&1
        """


rule run_cutesv:
    """cuteSV: signature-based SV caller, with a per-sample working directory."""
    input:
        unpack(caller_inputs)
    output:
        vcf = f"{S}/callers/cutesv.vcf"
    params:
        workdir     = f"{S}/callers/cutesv_work",
        min_sv_len  = config["sv_params"]["min_sv_len"],
        min_support = config["sv_params"]["min_support"]
    threads: config["threads"]["sv_callers"]
    log:
        f"{S}/logs/cutesv.log"
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
    """SVIM: split-read and coverage-based SV caller; variants.vcf is the output."""
    input:
        unpack(caller_inputs)
    output:
        vcf = f"{S}/callers/svim/variants.vcf"
    params:
        outdir           = f"{S}/callers/svim",
        min_sv_len       = config["sv_params"]["min_sv_len"],
        min_mapq         = config["sv_params"]["min_mapq"],
        cluster_max_dist = config["sv_params"]["cluster_max_distance"]
    log:
        f"{S}/logs/svim.log"
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
