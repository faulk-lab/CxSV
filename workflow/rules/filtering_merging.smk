# ═══════════════════════════════════════════════════════════════════════════════
# PHASE 3 — FILTERING & MERGING
# ═══════════════════════════════════════════════════════════════════════════════

rule filter_sniffles:
    """Filter Sniffles2 VCF: keep PASS calls with QUAL >= threshold."""
    input:  f"{RESULTS}/callers/{{sample}}_sniffles.vcf"
    output: f"{RESULTS}/callers/{{sample}}_sniffles_filtered.vcf"
    params: min_qual = config["filter_params"]["min_qual"]
    log:    f"logs/filter/{{sample}}_sniffles.log"
    shell:
        """
        bcftools view -f PASS -e 'QUAL < {params.min_qual}' \
            {input} -o {output} -O v > {log} 2>&1
        """


rule filter_cutesv:
    """Filter CuteSV VCF: keep PASS calls only."""
    input:  f"{RESULTS}/callers/{{sample}}_cutesv.vcf"
    output: f"{RESULTS}/callers/{{sample}}_cutesv_filtered.vcf"
    log:    f"logs/filter/{{sample}}_cutesv.log"
    shell:  "bcftools view -f PASS {input} -o {output} -O v > {log} 2>&1"


rule filter_svim:
    """
    Filter SVIM VCF: keep calls with QUAL >= threshold.
    SVIM uses QUAL as its confidence score (not the FILTER field),
    so the threshold is lower than for other callers.
    """
    input:  f"{RESULTS}/callers/{{sample}}_svim/variants.vcf"
    output: f"{RESULTS}/callers/{{sample}}_svim_filtered.vcf"
    params: min_qual = config["filter_params"]["svim_min_qual"]
    log:    f"logs/filter/{{sample}}_svim.log"
    shell:
        """
        bcftools view -e 'QUAL < {params.min_qual}' \
            {input} -o {output} -O v > {log} 2>&1
        """


rule merge_sv_callers:
    """
    Merge per-sample calls from all three callers with SURVIVOR.
    SUPP_VEC bit order: CuteSV=bit0, SVIM=bit1, Sniffles2=bit2.
    min_callers=1 keeps any call supported by at least one caller.
    """
    input:
        cutesv   = f"{RESULTS}/callers/{{sample}}_cutesv_filtered.vcf",
        svim     = f"{RESULTS}/callers/{{sample}}_svim_filtered.vcf",
        sniffles = f"{RESULTS}/callers/{{sample}}_sniffles_filtered.vcf"
    output:
        vcf      = f"{RESULTS}/per_sample/{{sample}}_merged_sv.vcf",
        filelist = f"{RESULTS}/per_sample/{{sample}}_vcf_list.txt"
    params:
        max_dist    = config["survivor_params"]["max_dist"],
        min_callers = config["survivor_params"]["min_callers"],
        min_sv_len  = config["sv_params"]["min_sv_len"]
    log: f"logs/merge/{{sample}}_merge.log"
    shell:
        """
        echo {input.cutesv}   >  {output.filelist}
        echo {input.svim}     >> {output.filelist}
        echo {input.sniffles} >> {output.filelist}
        SURVIVOR merge {output.filelist} {params.max_dist} \
            {params.min_callers} 1 0 0 {params.min_sv_len} {output.vcf} \
            > {log} 2>&1
        """


rule merge_all_samples:
    """
    Merge all per-sample VCFs into a single population-level VCF.
    Uses 'find' to collect all *_merged_sv.vcf files on disk so that
    this rule correctly includes all samples.
    """
    input:
        vcfs = expand(
            f"{RESULTS}/per_sample/{{sample}}_merged_sv.vcf", sample=SAMPLES
        )
    output:
        vcf      = f"{RESULTS}/final/long_read_filtered.vcf",
        filelist = f"{RESULTS}/final/all_samples_vcfs.txt"
    params:
        max_dist   = config["survivor_params"]["max_dist"],
        min_sv_len = config["sv_params"]["min_sv_len"]
    log: "logs/merge/all_samples.log"
    shell:
        """
        find {RESULTS}/per_sample -name "*_merged_sv.vcf" \
            | sort > {output.filelist}
        SURVIVOR merge {output.filelist} {params.max_dist} \
            1 1 0 0 {params.min_sv_len} {output.vcf} > {log} 2>&1
        echo "merge_all_samples: $(wc -l < {output.filelist}) samples merged" \
            >> {log}
        """


rule make_consensus_vcf:
    """
    Consensus VCF: variants supported by >= 2 of 3 callers.
    Higher-confidence call set — use this for downstream validation.
    """
    input:
        vcfs = expand(
            f"{RESULTS}/per_sample/{{sample}}_merged_sv.vcf", sample=SAMPLES
        )
    output:
        vcf      = f"{RESULTS}/final/consensus.vcf",
        filelist = f"{RESULTS}/final/consensus_vcf_list.txt"
    params:
        max_dist    = config["survivor_params"]["max_dist"],
        min_callers = config["survivor_params"]["min_callers_consensus"],
        min_sv_len  = config["sv_params"]["min_sv_len"]
    log: "logs/consensus_vcf.log"
    shell:
        """
        find {RESULTS}/per_sample -name "*_merged_sv.vcf" \
            | sort > {output.filelist}
        SURVIVOR merge {output.filelist} {params.max_dist} \
            {params.min_callers} 1 0 0 {params.min_sv_len} {output.vcf} \
            > {log} 2>&1
        """
