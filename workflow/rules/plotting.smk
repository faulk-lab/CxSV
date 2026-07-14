# ═══════════════════════════════════════════════════════════════════════════════
# PHASE 7 — VISUALIZATION
# ═══════════════════════════════════════════════════════════════════════════════

rule plot_sv_count_by_type:
    """Bar plot: SV count by type (DEL/DUP/INV/INS/BND/TRA)."""
    input:
        vcf_lr    = f"{RESULTS}/final/long_read_filtered.vcf",
        sv_counts = f"{RESULTS}/final/sv_counts.tsv"
    output: f"{RESULTS}/plots/sv_count_by_type.pdf"
    log:    f"{RESULTS}/logs/plots/sv_count_by_type.log"
    shell:
        "Rscript {SCRIPTS}/plot_sv_count_by_type.R {input.sv_counts} {output} \
            > {log} 2>&1"


rule plot_size_distribution:
    """Violin + boxplot: SV size distribution on log10 scale per type."""
    input:  f"{RESULTS}/final/long_read_filtered.vcf"
    output: f"{RESULTS}/plots/sv_size_distribution.pdf"
    log:    f"{RESULTS}/logs/plots/size_distribution.log"
    shell:
        "Rscript {SCRIPTS}/plot_size_distribution.R {input} {output} > {log} 2>&1"


rule plot_cxsv_venn:
    """
    Three-way Venn diagram: CxSV overlap between callers.
    Uses SURVIVOR SUPP_VEC from per-sample merged VCFs.
    SUPP_VEC bit order: CuteSV=bit0, SVIM=bit1, Sniffles2=bit2.
    Counts are pooled across all samples.
    """
    input:
        merged_vcfs = expand(
            f"{RESULTS}/per_sample/{{sample}}_merged_sv.vcf", sample=SAMPLES
        ),
        summary = f"{RESULTS}/final/cxsv_summary.tsv"
    output: f"{RESULTS}/plots/cxsv_venn.pdf"
    params:
        vcf_list = lambda wc, input: ",".join(input.merged_vcfs)
    log: f"{RESULTS}/logs/plots/cxsv_venn.log"
    shell:
        """
        Rscript {SCRIPTS}/plot_cxsv_venn.R \
            "{params.vcf_list}" \
            {input.summary} \
            {output} > {log} 2>&1
        """


rule plot_breakpoint_resolution:
    """
    Density plot: CIPOS interval width distribution per caller.
    Narrower CI = more precise breakpoint resolution.
    """
    input:
        sniffles = expand(
            f"{RESULTS}/callers/{{sample}}_sniffles_filtered.vcf", sample=SAMPLES
        ),
        cutesv = expand(
            f"{RESULTS}/callers/{{sample}}_cutesv_filtered.vcf", sample=SAMPLES
        ),
        svim = expand(
            f"{RESULTS}/callers/{{sample}}_svim_filtered.vcf", sample=SAMPLES
        )
    output: f"{RESULTS}/plots/breakpoint_resolution.pdf"
    params:
        sniffles_files = lambda wc, input: ",".join(input.sniffles),
        cutesv_files   = lambda wc, input: ",".join(input.cutesv),
        svim_files     = lambda wc, input: ",".join(input.svim)
    log: f"{RESULTS}/logs/plots/breakpoint_resolution.log"
    shell:
        """
        Rscript {SCRIPTS}/plot_breakpoint_resolution.R \
            {params.sniffles_files} \
            {params.cutesv_files} \
            {params.svim_files} \
            {output} > {log} 2>&1
        """


rule plot_genomic_context:
    """
    Bar plot: CxSV genomic context breakdown.
    Categories: Exonic / Promoter / Intronic /
                Repeat_intergenic / SegDup_intergenic / Intergenic.
    Annotation-only — no loci removed by context.
    Most CxSVs are expected in the grey-zone (intergenic/repeat/segdup).
    """
    input:  f"{RESULTS}/final/cxsv_master_table.tsv"
    output: f"{RESULTS}/plots/genomic_context.pdf"
    log:    f"{RESULTS}/logs/plots/genomic_context.log"
    shell:
        "Rscript {SCRIPTS}/plot_genomic_context.R {input} {output} > {log} 2>&1"


rule plot_complexity_class_heatmap:
    """
    Three-panel figure:
      Panel A — Complexity class frequency bar plot
      Panel B — SV type co-occurrence heatmap per complexity class
      Panel C — CxSV criteria frequency (C1-C5) bar plot
    """
    input:  f"{RESULTS}/final/cxsv_master_table.tsv"
    output: f"{RESULTS}/plots/complexity_class_heatmap.pdf"
    log:    f"{RESULTS}/logs/plots/complexity_class.log"
    shell:
        "Rscript {SCRIPTS}/plot_complexity_class.R {input} {output} > {log} 2>&1"


rule plot_chromosome_density:
    """
    Bar plot: SV density per chromosome (SVs per Mb).
    Normalized by CHM13v2.0 chromosome lengths.
    """
    input:  f"{RESULTS}/final/long_read_filtered.vcf"
    output: f"{RESULTS}/plots/chromosome_density.pdf"
    log:    f"{RESULTS}/logs/plots/chromosome_density.log"
    shell:
        "Rscript {SCRIPTS}/plot_chromosome_density.R {input} {output} > {log} 2>&1"
