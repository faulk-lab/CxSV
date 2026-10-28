# ═══════════════════════════════════════════════════════════════════════════════
# PANEL — PER-COHORT POPULATION CALL SETS AND QC PLOTS
#
# SV-level summaries of each background cohort. These are QC and description
# only: CxSV loci come from each sample separately (rules/per_sample.smk), not
# from the cohort merge, which is NOT used for discovery.
# Outputs under results/panel/cohorts/<cohort>/.
# ═══════════════════════════════════════════════════════════════════════════════

C = f"{RESULTS_ROOT}/cohorts/{{cohort}}"


def cohort_files(pattern):
    """Input function: one path per sample of the wildcard cohort."""
    return lambda wc: per_sample(pattern, wc.cohort)


def joined(name):
    """params helper: comma-join a named input list (for the R scripts)."""
    return lambda wc, input: ",".join(input[name])


rule cohort_population_vcf:
    """SURVIVOR merge of the cohort's per-sample call sets (>= 1 caller)."""
    input:  vcfs = cohort_files(f"{S}/merged_sv.vcf")
    output:
        vcf      = f"{C}/long_read_filtered.vcf",
        filelist = f"{C}/population_vcf_list.txt"
    params:
        max_dist   = config["survivor_params"]["max_dist"],
        min_sv_len = config["sv_params"]["min_sv_len"]
    log: f"{C}/logs/population_vcf.log"
    run:
        with open(output.filelist, "w") as fh:
            fh.write("\n".join(input.vcfs) + "\n")
        shell("SURVIVOR merge {output.filelist} {params.max_dist} 1 1 0 0 "
              "{params.min_sv_len} {output.vcf} > {log} 2>&1")


rule cohort_consensus_vcf:
    """Cohort SVs supported by >= survivor_params.min_callers_consensus callers."""
    input:  vcfs = cohort_files(f"{S}/merged_sv.vcf")
    output:
        vcf      = f"{C}/consensus.vcf",
        filelist = f"{C}/consensus_vcf_list.txt"
    params:
        max_dist    = config["survivor_params"]["max_dist"],
        min_callers = config["survivor_params"]["min_callers_consensus"],
        min_sv_len  = config["sv_params"]["min_sv_len"]
    log: f"{C}/logs/consensus_vcf.log"
    run:
        with open(output.filelist, "w") as fh:
            fh.write("\n".join(input.vcfs) + "\n")
        shell("SURVIVOR merge {output.filelist} {params.max_dist} {params.min_callers} 1 0 0 "
              "{params.min_sv_len} {output.vcf} > {log} 2>&1")


rule cohort_sv_counts:
    """SV counts by type in the cohort population VCF."""
    input:  f"{C}/long_read_filtered.vcf"
    output: f"{C}/sv_counts.tsv"
    log:    f"{C}/logs/sv_counts.log"
    shell:
        r"""
        bcftools query -f '%INFO/SVTYPE\n' {input} 2> {log} \
        | sort | uniq -c \
        | awk 'BEGIN{{print "Count\tSVTYPE"}} {{print $1"\t"$2}}' > {output}
        """


rule cohort_cxsv_tables:
    """Stack the cohort's per-sample CxSV master tables and summaries (Sample_ID first)."""
    input:
        tables    = cohort_files(f"{S}/cxsv/cxsv_master_table.tsv"),
        summaries = cohort_files(f"{S}/cxsv/cxsv_summary.tsv")
    output:
        table   = f"{C}/cxsv_master_table.tsv",
        summary = f"{C}/cxsv_summary.tsv"
    run:
        for paths, out, keep in ((input.tables, output.table, None),
                                 (input.summaries, output.summary, "IsCxSV")):
            header_written = False
            with open(out, "w") as fo:
                for path in paths:
                    sample = os.path.basename(os.path.dirname(os.path.dirname(path)))
                    with open(path) as fi:
                        header = fi.readline().rstrip("\n").split("\t")
                        if not header_written:
                            fo.write("Sample_ID\t" + "\t".join(header) + "\n")
                            header_written = True
                        k = header.index(keep) if keep else None
                        for line in fi:
                            f = line.rstrip("\n").split("\t")
                            if k is None or f[k] == "YES":
                                fo.write(sample + "\t" + "\t".join(f) + "\n")


rule plot_sv_count_by_type:
    input:  sv_counts = f"{C}/sv_counts.tsv"
    output: f"{C}/plots/sv_count_by_type.pdf"
    log:    f"{C}/logs/plots/sv_count_by_type.log"
    shell:  "Rscript {SCRIPTS}/plot_sv_count_by_type.R {input.sv_counts} {output} > {log} 2>&1"


rule plot_size_distribution:
    input:  f"{C}/long_read_filtered.vcf"
    output: f"{C}/plots/sv_size_distribution.pdf"
    log:    f"{C}/logs/plots/size_distribution.log"
    shell:  "Rscript {SCRIPTS}/plot_size_distribution.R {input} {output} > {log} 2>&1"


rule plot_chromosome_density:
    input:  f"{C}/long_read_filtered.vcf"
    output: f"{C}/plots/chromosome_density.pdf"
    log:    f"{C}/logs/plots/chromosome_density.log"
    shell:  "Rscript {SCRIPTS}/plot_chromosome_density.R {input} {output} > {log} 2>&1"


rule plot_breakpoint_resolution:
    """CIPOS interval width per caller (narrower = more precise breakpoints)."""
    input:
        sniffles = cohort_files(f"{S}/callers/sniffles_filtered.vcf"),
        cutesv   = cohort_files(f"{S}/callers/cutesv_filtered.vcf"),
        svim     = cohort_files(f"{S}/callers/svim_filtered.vcf")
    output: f"{C}/plots/breakpoint_resolution.pdf"
    params:
        sniffles = joined("sniffles"),
        cutesv   = joined("cutesv"),
        svim     = joined("svim")
    log: f"{C}/logs/plots/breakpoint_resolution.log"
    shell:
        "Rscript {SCRIPTS}/plot_breakpoint_resolution.R {params.sniffles} {params.cutesv} "
        "{params.svim} {output} > {log} 2>&1"


rule plot_cxsv_venn:
    """
    Caller concordance for calls inside CxSV loci, pooled over the cohort's
    samples. Each sample contributes only its own CxSV calls (cxsv_only.vcf).
    """
    input:
        vcfs    = cohort_files(f"{S}/cxsv/cxsv_only.vcf"),
        summary = f"{C}/cxsv_summary.tsv"
    output: f"{C}/plots/cxsv_venn.pdf"
    params: vcfs = joined("vcfs")
    log:    f"{C}/logs/plots/cxsv_venn.log"
    shell:
        'Rscript {SCRIPTS}/plot_cxsv_venn.R "{params.vcfs}" {input.summary} {output} > {log} 2>&1'


rule plot_genomic_context:
    """CxSV genomic context (annotation only; nothing is filtered by context)."""
    input:  f"{C}/cxsv_master_table.tsv"
    output: f"{C}/plots/genomic_context.pdf"
    log:    f"{C}/logs/plots/genomic_context.log"
    shell:  "Rscript {SCRIPTS}/plot_genomic_context.R {input} {output} > {log} 2>&1"


rule plot_complexity_class_heatmap:
    """Complexity class frequency, SV-type co-occurrence and C1-C5 criteria."""
    input:  f"{C}/cxsv_master_table.tsv"
    output: f"{C}/plots/complexity_class_heatmap.pdf"
    log:    f"{C}/logs/plots/complexity_class.log"
    shell:  "Rscript {SCRIPTS}/plot_complexity_class.R {input} {output} > {log} 2>&1"


COHORT_OUTPUTS = [
    "long_read_filtered.vcf", "consensus.vcf", "sv_counts.tsv", "cxsv_master_table.tsv",
    "plots/sv_count_by_type.pdf", "plots/sv_size_distribution.pdf",
    "plots/chromosome_density.pdf", "plots/breakpoint_resolution.pdf",
    "plots/cxsv_venn.pdf", "plots/genomic_context.pdf", "plots/complexity_class_heatmap.pdf",
]
