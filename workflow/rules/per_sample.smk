# ═══════════════════════════════════════════════════════════════════════════════
# PER-SAMPLE FILTERING, MERGING AND CxSV DISCOVERY  (shared — fingerprinted)
#
# CxSVs are discovered in each genome separately: a CxSV is a rearrangement in
# one person, so only that person's SVs are clustered and classified. The
# panel then counts carriers per locus (rules/panel_aggregate.smk); the patient
# workflow runs the very same rules on the patient (decision 2026-10-28, after
# pooled clustering across a cohort was found to chain SVs from different
# people into artificial loci covering ~12% of the genome at ~90 samples).
#
# Outputs per sample, under {S} = results/<workflow>/samples/<cohort>/<sample>:
#   merged_sv.vcf                  three-caller SURVIVOR merge
#   cxsv/cxsv_summary.tsv          every cluster, IsCxSV + criteria C1-C5
#   cxsv/cxsv_only.vcf             merged calls inside CxSV loci
#   cxsv/cxsv_master_table.tsv     annotated CxSV loci
#   cxsv/cxsv_loci.bed             Chr, Start, End, Exon_Disrupt
# ═══════════════════════════════════════════════════════════════════════════════

rule filter_sniffles:
    """Sniffles2: keep PASS calls with QUAL >= filter_params.min_qual."""
    input:  f"{S}/callers/sniffles.vcf"
    output: f"{S}/callers/sniffles_filtered.vcf"
    params: min_qual = config["filter_params"]["min_qual"]
    log:    f"{S}/logs/filter_sniffles.log"
    shell:
        "bcftools view -f PASS -e 'QUAL < {params.min_qual}' {input} -o {output} -O v > {log} 2>&1"


rule filter_cutesv:
    """cuteSV: keep PASS calls."""
    input:  f"{S}/callers/cutesv.vcf"
    output: f"{S}/callers/cutesv_filtered.vcf"
    log:    f"{S}/logs/filter_cutesv.log"
    shell:  "bcftools view -f PASS {input} -o {output} -O v > {log} 2>&1"


rule filter_svim:
    """SVIM: keep calls with QUAL >= filter_params.svim_min_qual (SVIM scores in QUAL)."""
    input:  f"{S}/callers/svim/variants.vcf"
    output: f"{S}/callers/svim_filtered.vcf"
    params: min_qual = config["filter_params"]["svim_min_qual"]
    log:    f"{S}/logs/filter_svim.log"
    shell:
        "bcftools view -e 'QUAL < {params.min_qual}' {input} -o {output} -O v > {log} 2>&1"


rule merge_sv_callers:
    """
    SURVIVOR merge of the three filtered call sets.
    SUPP_VEC bit order: CuteSV=bit0, SVIM=bit1, Sniffles2=bit2.
    """
    input:
        cutesv   = f"{S}/callers/cutesv_filtered.vcf",
        svim     = f"{S}/callers/svim_filtered.vcf",
        sniffles = f"{S}/callers/sniffles_filtered.vcf"
    output:
        vcf      = f"{S}/merged_sv.vcf",
        filelist = f"{S}/merged_vcf_list.txt"
    params:
        max_dist    = config["survivor_params"]["max_dist"],
        min_callers = config["survivor_params"]["min_callers"],
        min_sv_len  = config["sv_params"]["min_sv_len"]
    log: f"{S}/logs/merge_sv_callers.log"
    shell:
        """
        echo {input.cutesv}   >  {output.filelist}
        echo {input.svim}     >> {output.filelist}
        echo {input.sniffles} >> {output.filelist}
        SURVIVOR merge {output.filelist} {params.max_dist} \
            {params.min_callers} 1 0 0 {params.min_sv_len} {output.vcf} > {log} 2>&1
        """


rule sample_vcf_to_bed:
    """Merged VCF -> BED for clustering; each BND/TRA emits both ends."""
    input:
        vcf    = f"{S}/merged_sv.vcf",
        script = f"{SCRIPTS}/vcf_to_bed.py"
    output: f"{S}/cxsv/merged_sv.bed"
    params: window = config["cxsv_params"]["bnd_window"]
    log:    f"{S}/logs/vcf_to_bed.log"
    shell:
        "python {input.script} --vcf {input.vcf} --out {output} --window {params.window} 2> {log}"


rule sample_cluster_sv:
    """Cluster this sample's SVs within cxsv_params.cluster_distance."""
    input:  f"{S}/cxsv/merged_sv.bed"
    output:
        sorted_bed   = f"{S}/cxsv/merged_sv_sorted.bed",
        clusters_bed = f"{S}/cxsv/clusters.bed"
    params: cluster_dist = config["cxsv_params"]["cluster_distance"]
    log:    f"{S}/logs/cluster_sv.log"
    shell:
        """
        bedtools sort -i {input} > {output.sorted_bed} 2>> {log}
        bedtools cluster -d {params.cluster_dist} -i {output.sorted_bed} \
            > {output.clusters_bed} 2>> {log}
        """


rule sample_classify_cxsv:
    """
    Classify this sample's clusters against criteria C1-C5
    (workflow/scripts/classify_cxsv.py). Genomic context is never a criterion.
    """
    input:  f"{S}/cxsv/clusters.bed"
    output: f"{S}/cxsv/cxsv_summary.tsv"
    params:
        min_bp      = config["cxsv_params"]["min_breakpoints"],
        max_win     = config["cxsv_params"]["cluster_distance"],
        ca_min      = config["cxsv_params"]["chromoanag_min_window"],
        ca_max      = config["cxsv_params"]["chromoanag_max_window"],
        ca_min_sv   = config["cxsv_params"]["chromoanag_min_svs"],
        max_cluster = config["cxsv_params"]["max_cluster_size"]
    log: f"{S}/logs/classify_cxsv.log"
    shell:
        """
        python {SCRIPTS}/classify_cxsv.py \
            --bed               {input} \
            --output            {output} \
            --min-bp            {params.min_bp} \
            --max-window        {params.max_win} \
            --chromoanag-min    {params.ca_min} \
            --chromoanag-max    {params.ca_max} \
            --chromoanag-min-sv {params.ca_min_sv} \
            --max-cluster-size  {params.max_cluster} \
            2> {log}
        """


rule sample_extract_cxsv_vcf:
    """This sample's merged calls that fall inside its own CxSV loci."""
    input:
        vcf     = f"{S}/merged_sv.vcf",
        summary = f"{S}/cxsv/cxsv_summary.tsv"
    output: f"{S}/cxsv/cxsv_only.vcf"
    log:    f"{S}/logs/extract_cxsv_vcf.log"
    shell:
        "python {SCRIPTS}/extract_cxsv_vcf.py --vcf {input.vcf} --summary {input.summary} "
        "--output {output} > {log} 2>&1"


rule sample_annotate_cxsv:
    """
    Annotated CxSV master table and loci BED for this sample (genes, exons,
    repeats, segdups, complexity class). Annotation never removes a locus.
    """
    input:
        cxsv_summary = f"{S}/cxsv/cxsv_summary.tsv",
        vcf          = f"{S}/cxsv/cxsv_only.vcf",
        genes_bed    = config["annotation"]["genes_bed"],
        exons_bed    = config["annotation"]["exons_bed"],
        repeats_bed  = config["annotation"]["repeats_bed"],
        segdup_bed   = config["annotation"]["segdup_bed"]
    output:
        master_table = f"{S}/cxsv/cxsv_master_table.tsv",
        loci_bed     = f"{S}/cxsv/cxsv_loci.bed"
    log: f"{S}/logs/annotate_cxsv.log"
    shell:
        """
        python {SCRIPTS}/annotate_cxsv.py \
            --summary   {input.cxsv_summary} \
            --vcf       {input.vcf} \
            --genes     {input.genes_bed} \
            --exons     {input.exons_bed} \
            --repeats   {input.repeats_bed} \
            --segdups   {input.segdup_bed} \
            --out-table {output.master_table} \
            --out-bed   {output.loci_bed} \
            > {log} 2>&1
        """
