# ═══════════════════════════════════════════════════════════════════════════════
# PRIVATE CxSV FILTERING — query cohorts only
#
# Removes CxSV loci already present in the pooled public reference master
# list (built by scripts/build_public_reference.py from every "reference"
# cohort), leaving only loci private to this cohort — candidate
# pathogenic/novel calls. Rules defined here are only pulled into `rule all`
# when the active cohort's role is "query" (see workflow/Snakefile).
#
# Matching uses a slop window (public_reference.match_window, default 1000 bp)
# rather than exact overlap: CxSV loci are already clustered regions, so this
# only needs to absorb boundary differences between independently-run
# cohorts, not merge genuinely distinct loci.
# ═══════════════════════════════════════════════════════════════════════════════

PUBLIC_REFERENCE_BED = config["public_reference"]["master_bed"]
MATCH_WINDOW         = config["public_reference"].get("match_window", 1000)


rule private_cxsv_population:
    """
    Split this cohort's population-level CxSV loci into:
      - private_cxsv_population.bed        — no reference locus within
                                              {MATCH_WINDOW} bp (candidate
                                              pathogenic/private)
      - public_overlap_cxsv_population.bed — a reference locus within
                                              {MATCH_WINDOW} bp (likely
                                              benign; kept for QC)
    Requires the reference master list to already exist —
    run scripts/build_public_reference.py first.
    """
    input:
        query_bed     = f"{RESULTS}/final/cxsv_population.bed",
        reference_bed = PUBLIC_REFERENCE_BED
    output:
        private_bed = f"{RESULTS}/final/private_cxsv_population.bed",
        shared_bed  = f"{RESULTS}/final/public_overlap_cxsv_population.bed"
    params:
        window = MATCH_WINDOW
    log: f"{RESULTS}/logs/private_filter/population.log"
    shell:
        """
        bedtools window -w {params.window} -v \
            -a {input.query_bed} -b {input.reference_bed} \
            > {output.private_bed} 2>{log}
        bedtools window -w {params.window} -u \
            -a {input.query_bed} -b {input.reference_bed} \
            > {output.shared_bed} 2>>{log}
        echo "private_cxsv_population: $(wc -l < {output.private_bed}) private, $(wc -l < {output.shared_bed}) shared with reference" \
            >> {log}
        """


rule private_cxsv_master_table:
    """
    Annotated master table filtered down to private-only CxSV loci (rows
    matching a locus in private_cxsv_population.bed). Exact match here —
    the slop tolerance was already applied upstream when that BED was built.
    """
    input:
        master_table = f"{RESULTS}/final/cxsv_master_table.tsv",
        private_bed  = f"{RESULTS}/final/private_cxsv_population.bed"
    output:
        f"{RESULTS}/final/private_cxsv_master_table.tsv"
    log: f"{RESULTS}/logs/private_filter/master_table.log"
    shell:
        """
        python {SCRIPTS}/filter_table_by_bed.py \
            --table  {input.master_table} \
            --bed    {input.private_bed} \
            --mode   keep \
            --window 0 \
            --output {output} \
            > {log} 2>&1
        """


rule extract_private_cxsv_per_sample:
    """
    Per-sample private (non-reference) CxSV VCF — the primary patient
    deliverable. Filters each sample's own cxsv_only.vcf directly against
    the reference master BED, using the same slop tolerance as
    private_cxsv_population above.
    """
    input:
        cxsv_vcf      = f"{RESULTS}/per_sample/{{sample}}_cxsv_only.vcf",
        reference_bed = PUBLIC_REFERENCE_BED
    output:
        f"{RESULTS}/per_sample/{{sample}}_private_cxsv.vcf"
    params:
        window = MATCH_WINDOW
    log: f"{RESULTS}/logs/private_filter/{{sample}}.log"
    shell:
        """
        python {SCRIPTS}/filter_vcf_by_bed.py \
            --vcf    {input.cxsv_vcf} \
            --bed    {input.reference_bed} \
            --mode   exclude \
            --window {params.window} \
            --output {output} \
            > {log} 2>&1
        """
