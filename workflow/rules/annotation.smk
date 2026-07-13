# ═══════════════════════════════════════════════════════════════════════════════
# PHASE 5 — ANNOTATION
# ═══════════════════════════════════════════════════════════════════════════════

rule annotate_cxsv:
    """
    Build the 15-column CxSV master table and population BED file.
    Annotations added:
      - Gene overlap / nearest gene / distance to nearest gene
      - Exon overlap (yes/no)
      - Repeat element class (RepeatMasker)
      - Segmental duplication overlap fraction (0.0-1.0)
      - Genomic context label (Exonic/Intronic/Repeat_intergenic/etc.)
      - Complexity class (DEL+INV, BND_cluster, Chromothripsis_like, etc.)
    IMPORTANT: genomic context is annotation-only — no CxSVs are removed
    based on where they fall in the genome.
    """
    input:
        cxsv_summary = f"{RESULTS}/final/cxsv_summary.tsv",
        vcf          = f"{RESULTS}/final/cxsv_only.vcf",
        genes_bed    = config["annotation"]["genes_bed"],
        exons_bed    = config["annotation"]["exons_bed"],
        repeats_bed  = config["annotation"]["repeats_bed"],
        segdup_bed   = config["annotation"]["segdup_bed"]
    output:
        master_table = f"{RESULTS}/final/cxsv_master_table.tsv",
        pop_bed      = f"{RESULTS}/final/cxsv_population.bed"
    log: "logs/annotate_cxsv.log"
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
            --out-bed   {output.pop_bed} \
            > {log} 2>&1
        """
