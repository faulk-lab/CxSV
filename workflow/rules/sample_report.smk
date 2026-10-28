# ═══════════════════════════════════════════════════════════════════════════════
# PER-SAMPLE REPORT  (shared, NOT fingerprinted: reporting only)
# ═══════════════════════════════════════════════════════════════════════════════

rule sample_sv_summary:
    """
    Per-sample text report: SV counts by type per caller, CxSV counts by
    criterion and complexity class, and caller concordance (SURVIVOR SUPP_VEC;
    bit order CuteSV=bit0, SVIM=bit1, Sniffles2=bit2).
    """
    input:
        merged   = f"{S}/merged_sv.vcf",
        cutesv   = f"{S}/callers/cutesv_filtered.vcf",
        sniffles = f"{S}/callers/sniffles_filtered.vcf",
        svim     = f"{S}/callers/svim_filtered.vcf",
        cxsv_vcf = f"{S}/cxsv/cxsv_only.vcf",
        summary  = f"{S}/cxsv/cxsv_summary.tsv"
    output: f"{S}/sv_summary.txt"
    log:    f"{S}/logs/sv_summary.log"
    shell:
        """
        python {SCRIPTS}/summarize_sv.py \
            --merged   {input.merged} \
            --cutesv   {input.cutesv} \
            --sniffles {input.sniffles} \
            --svim     {input.svim} \
            --cxsv     {input.cxsv_vcf} \
            --summary  {input.summary} \
            --sample   {wildcards.sample} \
            --output   {output} \
            > {log} 2>&1
        """
