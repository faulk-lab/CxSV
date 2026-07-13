# ═══════════════════════════════════════════════════════════════════════════════
# PHASE 6 — SUMMARY STATS
# ═══════════════════════════════════════════════════════════════════════════════

rule calculate_sv_counts:
    """Tabulate SV counts by type (DEL/DUP/INV/INS/BND/TRA) across all samples."""
    input:  f"{RESULTS}/final/long_read_filtered.vcf"
    output: f"{RESULTS}/final/sv_counts.tsv"
    log:    "logs/sv_counts.log"
    shell:
        r"""
        bcftools view {input} | grep -v "^#" \
        | awk -F'\t' '{{
            n = split($8, info, ";");
            for (i = 1; i <= n; i++) {{
                if (info[i] ~ /^SVTYPE=/) {{
                    split(info[i], kv, "=");
                    print kv[2];
                    break
                }}
            }}
        }}' \
        | sort | uniq -c \
        | awk 'BEGIN{{print "Count\tSVTYPE"}} {{print $1"\t"$2}}' \
        > {output} 2>{log}
        """


rule summarize_sv_overlap:
    """
    Per-sample summary report with three sections:
      Section 1 — SV counts by type per caller (Sniffles2 / CuteSV / SVIM / Merged)
      Section 2 — CxSV count + breakdown by criteria (C1-C5) + complexity class
      Section 3 — Caller concordance from SURVIVOR SUPP_VEC
                  (bit order: CuteSV=bit0, SVIM=bit1, Sniffles2=bit2)
    """
    input:
        merged   = f"{RESULTS}/per_sample/{{sample}}_merged_sv.vcf",
        cutesv   = f"{RESULTS}/callers/{{sample}}_cutesv_filtered.vcf",
        sniffles = f"{RESULTS}/callers/{{sample}}_sniffles_filtered.vcf",
        svim     = f"{RESULTS}/callers/{{sample}}_svim_filtered.vcf",
        cxsv_vcf = f"{RESULTS}/per_sample/{{sample}}_cxsv_only.vcf",
        summary  = f"{RESULTS}/final/cxsv_summary.tsv"
    output: f"{RESULTS}/per_sample/{{sample}}_sv_summary.txt"
    log:    f"logs/summary/{{sample}}.log"
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
