# ═══════════════════════════════════════════════════════════════════════════════
# PHASE 1 — BAM INDEXING
# Creates the .bam.bai index required by all three SV callers.
# samtools index takes ~5 min per BAM and produces a ~49 MB .bai file.
# ═══════════════════════════════════════════════════════════════════════════════

rule index_bam:
    """
    Index each BAM with samtools.
    Input:  {data_dir}/{sample}.bam
    Output: {data_dir}/{sample}.bam.bai
    All three callers depend on this rule before they can start.
    """
    input:
        bam = f"{DATA}/{{sample}}.bam"
    output:
        bai = f"{DATA}/{{sample}}.bam.bai"
    threads: config["threads"]["samtools"]
    log:
        f"{RESULTS}/logs/index_bam/{{sample}}.log"
    shell:
        "samtools index -@ {threads} {input.bam} > {log} 2>&1"
