# ═══════════════════════════════════════════════════════════════════════════════
# PHASE 4 — CxSV DISCOVERY
# ═══════════════════════════════════════════════════════════════════════════════

rule vcf_to_bed:
    """
    Convert the population VCF to BED format for clustering.
    Uses vcf_to_bed.py which handles BND/TRA records correctly —
    each BND emits two BED records (origin + mate position) so
    inter-chromosomal rearrangements are captured in clustering.
    """
    input:
        vcf    = f"{RESULTS}/final/long_read_filtered.vcf",
        script = f"{SCRIPTS}/vcf_to_bed.py"
    output:
        bed = f"{RESULTS}/final/merged_sv.bed"
    params:
        window = config["cxsv_params"]["bnd_window"]
    log: f"{RESULTS}/logs/vcf_to_bed.log"
    shell:
        """
        python {input.script} \
            --vcf    {input.vcf} \
            --out    {output.bed} \
            --window {params.window} \
            2>{log}
        echo "vcf_to_bed: $(wc -l < {output.bed}) records written" >> {log}
        """


rule cluster_sv:
    """
    Sort the BED and cluster SVs within the 50 kb window.
    Produces clusters.bed with a cluster ID appended as the last column.
    Split into two steps so each command's stderr is captured separately.
    """
    input:  f"{RESULTS}/final/merged_sv.bed"
    output:
        sorted_bed   = f"{RESULTS}/final/merged_sv_sorted.bed",
        clusters_bed = f"{RESULTS}/final/clusters.bed"
    params:
        cluster_dist = config["cxsv_params"]["cluster_distance"]
    log: f"{RESULTS}/logs/cluster_sv.log"
    shell:
        """
        bedtools sort -i {input} \
            > {output.sorted_bed} 2>>{log}
        bedtools cluster \
            -d {params.cluster_dist} \
            -i {output.sorted_bed} \
            > {output.clusters_bed} 2>>{log}
        echo "cluster_sv: $(wc -l < {output.clusters_bed}) records" >> {log}
        """


rule classify_cxsv:
    """
    Classify SV clusters as CxSV using all five project criteria.
    See workflow/scripts/classify_cxsv.py for full documentation.

    C1 — >= 3 unique breakpoints within 50 kb window
    C2 — Nested SV (one interval wholly inside another)
    C3 — Overlapping SV types (physically share coordinates)
    C4 — Copy-number shift (DEL/DUP) + orientation change (INV/BND/TRA)
    C5 — Chromoanagenesis cluster: >= 3 SVs in 10-50 kb window

    Outputs a CxSV_Criteria column listing which criteria fired.
    CxSVs are NEVER filtered by genomic context — intergenic, repeat,
    and segdup loci are retained (most CxSVs are expected there).
    """
    input:  f"{RESULTS}/final/clusters.bed"
    output: f"{RESULTS}/final/cxsv_summary.tsv"
    log:    f"{RESULTS}/logs/classify_cxsv.log"
    params:
        min_bp      = config["cxsv_params"]["min_breakpoints"],
        max_win     = config["cxsv_params"]["cluster_distance"],
        ca_min      = config["cxsv_params"]["chromoanag_min_window"],
        ca_max      = config["cxsv_params"]["chromoanag_max_window"],
        ca_min_sv   = config["cxsv_params"]["chromoanag_min_svs"],
        max_cluster = config["cxsv_params"]["max_cluster_size"]
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
            2>{log}
        echo "classify_cxsv: $(grep -c 'YES' {output} || true) CxSV loci" \
            >> {log}
        """


rule count_cxsv:
    """
    Count CxSV loci and break down by:
      - Inter vs intra-chromosomal
      - Which of C1-C5 criteria fired (loci can match multiple)
    Column layout from classify_cxsv.py:
      col 9  = Has_Inter
      col 10 = IsCxSV
      col 11 = CxSV_Criteria
    """
    input:  f"{RESULTS}/final/cxsv_summary.tsv"
    output: f"{RESULTS}/final/cxsv_count.txt"
    log:    f"{RESULTS}/logs/count_cxsv.log"
    shell:
        r"""
        awk 'BEGIN{{OFS="\t"}}
             NR==1{{next}}
             $10=="YES"{{
                 total++;
                 if ($9=="YES") inter++; else intra++;
                 n=split($11,crit,",");
                 for (i=1;i<=n;i++) crit_count[crit[i]]++;
             }}
             END{{
                 print "Total CxSV loci:", total+0;
                 print "  Inter-chromosomal:", inter+0;
                 print "  Intra-chromosomal only:", intra+0;
                 print "";
                 print "Criteria breakdown (loci can match multiple):";
                 for (c in crit_count) print "  "c":", crit_count[c];
             }}' {input} > {output} 2>{log}
        """


rule extract_cxsv_vcf:
    """Extract the population-level CxSV-only VCF."""
    input:
        vcf     = f"{RESULTS}/final/long_read_filtered.vcf",
        summary = f"{RESULTS}/final/cxsv_summary.tsv"
    output: f"{RESULTS}/final/cxsv_only.vcf"
    log: f"{RESULTS}/logs/extract_cxsv_vcf.log"
    shell:
        """
        python {SCRIPTS}/extract_cxsv_vcf.py \
            --vcf     {input.vcf} \
            --summary {input.summary} \
            --output  {output} \
            > {log} 2>&1
        """


rule extract_per_sample_cxsv:
    """Extract a CxSV-only VCF for each individual sample."""
    input:
        vcf     = f"{RESULTS}/per_sample/{{sample}}_merged_sv.vcf",
        summary = f"{RESULTS}/final/cxsv_summary.tsv"
    output: f"{RESULTS}/per_sample/{{sample}}_cxsv_only.vcf"
    log: f"{RESULTS}/logs/cxsv_vcf/{{sample}}.log"
    shell:
        """
        python {SCRIPTS}/extract_cxsv_vcf.py \
            --vcf     {input.vcf} \
            --summary {input.summary} \
            --output  {output} \
            > {log} 2>&1
        """
