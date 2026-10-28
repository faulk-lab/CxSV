# ═══════════════════════════════════════════════════════════════════════════════
# PER-SAMPLE JUNCTION TABLE  (shared — fingerprinted)
#
# Every merged SV call of a sample as a breakpoint junction (chrom1 pos1 chrom2
# pos2 svtype svlen), with duplicate and reciprocal breakend records collapsed
# and each junction's CxSV locus, if any (workflow/lib/cxsv_junctions.py,
# docs/JUNCTIONS.md). The panel clusters these tables into junction sites; the
# patient workflow matches its CxSV loci against those sites junction by
# junction (roadmap item 3). Both sides must build the table the same way, so
# this file, the script, the library and junction_params are fingerprinted.
#
#   {S}/junctions/junctions.tsv
# ═══════════════════════════════════════════════════════════════════════════════

rule sample_junctions:
    input:
        vcf      = f"{S}/merged_sv.vcf",
        clusters = f"{S}/cxsv/clusters.bed",
        summary  = f"{S}/cxsv/cxsv_summary.tsv",
        fai      = f"{REFERENCE}.fai",
        script   = f"{SCRIPTS}/extract_junctions.py"
    output: f"{S}/junctions/junctions.tsv"
    params:
        distance  = config["junction_params"]["collapse_distance"],
        ins_ratio = config["junction_params"]["ins_len_ratio"]
    log: f"{S}/logs/extract_junctions.log"
    shell:
        "python {input.script} --vcf {input.vcf} --clusters {input.clusters} "
        "--summary {input.summary} --fai {input.fai} --collapse-distance {params.distance} "
        "--ins-len-ratio {params.ins_ratio} --out {output} 2> {log}"
