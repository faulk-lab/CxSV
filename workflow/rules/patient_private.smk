# ═══════════════════════════════════════════════════════════════════════════════
# PATIENT — PANEL COMPATIBILITY, PRIVATE/PUBLIC LABELLING, PER-SAMPLE PLOTS
#
# Needs only a panel release directory (config: panel_release), never the
# background reads, so it can run inside a restricted environment.
#
# 1. check_panel_compatibility — this run's calling fingerprint must equal the
#    panel's (same parameters, reference, tool versions and loci code). A
#    mismatch stops the run: otherwise differences in how the two sides were
#    called would surface as false "private" CxSVs. Override only for
#    exploratory work with allow_fingerprint_mismatch: true (the mismatch is
#    then written into panel_compatibility.txt and every output carries it).
# 2. classify_patient_cxsv — carrier count, frequency and public/private label
#    for every patient CxSV locus.
# ═══════════════════════════════════════════════════════════════════════════════

PANEL_RELEASE = config.get("panel_release")
if not PANEL_RELEASE:
    raise ValueError("[patient] config['panel_release'] must point at a panel release directory")
MATCH = config.get("panel_match", {})


rule check_panel_compatibility:
    """
    Compare this run's calling fingerprint with the panel release's. The full
    report goes to logs/panel_compatibility.log, which survives a failure
    (Snakemake deletes a failed job's outputs, not its log); the output file
    is a copy written only when the run may proceed.
    """
    input:
        ours   = f"{RESULTS_ROOT}/calling_fingerprint.json",
        theirs = f"{PANEL_RELEASE}/calling_fingerprint.json"
    output: f"{RESULTS_ROOT}/panel_compatibility.txt"
    log:    f"{RESULTS_ROOT}/logs/panel_compatibility.log"
    run:
        with open(input.ours) as fh:
            ours = json.load(fh)
        with open(input.theirs) as fh:
            theirs = json.load(fh)
        diffs = cxsv_fingerprint.compare(ours, theirs)
        report = [f"this run fingerprint:  {ours['sha256']}",
                  f"panel fingerprint:     {theirs.get('sha256')}",
                  f"panel release:         {PANEL_RELEASE}"]
        if not diffs:
            report.append("status: COMPATIBLE")
        else:
            report.append(f"status: MISMATCH ({len(diffs)} difference(s))")
            report += [f"  {d}" for d in diffs]
            if config.get("allow_fingerprint_mismatch", False):
                report.append("allow_fingerprint_mismatch is set: continuing for exploratory use; "
                              "public/private labels are NOT comparable with the panel")
        text = "\n".join(report) + "\n"
        with open(log[0], "w") as fh:
            fh.write(text)
        if diffs and not config.get("allow_fingerprint_mismatch", False):
            raise ValueError(
                f"this run is not comparable with the panel at {PANEL_RELEASE}: "
                f"{len(diffs)} difference(s), first: {diffs[0]}. Full list: {log[0]}. "
                "Rebuild the panel with this configuration, or set "
                "allow_fingerprint_mismatch: true for exploratory use only.")
        with open(output[0], "w") as fh:
            fh.write(text)


rule classify_patient_cxsv:
    input:
        table   = f"{S}/cxsv/cxsv_master_table.tsv",
        calls   = f"{PANEL_RELEASE}/panel_cxsv_calls.tsv",
        summary = f"{PANEL_RELEASE}/panel_summary.tsv",
        compat  = f"{RESULTS_ROOT}/panel_compatibility.txt",
        script  = f"{SCRIPTS}/classify_patient_cxsv.py"
    output:
        all_loci = f"{S}/private/cxsv_panel_status.tsv",
        private  = f"{S}/private/private_cxsv_master_table.tsv",
        bed      = f"{S}/private/private_cxsv_loci.bed"
    params:
        window       = MATCH.get("window", 1000),
        min_carriers = MATCH.get("min_carriers", 1)
    log: f"{S}/logs/classify_patient_cxsv.log"
    shell:
        """
        python {input.script} --master-table {input.table} --panel-calls {input.calls} \
            --panel-summary {input.summary} --window {params.window} \
            --min-carriers {params.min_carriers} --out-all {output.all_loci} \
            --out-private {output.private} --out-bed {output.bed} 2> {log}
        """


rule patient_private_vcf:
    """The patient's merged SV records inside private CxSV loci."""
    input:
        vcf = f"{S}/cxsv/cxsv_only.vcf",
        bed = f"{S}/private/private_cxsv_loci.bed"
    output: f"{S}/private/private_cxsv.vcf"
    log:    f"{S}/logs/private_cxsv_vcf.log"
    shell:
        "python {SCRIPTS}/filter_vcf_by_bed.py --vcf {input.vcf} --bed {input.bed} "
        "--mode keep --window 0 --output {output} > {log} 2>&1"


rule patient_plot_genomic_context:
    input:  f"{S}/cxsv/cxsv_master_table.tsv"
    output: f"{S}/plots/genomic_context.pdf"
    log:    f"{S}/logs/plots/genomic_context.log"
    shell:  "Rscript {SCRIPTS}/plot_genomic_context.R {input} {output} > {log} 2>&1"


rule patient_plot_complexity_class:
    input:  f"{S}/cxsv/cxsv_master_table.tsv"
    output: f"{S}/plots/complexity_class_heatmap.pdf"
    log:    f"{S}/logs/plots/complexity_class.log"
    shell:  "Rscript {SCRIPTS}/plot_complexity_class.R {input} {output} > {log} 2>&1"
