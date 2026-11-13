# ═══════════════════════════════════════════════════════════════════════════════
# PATIENT — CLINICIAN REPORT  (roadmap item 16l, panel-v1.3; NOT fingerprinted: reporting only)
#
# The 2-page clinician-facing research report of each patient sample, built from this
# run's outputs only (scripts/clinician_report.py; fixed prose in config/clinician_report.yaml):
#
#   report/<sample>.clinician_report.docx   editable report
#   report/<sample>.clinician_report.pdf    the same report as PDF
#   report/<sample>.tier1_points.png        Tier 1 points graphic embedded in both
#                                           (a one-line note when the sample has no Tier 1)
#   report/<sample>.reconstruction/         one structural-reconstruction panel per Tier 1
#                                           locus (NN_<locus>.png; empty without Tier 1)
#   report/<sample>.tier1_reconstruction.tsv  segments, junctions with support, derivative,
#                                           alternatives and ambiguity of each Tier 1 locus
#
# The reconstruction (lib/cxsv_recon.py, drawn by lib/cxsv_recon_figure.py) uses only this
# run's junction table, filtered caller VCFs, merged CxSV VCF, locus read depth, CN segments
# and the genes BED: no alignment access, so it runs offline in the bundle.
#
# Needs python-docx and reportlab (envs/cxsv-pipeline.pip.txt; container image from panel-v1.3).
# Bundle identity comes from CXSV_BUNDLE_ID / CXSV_COMMIT (set by container/run_cxsv.sh).
# ═══════════════════════════════════════════════════════════════════════════════

import json as _json

REPORT_TEXT = config.get("clinician_report") or {}


rule clinician_report_text:
    """The report prose this run uses (config/clinician_report.yaml after site overlays)."""
    output: f"{RESULTS_ROOT}/clinician_report_text.json"
    run:
        if not REPORT_TEXT:
            raise ValueError("[clinician_report] config has no clinician_report section "
                             "(config/clinician_report.yaml is loaded by workflow/patient/Snakefile)")
        with open(output[0], "w") as fh:
            _json.dump({"clinician_report": REPORT_TEXT}, fh, indent=2, sort_keys=True)


rule clinician_report:
    """Clinician-facing research report (docx + pdf) with the Tier 1 points graphic and the
    structural reconstruction of each Tier 1 locus."""
    input:
        scores  = f"{S}/scoring/cxsv_scores.tsv",
        quarantine = f"{S}/scoring/quarantined_loci.tsv",
        params  = f"{S}/scoring/scoring_params.json",
        stats   = f"{S}/cn/stats.json",
        cn_only = f"{S}/cn/cn_only.tsv",
        libqc   = f"{S}/qc/library_qc.tsv",
        depth   = f"{S}/cn/depth.regions.bed.gz",      # cn_depth.sh also writes depth.mosdepth.summary.txt
        compat  = f"{RESULTS_ROOT}/panel_compatibility.txt",
        text    = f"{RESULTS_ROOT}/clinician_report_text.json",
        junctions = f"{S}/junctions/junctions.tsv",
        cxsv_vcf  = f"{S}/cxsv/cxsv_only.vcf",
        sniffles  = f"{S}/callers/sniffles_filtered.vcf",
        cutesv    = f"{S}/callers/cutesv_filtered.vcf",
        svim      = f"{S}/callers/svim_filtered.vcf",
        locus_depth = f"{S}/scoring/locus_depth.tsv",
        cn_segments = f"{S}/cn/segments.tsv",
        genes   = config["annotation"]["genes_bed"],
        script  = f"{SCRIPTS}/clinician_report.py",
        figlib  = f"{SCRIPTS}/../lib/cxsv_report_figure.py",
        reconlib = f"{SCRIPTS}/../lib/cxsv_recon.py",
        recondraw = f"{SCRIPTS}/../lib/cxsv_recon_figure.py"
    output:
        docx   = f"{S}/report/{{sample}}.clinician_report.docx",
        pdf    = f"{S}/report/{{sample}}.clinician_report.pdf",
        figure = f"{S}/report/{{sample}}.tier1_points.png",
        recon_dir = directory(f"{S}/report/{{sample}}.reconstruction"),
        recon_tsv = f"{S}/report/{{sample}}.tier1_reconstruction.tsv"
    params:
        summary = f"{S}/cn/depth.mosdepth.summary.txt"
    log: f"{S}/logs/clinician_report.log"
    shell:
        "MPLCONFIGDIR=${{TMPDIR:-/tmp}}/mpl_$$ python {input.script} --sample {wildcards.sample} "
        "--scores {input.scores} --quarantine {input.quarantine} --params {input.params} --cn-stats {input.stats} --cn-only {input.cn_only} "
        "--library-qc {input.libqc} --depth-summary {params.summary} --compatibility {input.compat} "
        "--text {input.text} --out-docx {output.docx} --out-pdf {output.pdf} --out-figure {output.figure} "
        "--junctions {input.junctions} --cxsv-vcf {input.cxsv_vcf} --caller-vcf Sniffles2={input.sniffles} "
        "--caller-vcf cuteSV={input.cutesv} --caller-vcf SVIM={input.svim} --locus-depth {input.locus_depth} "
        "--cn-segments {input.cn_segments} --genes {input.genes} "
        "--out-recon-dir {output.recon_dir} --out-recon-tsv {output.recon_tsv} "
        "> {log} 2>&1"
