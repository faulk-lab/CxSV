# ═══════════════════════════════════════════════════════════════════════════════
# PATIENT — PANEL COMPATIBILITY, PRIVATE/PUBLIC LABELLING, PER-SAMPLE PLOTS
#
# Needs only a panel release directory (config: panel_release), never the
# background reads, so it can run inside a restricted environment.
#
# 0. The release must match its release_manifest.tsv (every listed file present
#    with its size and sha256, nothing unlisted; workflow/lib/cxsv_release.py).
#    Checked three times: when the workflow is parsed (a release that does not
#    verify stops the run before any job starts, dry runs included), by rule
#    verify_panel_release (the record, results/panel_release_check.txt), and
#    by classify_patient_cxsv right before it reads the release (a file changed
#    during a run, even with its mtime kept, is caught too).
# 1. check_panel_compatibility — this run's calling fingerprint must equal the
#    panel's (same parameters, reference, tool versions and loci code). A
#    mismatch stops the run: otherwise differences in how the two sides were
#    called would surface as false "private" CxSVs. Override only for
#    exploratory work with allow_fingerprint_mismatch: true (the mismatch is
#    then written into panel_compatibility.txt and every output carries it).
# 2. classify_patient_cxsv — breakpoint-level matching (roadmap item 3): each
#    patient CxSV locus's junctions are matched against the panel's junction
#    sites; a panel individual carries the locus when it carries at least
#    panel_match.min_junction_fraction of them. Status private /
#    shared_junctions / public, carrier count and frequency by platform and
#    superpopulation, plus region-overlap carriers as a hotspot flag. Every
#    patient junction (also outside CxSV loci) gets its panel carriers too.
#    Next to each carrier count: Panel_Callable / panel_callable, the panel
#    individuals callable at the locus's junctions (both breakpoints of every
#    junction), from the release's callability masks — so a private locus in
#    a region the panel could not see is not mistaken for a rare one.
# ═══════════════════════════════════════════════════════════════════════════════

import cxsv_release            # noqa: E402

PANEL_RELEASE = config.get("panel_release")
if not PANEL_RELEASE:
    raise ValueError("[patient] config['panel_release'] must point at a panel release directory")
PANEL_RELEASE = os.path.abspath(PANEL_RELEASE)
MATCH = config.get("panel_match", {})

_release_check = cxsv_release.verify(PANEL_RELEASE)
if not _release_check["ok"]:
    raise ValueError(
        f"[patient] panel release {PANEL_RELEASE} does not verify against its "
        f"{cxsv_release.MANIFEST}; refusing it ({len(_release_check['problems'])} problem(s)):\n  "
        + "\n  ".join(_release_check["problems"][:20]))


def release_listed_files():
    """Files the release manifest lists (and that exist): inputs of the check, so
    a changed file re-triggers it. Missing ones are reported by the check itself."""
    try:
        rows = cxsv_release.read_manifest(PANEL_RELEASE)
    except (OSError, ValueError):
        return []
    return [f"{PANEL_RELEASE}/{r['file']}" for r in rows if os.path.isfile(f"{PANEL_RELEASE}/{r['file']}")]


rule verify_panel_release:
    """Refuse a panel release whose files do not match its manifest."""
    input:
        manifest = f"{PANEL_RELEASE}/{cxsv_release.MANIFEST}",
        files    = release_listed_files(),
        script   = f"{SCRIPTS}/verify_panel_release.py"
    output: f"{RESULTS_ROOT}/panel_release_check.txt"
    log:    f"{RESULTS_ROOT}/logs/panel_release_check.log"
    shell:
        "python {input.script} {PANEL_RELEASE} --out {output} > {log} 2>&1 "
        "|| {{ cat {log} >&2; exit 1; }}"


rule check_panel_compatibility:
    """
    Compare this run's calling fingerprint with the panel release's. The full
    report goes to logs/panel_compatibility.log, which survives a failure
    (Snakemake deletes a failed job's outputs, not its log); the output file
    is a copy written only when the run may proceed.
    """
    input:
        ours   = f"{RESULTS_ROOT}/calling_fingerprint.json",
        theirs = f"{PANEL_RELEASE}/calling_fingerprint.json",
        check  = f"{RESULTS_ROOT}/panel_release_check.txt"
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
        table       = f"{S}/cxsv/cxsv_master_table.tsv",
        junctions   = f"{S}/junctions/junctions.tsv",
        carriers    = f"{PANEL_RELEASE}/panel_junction_carriers.tsv",
        sites       = f"{PANEL_RELEASE}/panel_junction_sites.tsv",
        individuals = f"{PANEL_RELEASE}/panel_individuals.tsv",
        calls       = f"{PANEL_RELEASE}/panel_cxsv_calls.tsv",
        release     = f"{PANEL_RELEASE}/{cxsv_release.METADATA}",
        check       = f"{RESULTS_ROOT}/panel_release_check.txt",
        compat      = f"{RESULTS_ROOT}/panel_compatibility.txt",
        script      = f"{SCRIPTS}/classify_patient_cxsv.py"
    output:
        all_loci  = f"{S}/private/cxsv_panel_status.tsv",
        private   = f"{S}/private/private_cxsv_master_table.tsv",
        bed       = f"{S}/private/private_cxsv_loci.bed",
        junctions = f"{S}/private/junction_panel_status.tsv"
    params:
        tolerance     = MATCH.get("tolerance", 500),
        ins_ratio     = MATCH.get("ins_len_ratio", 0.5),
        min_fraction  = MATCH.get("min_junction_fraction", 0.75),
        min_carriers  = MATCH.get("min_carriers", 1),
        window        = MATCH.get("hotspot_window", 1000)
    log: f"{S}/logs/classify_patient_cxsv.log"
    shell:
        """
        python {input.script} --master-table {input.table} --junctions {input.junctions} \
            --panel-carriers {input.carriers} --panel-sites {input.sites} \
            --panel-individuals {input.individuals} --panel-calls {input.calls} \
            --tolerance {params.tolerance} --ins-len-ratio {params.ins_ratio} \
            --min-junction-fraction {params.min_fraction} --min-carriers {params.min_carriers} \
            --window {params.window} --release-dir {PANEL_RELEASE} \
            --out-all {output.all_loci} --out-private {output.private} \
            --out-bed {output.bed} --out-junctions {output.junctions} 2> {log}
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
