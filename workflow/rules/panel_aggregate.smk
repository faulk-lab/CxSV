# ═══════════════════════════════════════════════════════════════════════════════
# PANEL — CARRIER TABLE AND PANEL RELEASE
#
# Stacks every background sample's CxSV loci and counts distinct carrier
# individuals per locus (workflow/scripts/build_panel_loci.py). The release
# directory is what the patient workflow needs, and all it needs: no reads,
# no per-sample VCFs.
#
# release/
#   panel_cxsv_calls.tsv      one row per (sample, CxSV locus): the match table
#   panel_cxsv_loci.tsv       merged loci with carrier counts (summary/plots)
#   panel_summary.tsv         samples / individuals / calls per cohort, platform
#   calling_fingerprint.json  what a patient run must reproduce
#   tool_versions.tsv
#   release_manifest.tsv      sha256 of every file above
# Versioned, self-describing releases and the offline bundle: roadmap item 2.
# ═══════════════════════════════════════════════════════════════════════════════

import hashlib
import shutil

P   = f"{RESULTS_ROOT}/panel"
REL = f"{RESULTS_ROOT}/release"


rule panel_sample_manifest:
    """Which samples (and which individuals/platforms) the panel is built from."""
    input:  per_sample(f"{S}/cxsv/cxsv_master_table.tsv")
    output: f"{P}/sample_manifest.tsv"
    run:
        with open(output[0], "w") as fh:
            fh.write("cohort\tsample_id\tindividual_id\tplatform\tcoverage_x\tmaster_table\n")
            for (cohort, sample), table in zip(KEYS, input):
                r = SAMPLE_TABLE[(cohort, sample)]
                fh.write(f"{cohort}\t{sample}\t{r.individual_id}\t{r.platform}\t{r.coverage_x}\t{table}\n")


rule build_panel:
    input:
        manifest = f"{P}/sample_manifest.tsv",
        script   = f"{SCRIPTS}/build_panel_loci.py"
    output:
        calls   = f"{P}/panel_cxsv_calls.tsv",
        loci    = f"{P}/panel_cxsv_loci.tsv",
        sharing = f"{P}/panel_cxsv_sharing.bed",
        summary = f"{P}/panel_summary.tsv"
    log: f"{P}/logs/build_panel.log"
    shell:
        """
        python {input.script} --manifest {input.manifest} --calls {output.calls} \
            --loci {output.loci} --sharing {output.sharing} --summary {output.summary} \
            2> {log}
        """


rule plot_panel_sharing:
    """How many individuals share each panel locus, and where they are."""
    input:  f"{P}/panel_cxsv_sharing.bed"
    output: f"{P}/plots/cxsv_population_sharing.pdf"
    log:    f"{P}/logs/plot_sharing.log"
    shell:  "Rscript {SCRIPTS}/plot_cxsv_population_sharing.R {input} {output} > {log} 2>&1"


RELEASE_FILES = {
    "panel_cxsv_calls.tsv":     f"{P}/panel_cxsv_calls.tsv",
    "panel_cxsv_loci.tsv":      f"{P}/panel_cxsv_loci.tsv",
    "panel_summary.tsv":        f"{P}/panel_summary.tsv",
    "calling_fingerprint.json": f"{RESULTS_ROOT}/calling_fingerprint.json",
    "tool_versions.tsv":        f"{RESULTS_ROOT}/tool_versions.tsv",
}


rule panel_release:
    input:  list(RELEASE_FILES.values())
    output:
        files    = [f"{REL}/{name}" for name in RELEASE_FILES],
        manifest = f"{REL}/release_manifest.tsv"
    run:
        with open(output.manifest, "w") as mf:
            mf.write("file\tbytes\tsha256\n")
            for name, src in RELEASE_FILES.items():
                dst = f"{REL}/{name}"
                shutil.copyfile(src, dst)
                h = hashlib.sha256()
                with open(dst, "rb") as fh:
                    for block in iter(lambda: fh.read(1 << 20), b""):
                        h.update(block)
                mf.write(f"{name}\t{os.path.getsize(dst)}\t{h.hexdigest()}\n")
