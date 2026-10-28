# ═══════════════════════════════════════════════════════════════════════════════
# PANEL — CARRIER TABLES, JUNCTION SITES AND PANEL RELEASE
#
# 1. build_panel: stacks every background sample's CxSV loci and counts
#    distinct carrier individuals per locus (scripts/build_panel_loci.py).
# 2. build_panel_sites: clusters every sample's junctions (whole SV set,
#    rules/junctions.smk) across individuals into junction sites with carrier
#    individuals and allele frequencies by platform and superpopulation
#    (scripts/build_panel_sites.py; roadmap item 3). The patient workflow
#    matches its CxSV loci against these sites junction by junction.
# 3. panel_sniffles_population: Sniffles2 population-mode joint genotyping of
#    the per-sample .snf files; its per-site allele frequencies are kept
#    alongside and compared with the clustered sites
#    (scripts/sniffles_population_sites.py).
#
# The release directory is what the patient workflow needs, and all it needs:
# no reads, no per-sample VCFs.
#
# release/
#   panel_cxsv_calls.tsv            one row per (sample, CxSV locus): region-overlap
#                                   (hotspot) information for patient loci
#   panel_cxsv_loci.tsv             merged loci with carrier counts (summary/plots)
#   panel_summary.tsv               samples / individuals / calls per cohort, platform
#   panel_junction_sites.tsv        junction sites, shared junction-site schema
#   panel_junction_carriers.tsv     the same sites with envelopes + carrier individuals
#   panel_individuals.tsv           sex, superpopulation, platforms per individual
#   sniffles_population_sites.tsv   Sniffles2 population-mode sites, same schema
#   panel_frequency_comparison*.tsv Sniffles2 vs clustered carriers per site + summary
#   calling_fingerprint.json        what a patient run must reproduce
#   tool_versions.tsv
#   release_manifest.tsv            sha256 of every file above
# Versioned, self-describing releases and the offline bundle: roadmap item 2.
# ═══════════════════════════════════════════════════════════════════════════════

import hashlib
import shutil

P   = f"{RESULTS_ROOT}/panel"
REL = f"{RESULTS_ROOT}/release"


PANEL_SITES = config.get("panel_sites", {})
PANEL_INDIVIDUALS = config.get("panel_individuals") or ""


# Release label written into every junction-site row: panel_release_id when
# set, else "fp-" + the first 12 hex digits of the calling fingerprint
# (resolved by the scripts from --fingerprint).
PANEL_VERSION_ARG = (f"--version {config['panel_release_id']}" if config.get("panel_release_id")
                     else f"--fingerprint {RESULTS_ROOT}/calling_fingerprint.json")


rule panel_sample_manifest:
    """Which samples (and which individuals/platforms) the panel is built from."""
    input:
        tables    = per_sample(f"{S}/cxsv/cxsv_master_table.tsv"),
        junctions = per_sample(f"{S}/junctions/junctions.tsv")
    output: f"{P}/sample_manifest.tsv"
    run:
        with open(output[0], "w") as fh:
            fh.write("cohort\tsample_id\tindividual_id\tplatform\tcoverage_x\tmaster_table\tjunctions\n")
            for (cohort, sample), table, junc in zip(KEYS, input.tables, input.junctions):
                r = SAMPLE_TABLE[(cohort, sample)]
                fh.write(f"{cohort}\t{sample}\t{r.individual_id}\t{r.platform}\t{r.coverage_x}\t"
                         f"{table}\t{junc}\n")


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


rule build_panel_sites:
    """Cluster all panel junctions into sites (shared junction-site schema)."""
    input:
        manifest    = f"{P}/sample_manifest.tsv",
        fingerprint = f"{RESULTS_ROOT}/calling_fingerprint.json",
        fai         = f"{REFERENCE}.fai",
        script      = f"{SCRIPTS}/build_panel_sites.py"
    output:
        sites       = f"{P}/panel_junction_sites.tsv",
        carriers    = f"{P}/panel_junction_carriers.tsv",
        individuals = f"{P}/panel_individuals.tsv"
    params:
        tolerance   = PANEL_SITES.get("tolerance", 500),
        ins_ratio   = PANEL_SITES.get("ins_len_ratio", 0.5),
        individuals = f"--individuals {PANEL_INDIVIDUALS}" if PANEL_INDIVIDUALS else "",
        min_svlen   = config["sv_params"]["min_sv_len"],
        version     = PANEL_VERSION_ARG
    log: f"{P}/logs/build_panel_sites.log"
    shell:
        """
        python {input.script} --manifest {input.manifest} {params.individuals} --fai {input.fai} \
            --tolerance {params.tolerance} --ins-len-ratio {params.ins_ratio} \
            --min-svlen {params.min_svlen} --source cxsv_panel {params.version} --sites {output.sites} \
            --carriers {output.carriers} --panel-individuals {output.individuals} 2> {log}
        """


def sniffles_ids():
    """Unique Sniffles2 sample ids: the sample id, or cohort.sample when a
    sample id occurs in more than one cohort."""
    dup = {s for s in SAMPLE_OF_KEY if SAMPLE_OF_KEY.count(s) > 1}
    return [f"{c}.{s}" if s in dup else s for c, s in KEYS]


rule panel_sniffles_population:
    """
    Sniffles2 population mode: joint genotyping of every panel sample's .snf.
    sniffles_population.extra passes options for very large panels (e.g.
    --combine-max-inmemory-results / --no-sort; see sniffles --expert-help).
    """
    input:  snfs = per_sample(f"{S}/callers/sniffles.snf")
    output:
        tsv = f"{P}/sniffles_population/snf_inputs.tsv",
        map = f"{P}/sniffles_population/samples.tsv",
        vcf = f"{P}/sniffles_population/sniffles_population.vcf"
    params: extra = config.get("sniffles_population", {}).get("extra", "")
    threads: config["threads"]["sv_callers"]
    log: f"{P}/logs/sniffles_population.log"
    run:
        ids = sniffles_ids()
        with open(output.tsv, "w") as ft, open(output.map, "w") as fm:
            fm.write("snf\tsniffles_id\tindividual_id\n")
            for (cohort, sample), snf, sid in zip(KEYS, input.snfs, ids):
                ft.write(f"{os.path.abspath(snf)}\t{sid}\n")
                fm.write(f"{snf}\t{sid}\t{SAMPLE_TABLE[(cohort, sample)].individual_id}\n")
        shell("sniffles --input {output.tsv} --vcf {output.vcf} --threads {threads} "
              "{params.extra} > {log} 2>&1")


rule panel_sniffles_sites:
    """Sniffles2 population sites in the shared schema + comparison with clustered sites."""
    input:
        vcf         = f"{P}/sniffles_population/sniffles_population.vcf",
        map         = f"{P}/sniffles_population/samples.tsv",
        individuals = f"{P}/panel_individuals.tsv",
        carriers    = f"{P}/panel_junction_carriers.tsv",
        fingerprint = f"{RESULTS_ROOT}/calling_fingerprint.json",
        fai         = f"{REFERENCE}.fai",
        script      = f"{SCRIPTS}/sniffles_population_sites.py"
    output:
        sites   = f"{P}/sniffles_population_sites.tsv",
        compare = f"{P}/panel_frequency_comparison.tsv",
        summary = f"{P}/panel_frequency_comparison_summary.tsv"
    params:
        tolerance = PANEL_SITES.get("tolerance", 500),
        ins_ratio = PANEL_SITES.get("ins_len_ratio", 0.5),
        min_svlen = config["sv_params"]["min_sv_len"],
        version   = PANEL_VERSION_ARG
    log: f"{P}/logs/sniffles_population_sites.log"
    shell:
        """
        python {input.script} --vcf {input.vcf} --samples {input.map} \
            --panel-individuals {input.individuals} --carriers {input.carriers} --fai {input.fai} \
            --tolerance {params.tolerance} --ins-len-ratio {params.ins_ratio} --min-svlen {params.min_svlen} {params.version} \
            --sites {output.sites} --compare {output.compare} --summary {output.summary} 2> {log}
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
    "panel_junction_sites.tsv":    f"{P}/panel_junction_sites.tsv",
    "panel_junction_carriers.tsv": f"{P}/panel_junction_carriers.tsv",
    "panel_individuals.tsv":       f"{P}/panel_individuals.tsv",
    "sniffles_population_sites.tsv":          f"{P}/sniffles_population_sites.tsv",
    "panel_frequency_comparison.tsv":         f"{P}/panel_frequency_comparison.tsv",
    "panel_frequency_comparison_summary.tsv": f"{P}/panel_frequency_comparison_summary.tsv",
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
