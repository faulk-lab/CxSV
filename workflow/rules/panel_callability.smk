# ═══════════════════════════════════════════════════════════════════════════════
# PANEL — CALLABILITY  (panel builder only; roadmap item 2)
#
# Without callability, "absent from the panel" cannot be told apart from "the
# panel could not have seen it there". Per panel sample, mosdepth quantizes the
# depth of reads with MAPQ >= sv_params.min_mapq (the reads the callers use)
# into LOW / CALLABLE at panel_callability.min_depth; the CALLABLE intervals are
# kept as a gzipped BED mask (gzip -n: byte-reproducible). The masks ship in the
# release, so the patient workflow can ask how many panel individuals were
# callable at ANY junction, including private ones that have no panel site.
# Definitions: workflow/lib/cxsv_callability.py, docs/RELEASE_AND_BUNDLE.md.
#
#   1. callable_mask            per sample: {S}/callability/callable.bed.gz
#   2. panel_callability_masks  collects them as {P}/callability/ + samples.tsv
#   3. panel_site_callability   per junction site: n panel samples / individuals
#                               callable at both breakpoints
# Not fingerprinted: callability describes the panel's coverage, it does not
# decide which loci a genome has. Its parameters are recorded in release.json.
# ═══════════════════════════════════════════════════════════════════════════════

import shutil

import cxsv_callability       # noqa: E402

CALLABILITY = config.get("panel_callability", {}) or {}
CALL_MIN_DEPTH = int(CALLABILITY.get("min_depth", 10))
CALL_FLANK     = int(CALLABILITY.get("flank", 100))
CALL_THREADS   = int(CALLABILITY.get("threads", 4))
CALL_CRAM31    = str(CALLABILITY.get("cram31", "transcode"))      # transcode | fail
if CALL_CRAM31 not in ("transcode", "fail"):
    raise ValueError(f"panel_callability.cram31 must be transcode or fail, not {CALL_CRAM31!r}")
CALL_MAPQ      = int(config["sv_params"]["min_mapq"])
CALL_DIR       = f"{RESULTS_ROOT}/panel/callability"
CALL_MASKS     = {k: cxsv_callability.mask_name(*k) for k in KEYS}


rule callable_mask:
    """mosdepth quantized depth -> CALLABLE intervals of one panel sample
    (scripts/callable_mask.sh: CRAM 3.1 transcoded for mosdepth or refused, per
    panel_callability.cram31; decode errors or a 0-bp mask fail the step)."""
    input:
        unpack(caller_inputs),
        ref    = REFERENCE,
        script = f"{SCRIPTS}/callable_mask.sh"
    output:
        mask    = f"{S}/callability/callable.bed.gz",
        summary = f"{S}/callability/mosdepth.mosdepth.summary.txt"
    params:
        prefix    = f"{S}/callability/mosdepth",
        min_depth = CALL_MIN_DEPTH,
        mapq      = CALL_MAPQ,
        cram31    = CALL_CRAM31
    threads: CALL_THREADS
    log: f"{S}/logs/callable_mask.log"
    shell:
        HTS_ENV + """
        bash {input.script} {input.aln} {input.ref} {params.prefix} {output.mask} \
            {params.min_depth} {params.mapq} {threads} {params.cram31} > {log} 2>&1
        """


rule panel_callability_masks:
    """Collect every sample's mask under one directory with samples.tsv (release layout)."""
    input:
        masks = per_sample(f"{S}/callability/callable.bed.gz")
    output:
        samples = f"{CALL_DIR}/{cxsv_callability.SAMPLES_TSV}",
        masks   = [f"{CALL_DIR}/{CALL_MASKS[k]}" for k in KEYS]
    run:
        rows = []
        for (cohort, sample), src in zip(KEYS, input.masks):
            dst = f"{CALL_DIR}/{CALL_MASKS[(cohort, sample)]}"
            shutil.copyfile(src, dst)
            r = SAMPLE_TABLE[(cohort, sample)]
            bp = cxsv_callability.Mask.from_bed(dst).callable_bp()
            rows.append([CALL_MASKS[(cohort, sample)], cohort, sample, r.individual_id, r.platform, bp])
        with open(output.samples, "w") as fh:
            fh.write("\t".join(cxsv_callability.SAMPLES_COLUMNS) + "\n")
            for row in rows:
                fh.write("\t".join(map(str, row)) + "\n")


rule panel_site_callability:
    input:
        sites       = f"{RESULTS_ROOT}/panel/panel_junction_sites.tsv",
        individuals = f"{RESULTS_ROOT}/panel/panel_individuals.tsv",
        samples     = f"{CALL_DIR}/{cxsv_callability.SAMPLES_TSV}",
        masks       = [f"{CALL_DIR}/{CALL_MASKS[k]}" for k in KEYS],
        script      = f"{SCRIPTS}/panel_site_callability.py"
    output:
        sites   = f"{RESULTS_ROOT}/panel/panel_junction_callability.tsv",
        summary = f"{RESULTS_ROOT}/panel/panel_callability_summary.tsv"
    params: flank = CALL_FLANK
    log: f"{RESULTS_ROOT}/panel/logs/panel_site_callability.log"
    shell:
        """
        python {input.script} --sites {input.sites} --callability-dir {CALL_DIR} \
            --panel-individuals {input.individuals} --flank {params.flank} \
            --out {output.sites} --summary {output.summary} 2> {log}
        """
