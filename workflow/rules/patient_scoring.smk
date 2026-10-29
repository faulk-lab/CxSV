# ═══════════════════════════════════════════════════════════════════════════════
# PATIENT — EXTERNAL CATALOGUE HITS AND SCORING  (roadmap items 4 and 5)
#
# Not fingerprinted: neither step decides which loci exist. Settings in
# config/scoring.yaml (docs/SCORING.md, docs/EXTERNAL_CATALOGUES.md).
#
#   private/cxsv_catalogue_hits.tsv     per locus: junctions matched / total, max AF,
#                                       contains yes/no for every catalogue
#   private/junction_catalogue_hits.tsv per locus junction: matched sites, max AF
#   private/catalogue_check.tsv         every catalogue table verified against its manifest
#   scoring/locus_segments.bed, scoring/locus_depth.tsv   read depth of locus segments
#   scoring/scoring_params.json         the weights, tiers and patient info used
#   scoring/cxsv_scores.tsv             ranked loci with every feature's value and points
#   scoring/cxsv_score_features.tsv     the same in long format
#   scoring/cxsv_score_report.md        human-readable per-call report
# ═══════════════════════════════════════════════════════════════════════════════

import csv as _csv
import json as _json

CATS = config.get("external_catalogues") or {}
CAT_MATCH = CATS.get("match", {})
SCORING_RES = config.get("scoring_resources") or {}
CN = config["copy_number"]


def catalogue_specs():
    """[(name, version, table, manifest)] from external_catalogues.root/tables."""
    root = CATS.get("root")
    if not root:
        return []
    out = []
    for t in CATS.get("tables") or []:
        d = os.path.join(root, t["name"], t["version"])
        out.append((t["name"], t["version"], os.path.join(d, f"{t['name']}.{t['version']}.junction_sites.tsv.gz"),
                    os.path.join(d, f"{t['name']}.{t['version']}.manifest.tsv")))
    return out


def scoring_resource_files():
    return {k: v for k, v in SCORING_RES.items() if v}


def patient_info(cohort, sample):
    """sex, hpo_terms, mother, father, parents_affected: sample-sheet columns, then
    config scoring 'patients' entries (which win)."""
    info = {}
    spec = COHORT_SPECS.get(cohort) or {}
    sheet = spec.get("samples") if isinstance(spec, dict) else None
    if sheet and os.path.exists(sheet):
        with open(sheet, newline="") as fh:
            for r in _csv.DictReader(fh, delimiter="\t"):
                if (r.get("sample_id") or "").strip() == sample:
                    for k in ("sex", "hpo_terms", "mother", "father", "parents_affected"):
                        if (r.get(k) or "").strip():
                            info[k] = r[k].strip()
    info.update((config.get("patients") or {}).get(sample) or {})
    terms = info.get("hpo_terms") or []
    if isinstance(terms, str):
        terms = [t.strip() for t in terms.replace(";", ",").split(",") if t.strip()]
    info["hpo_terms"] = sorted(set(terms))
    pa = info.get("parents_affected", False)
    info["parents_affected"] = pa if isinstance(pa, bool) else str(pa).lower() in ("yes", "true", "1")
    return info


def parent_junctions(wildcards):
    info = patient_info(wildcards.cohort, wildcards.sample)
    out = {}
    for role in ("mother", "father"):
        pid = info.get(role)
        if pid:
            if (wildcards.cohort, pid) not in SAMPLE_TABLE:
                raise ValueError(f"[scoring] {wildcards.sample}: {role} '{pid}' is not a sample of run "
                                 f"'{wildcards.cohort}'; parents must be in the same patient run")
            out[role] = f"{SAMPLES_ROOT}/{wildcards.cohort}/{pid}/junctions/junctions.tsv"
    return out


rule catalogue_hits:
    """Which external catalogue contains each patient CxSV locus (item-3 matcher)."""
    input:
        junctions = f"{S}/junctions/junctions.tsv",
        loci      = f"{S}/private/cxsv_panel_status.tsv",
        tables    = [t for _, _, t, _ in catalogue_specs()],
        manifests = [m for _, _, _, m in catalogue_specs()],
        script    = f"{SCRIPTS}/match_catalogues.py"
    output:
        loci      = f"{S}/private/cxsv_catalogue_hits.tsv",
        junctions = f"{S}/private/junction_catalogue_hits.tsv",
        check     = f"{S}/private/catalogue_check.tsv"
    params:
        cats      = " ".join(f"--catalogue {n}:{v}:{t}:{m}" for n, v, t, m in catalogue_specs()),
        tolerance = CAT_MATCH.get("tolerance", 500),
        ins_ratio = CAT_MATCH.get("ins_len_ratio", 0.5),
        fraction  = CAT_MATCH.get("min_junction_fraction", 0.75)
    log: f"{S}/logs/catalogue_hits.log"
    shell:
        "python {input.script} --junctions {input.junctions} --loci {input.loci} {params.cats} "
        "--tolerance {params.tolerance} --ins-len-ratio {params.ins_ratio} "
        "--min-junction-fraction {params.fraction} --out-loci {output.loci} "
        "--out-junctions {output.junctions} --out-check {output.check} 2> {log}"


rule locus_depth:
    """Mean read depth (MAPQ >= copy_number.min_mapq) over each locus's segments and flanks."""
    input:
        unpack(caller_inputs),
        junctions = f"{S}/junctions/junctions.tsv",
        loci      = f"{S}/private/cxsv_panel_status.tsv",
        fai       = f"{REFERENCE}.fai",
        script    = f"{SCRIPTS}/locus_segments.py"
    output:
        bed   = f"{S}/scoring/locus_segments.bed",
        depth = f"{S}/scoring/locus_depth.tsv"
    params:
        collapse  = config["junction_params"]["collapse_distance"],
        composite = config["cxsv_params"]["composite_fraction"],
        flank = CN["flank_bp"], gap = CN["gap_bp"], margin = CN["margin_bp"], mapq = CN["min_mapq"]
    log: f"{S}/logs/locus_depth.log"
    shell:
        HTS_ENV + """
        python {input.script} --junctions {input.junctions} --loci {input.loci} --fai {input.fai} \
            --collapse-distance {params.collapse} --composite-fraction {params.composite} \
            --flank {params.flank} --gap {params.gap} --margin {params.margin} --out {output.bed} 2> {log}
        if [ -s {output.bed} ]; then
            samtools bedcov -Q {params.mapq} -j --reference {REFERENCE} {output.bed} {input.aln} \
                > {output.depth} 2>> {log}
        else
            : > {output.depth}
        fi
        """


rule scoring_params:
    """The weights, tiers, resources and patient information this run scores with."""
    input:
        release = f"{PANEL_RELEASE}/{cxsv_release.METADATA}"
    output: f"{S}/scoring/scoring_params.json"
    run:
        with open(input.release) as fh:
            rel = _json.load(fh)
        body = {
            "weights": config["weights"], "tiers": config["tiers"], "copy_number": config["copy_number"],
            "clustering": config["clustering"], "foldback": config["foldback"],
            "scoring_resources": scoring_resource_files(),
            "panel_match": {"tolerance": MATCH.get("tolerance", 500),
                            "ins_len_ratio": MATCH.get("ins_len_ratio", 0.5),
                            "min_junction_fraction": MATCH.get("min_junction_fraction", 0.75)},
            "collapse_distance": config["junction_params"]["collapse_distance"],
            "composite_fraction": config["cxsv_params"]["composite_fraction"],
            "panel_release_id": rel.get("release_id"),
            "patient": patient_info(wildcards.cohort, wildcards.sample),
        }
        with open(output[0], "w") as fh:
            _json.dump(body, fh, indent=2, sort_keys=True)


rule score_cxsv:
    """Explainable points per locus, total, tier, report (scripts/score_cxsv.py)."""
    input:
        unpack(lambda wc: {f"parent_{k}": v for k, v in parent_junctions(wc).items()}),
        status    = f"{S}/private/cxsv_panel_status.tsv",
        cats      = f"{S}/private/cxsv_catalogue_hits.tsv",
        check     = f"{S}/private/catalogue_check.tsv",
        junctions = f"{S}/junctions/junctions.tsv",
        jstatus   = f"{S}/private/junction_panel_status.tsv",
        vcf       = f"{S}/cxsv/cxsv_only.vcf",
        sniffles  = f"{S}/callers/sniffles_filtered.vcf",
        cutesv    = f"{S}/callers/cutesv_filtered.vcf",
        svim      = f"{S}/callers/svim_filtered.vcf",
        depth     = f"{S}/scoring/locus_depth.tsv",
        params    = f"{S}/scoring/scoring_params.json",
        genes     = config["annotation"]["genes_bed"],
        exons     = config["annotation"]["exons_bed"],
        repeats   = config["annotation"]["repeats_bed"],
        segdups   = config["annotation"]["segdup_bed"],
        fai       = f"{REFERENCE}.fai",
        resources = list(scoring_resource_files().values()),
        script    = f"{SCRIPTS}/score_cxsv.py"
    output:
        tsv      = f"{S}/scoring/cxsv_scores.tsv",
        features = f"{S}/scoring/cxsv_score_features.tsv",
        report   = f"{S}/scoring/cxsv_score_report.md"
    params:
        parents = lambda wc: " ".join(f"--parent {k}={v}" for k, v in parent_junctions(wc).items())
    log: f"{S}/logs/score_cxsv.log"
    shell:
        "python {input.script} --status {input.status} --catalogue-hits {input.cats} "
        "--catalogue-check {input.check} --junctions {input.junctions} --junction-status {input.jstatus} "
        "--vcf {input.vcf} --caller-vcf Sniffles2={input.sniffles} --caller-vcf cuteSV={input.cutesv} "
        "--caller-vcf SVIM={input.svim} --depth {input.depth} --params {input.params} "
        "--genes {input.genes} --exons {input.exons} --repeats {input.repeats} --segdups {input.segdups} "
        "--fai {input.fai} {params.parents} --sample {wildcards.sample} "
        "--out-tsv {output.tsv} --out-features {output.features} --out-report {output.report} 2> {log}"
