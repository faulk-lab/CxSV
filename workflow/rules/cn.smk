# ═══════════════════════════════════════════════════════════════════════════════
# COPY-NUMBER TRACK  (shared — fingerprinted; roadmap item 15, docs/CN_LINKAGE.md)
#
# A second discovery track beside junction clustering: copy-number segments from
# binned read depth against a platform panel-of-normals (PoN), with every SV record
# whose breakend falls at a segment boundary attached (distance-free), become
# copy-number-linked CxSV loci (DUP-TRP/INV-DUP and other CN-linked classes). It adds
# loci and removes none; the scorer is unchanged.
#
#   cn_depth            per sample: {S}/cn/depth.regions.bed.gz (mosdepth 1 kb, MAPQ;
#                       cn_params.depth), or the precomputed full-panel pass
#                       (config cn_depth_precomputed: <root>/<cohort>/<sample>/)
#   panel_cn_pon        panel: {P}/cn/pon/<stratum>.npz from every panel library not
#                       in cn_params.pon_exclude (library-QC flags), + pon_members.tsv;
#                       chrX one reference per sex call, chrY males only (ambiguous
#                       sex calls: autosomes only; lib/cxsv_cn.classify_sex)
#   cn_segment          per sample: {S}/cn/segments.tsv, stats.json against the PoN of
#                       its stratum (panel: {P}/cn/pon; patient: the release's cn/pon)
#   panel_cn_segments   panel: {P}/cn/panel_cn_segments.tsv (every library's segments;
#                       the segment-level panel frequency of the linkage step)
#   sample_cn_linkage   per sample: junction-track discovery + CN loci ->
#                       cxsv/{clusters.bed,cxsv_summary.tsv,cxsv_only.vcf},
#                       cn/{cn_candidates,cn_loci,cn_only}.tsv
#
# Stratum of a sample (cn_stratum): sample-sheet column cn_stratum, then qc_group
# (the library-QC calibration group: the same four strata), then patients.<sample>.
# qc_group, then the cohort spec's qc_group, then cn_params.chemistry_stratum[sheet
# column chemistry], then cn_params.default_stratum[platform].
# ═══════════════════════════════════════════════════════════════════════════════

import csv as _cn_csv
import json as _cn_json

CNP = config["cn_params"]
CN_PRE = config.get("cn_depth_precomputed") or "none"        # site setting, not fingerprinted
_CN_SHEET_CACHE = {}


def _cn_path(p):
    return p if os.path.isabs(p) else os.path.join(PROJECT_DIR, p)


def _cn_sheet_row(cohort, sample):
    spec = COHORT_SPECS.get(cohort) or {}
    sheet = spec.get("samples") if isinstance(spec, dict) else None
    if not sheet or not os.path.exists(sheet):
        return {}
    if sheet not in _CN_SHEET_CACHE:
        with open(sheet, newline="") as fh:
            _CN_SHEET_CACHE[sheet] = {(r.get("sample_id") or "").strip(): r for r in _cn_csv.DictReader(fh, delimiter="\t")}
    return _CN_SHEET_CACHE[sheet].get(sample, {})


def _cn_strata_table():
    """config cn_strata_table: TSV cohort, sample_id, cn_stratum (run configuration, e.g.
    calibration/cn_linkage/panel_strata_v1.tsv for the 1kGP panel, whose sample sheets carry
    no chemistry column). Empty when unset."""
    path = config.get("cn_strata_table")
    if not path:
        return {}
    path = _cn_path(path)
    with open(path, newline="") as fh:
        return {(r["cohort"], r["sample_id"]): r["cn_stratum"].strip()
                for r in _cn_csv.DictReader(fh, delimiter="\t") if (r.get("cn_stratum") or "").strip()}


CN_STRATA_TABLE = _cn_strata_table()


def cn_stratum_explicit(cohort, sample):
    """PoN stratum set by the site, '' when none: sheet cn_stratum / qc_group, then
    cn_strata_table, patients.<sample>.qc_group, cohort-spec qc_group,
    cn_params.chemistry_stratum[sheet chemistry]."""
    row = _cn_sheet_row(cohort, sample)
    for col in ("cn_stratum", "qc_group"):
        if (row.get(col) or "").strip():
            return row[col].strip()
    if (cohort, sample) in CN_STRATA_TABLE:
        return CN_STRATA_TABLE[(cohort, sample)]
    g = ((config.get("patients") or {}).get(sample) or {}).get("qc_group")
    if g:
        return str(g)
    spec = COHORT_SPECS.get(cohort) or {}
    if isinstance(spec, dict) and spec.get("qc_group"):
        return str(spec["qc_group"])
    chem = (row.get("chemistry") or "").strip()
    if chem and chem in (CNP.get("chemistry_stratum") or {}):
        return CNP["chemistry_stratum"][chem]
    return ""


def cn_default_stratum(cohort, sample):
    plat = (SAMPLE_TABLE[(cohort, sample)].platform or "").upper()
    s = (CNP.get("default_stratum") or {}).get(plat)
    if not s:
        raise ValueError(f"[cn] no copy-number stratum for {cohort}/{sample} (platform {plat!r}); "
                         "set sample-sheet column cn_stratum or cn_params.default_stratum")
    return s


def cn_stratum_of(cohort, sample):
    """Panel: the explicit stratum, else cn_params.default_stratum[platform]. Patient runs
    resolve it at run time instead (rule cn_segment): explicit, else the library-QC group
    detected from the alignment (qc/library_qc.tsv, lib/cxsv_libgroup.py), else the
    platform default."""
    return cn_stratum_explicit(cohort, sample) or cn_default_stratum(cohort, sample)


def cn_stratum(wildcards):
    return cn_stratum_of(wildcards.cohort, wildcards.sample)


CN_STRATA = sorted({cn_stratum_of(c, s) for c, s in KEYS}) if WORKFLOW_NAME == "panel" else []


def cn_params_body():
    return _cn_json.dumps({"cn_params": CNP, "cxsv_params": config["cxsv_params"],
                           "junction_params": config["junction_params"]}, sort_keys=True)


rule cn_params_json:
    """The fingerprinted cn_params block plus the discovery parameters the linkage
    step re-uses (one file, so a change re-runs the track)."""
    output: f"{RESULTS_ROOT}/cn_params.json"
    params: body = lambda wc: cn_params_body()
    run:
        with open(output[0], "w") as fh:
            fh.write(params.body + "\n")


rule cn_depth:
    """1 kb binned depth of one sample (scripts/cn_depth.sh)."""
    input:
        unpack(caller_inputs),
        ref    = REFERENCE,
        script = f"{SCRIPTS}/cn_depth.sh"
    output:
        depth  = f"{S}/cn/depth.regions.bed.gz",
        sha    = f"{S}/cn/depth.sha256"
    params:
        out  = f"{S}/cn",
        bin  = CNP["depth"]["bin_bp"],
        mapq = CNP["depth"]["min_mapq"],
        pre  = lambda wc: f"{CN_PRE}/{wc.cohort}/{wc.sample}" if CN_PRE != "none" else "none"
    threads: 4
    log: f"{S}/logs/cn_depth.log"
    shell:
        HTS_ENV + "bash {input.script} {input.aln} {input.ref} {params.out} {params.bin} {params.mapq} "
        "{threads} {params.pre} > {log} 2>&1"


def cn_pon_input(wildcards):
    """Panel: the stratum's PoN. Patient: every PoN of the release (the stratum is
    resolved at run time; cn_segment.py picks <stratum>.npz)."""
    if WORKFLOW_NAME == "panel":
        return [f"{RESULTS_ROOT}/panel/cn/pon/{cn_stratum(wildcards)}.npz"]
    d = f"{PANEL_RELEASE}/cn/pon"
    have = sorted(x for x in os.listdir(d) if x.endswith(".npz")) if os.path.isdir(d) else []
    if not have:
        raise ValueError(f"[cn] panel release {PANEL_RELEASE} has no copy-number PoN (cn/pon/): it predates "
                         "the copy-number track (roadmap item 15; panel-v1.2 or later is needed)")
    st = cn_stratum_explicit(wildcards.cohort, wildcards.sample)
    if st and f"{st}.npz" not in have:
        raise ValueError(f"[cn] the panel release has no PoN for stratum {st!r}; set the sample's "
                         f"cn_stratum (or qc_group) to one of {[x[:-4] for x in have]}")
    return [f"{d}/{x}" for x in have]


def cn_segment_extra(wildcards):
    if WORKFLOW_NAME == "panel":
        return {}
    return {"libqc": f"{S}/qc/library_qc.tsv"}


def cn_panel_segments(wildcards=None):
    if WORKFLOW_NAME == "panel":
        return f"{RESULTS_ROOT}/panel/cn/panel_cn_segments.tsv"
    return f"{PANEL_RELEASE}/cn/panel_cn_segments.tsv"


if WORKFLOW_NAME == "panel":
    rule panel_cn_pon:
        """PoN per stratum from the panel libraries' depth (scripts/cn_pon.py)."""
        input:
            depth   = per_sample(f"{S}/cn/depth.regions.bed.gz"),
            params  = f"{RESULTS_ROOT}/cn_params.json",
            script  = f"{SCRIPTS}/cn_pon.py"
        output:
            pons    = [f"{RESULTS_ROOT}/panel/cn/pon/{s}.npz" for s in CN_STRATA],
            members = f"{RESULTS_ROOT}/panel/cn/pon/pon_members.tsv",
            table   = f"{RESULTS_ROOT}/panel/cn/pon_input.tsv"
        threads: 16
        log: f"{RESULTS_ROOT}/panel/logs/cn_pon.log"
        run:
            with open(output.table, "w") as fh:
                fh.write("key\tindividual_id\tstratum\tdepth\n")
                for c, s in KEYS:
                    fh.write(f"{c}/{s}\t{SAMPLE_TABLE[(c, s)].individual_id}\t{cn_stratum_of(c, s)}\t"
                             f"{SAMPLES_ROOT}/{c}/{s}/cn/depth.regions.bed.gz\n")
            shell("python {input.script} --members {output.table} "
                  "--params {input.params} --out-dir {RESULTS_ROOT}/panel/cn/pon --threads {threads} 2> {log}")


rule cn_segment:
    """Copy-number segments of one sample against its stratum's PoN (scripts/cn_segment.py)."""
    input:
        unpack(cn_segment_extra),
        depth  = f"{S}/cn/depth.regions.bed.gz",
        pon    = cn_pon_input,
        params = f"{RESULTS_ROOT}/cn_params.json",
        script = f"{SCRIPTS}/cn_segment.py"
    output:
        segs  = f"{S}/cn/segments.tsv",
        stats = f"{S}/cn/stats.json"
    params:
        pon     = lambda wc, input: input.pon[0] if WORKFLOW_NAME == "panel" else os.path.dirname(input.pon[0]),
        stratum = lambda wc: cn_stratum(wc) if WORKFLOW_NAME == "panel" else (cn_stratum_explicit(wc.cohort, wc.sample) or "auto"),
        default = lambda wc: cn_default_stratum(wc.cohort, wc.sample),
        libqc   = lambda wc, input: f"--library-qc {input.libqc}" if WORKFLOW_NAME != "panel" else ""
    log: f"{S}/logs/cn_segment.log"
    shell:
        "python {input.script} --depth {input.depth} --pon {params.pon} --params {input.params} "
        "--stratum {params.stratum} --default-stratum {params.default} {params.libqc} "
        "--out {output.segs} --stats {output.stats} 2> {log}"


if WORKFLOW_NAME == "panel":
    rule panel_cn_segments:
        """Every panel library's segments in one table (segment-level panel frequency). The
        sex column is the library's sex call (male / female / ambiguous): chrX / chrY
        segments count only towards samples of the same sex (cxsv_cnlink.cn_public)."""
        input:
            segs  = per_sample(f"{S}/cn/segments.tsv"),
            stats = per_sample(f"{S}/cn/stats.json")
        output: f"{RESULTS_ROOT}/panel/cn/panel_cn_segments.tsv"
        run:
            tmp = output[0] + ".tmp"
            with open(tmp, "w") as fo:
                head = None
                for (c, s), f, sf in zip(KEYS, input.segs, input.stats):
                    st = _cn_json.load(open(sf))
                    with open(f) as fh:
                        h = fh.readline().rstrip("\n")
                        if head is None:
                            head = h
                            fo.write("key\tindividual_id\tstratum\tsex\t" + h + "\n")
                        for line in fh:
                            fo.write(f"{c}/{s}\t{SAMPLE_TABLE[(c, s)].individual_id}\t{st['stratum']}\t"
                                     f"{st.get('sex_call', st['sex'])}\t{line}")
                if head is None:
                    fo.write("key\tindividual_id\tstratum\tsex\tchrom\tstart\tend\tstate\tploidy\n")
            os.replace(tmp, output[0])


rule sample_cn_linkage:
    """Junction-track discovery + copy-number-linked loci (scripts/cn_linkage.py)."""
    input:
        vcf      = f"{S}/merged_sv.vcf",
        clusters = f"{S}/cxsv/junction_clusters.bed",
        summary  = f"{S}/cxsv/junction_summary.tsv",
        cxsv_vcf = f"{S}/cxsv/junction_only.vcf",
        single   = f"{S}/cxsv/single_caller_filtered.tsv",
        segs     = f"{S}/cn/segments.tsv",
        stats    = f"{S}/cn/stats.json",
        panel    = cn_panel_segments,
        params   = f"{RESULTS_ROOT}/cn_params.json",
        fai      = f"{REFERENCE}.fai",
        script   = f"{SCRIPTS}/cn_linkage.py",
        lib      = [f"{WORKFLOW_DIR}/lib/cxsv_cnlink.py"]
    output:
        clusters   = f"{S}/cxsv/clusters.bed",
        summary    = f"{S}/cxsv/cxsv_summary.tsv",
        cxsv_vcf   = f"{S}/cxsv/cxsv_only.vcf",
        candidates = f"{S}/cn/cn_candidates.tsv",
        loci       = f"{S}/cn/cn_loci.tsv",
        cn_only    = f"{S}/cn/cn_only.tsv"
    params:
        work       = f"{S}/cn/work",
        individual = lambda wc: SAMPLE_TABLE[(wc.cohort, wc.sample)].individual_id
    log: f"{S}/logs/cn_linkage.log"
    shell:
        "python {input.script} --vcf {input.vcf} --clusters {input.clusters} --summary {input.summary} "
        "--cxsv-vcf {input.cxsv_vcf} --single {input.single} --segments {input.segs} "
        "--panel-segments {input.panel} --individual {params.individual:q} --stats {input.stats} --fai {input.fai} "
        "--calling {input.params} --scripts {SCRIPTS} --work {params.work} "
        "--out-clusters {output.clusters} --out-summary {output.summary} --out-vcf {output.cxsv_vcf} "
        "--out-candidates {output.candidates} --out-loci {output.loci} --out-cn-only {output.cn_only} 2> {log}"
