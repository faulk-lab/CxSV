#!/usr/bin/env python3
"""Score a patient run's discovery with the copy-number-linked loci added (tests C and E).

The downstream rules of workflow/rules/patient_private.smk and patient_scoring.smk,
run script by script with the run's own parameters, on the augmented discovery tree
(link_all.py OUT/samples/<cohort>/<sample>):

  classify_patient_cxsv.py   against the panel release (privacy; item-3 matcher)
  match_catalogues.py        external catalogues
  locus_segments.py + samtools bedcov -Q min_mapq -j on the run's alignment; depth rows of
                             loci whose id is unchanged are copied from the run
                             (same junctions -> same segments), only new loci are measured
  score_cxsv.py              the run's scoring_params.json, caller VCFs, library QC,
                             SVIM check; no weight changes
Outputs go into the augmented tree (private/, scoring/), so score_pos_control.py can
read it like a run. Run inside the pipeline image (python, samtools; env REF_PATH for CRAM).

  patient_rescore.py --run RUN_SAMPLE_DIR --aug AUG_SAMPLE_DIR --aln ALIGNMENT \
      --release RELEASE --code CODE --ref hs1.fa [--catalogues ROOT]
"""
import argparse
import csv
import json
import os
import subprocess

import yaml


def sh(cmd, log):
    with open(log, "a") as fh:
        fh.write("$ " + cmd + "\n")
        fh.flush()
        subprocess.run(cmd, shell=True, check=True, stdout=fh, stderr=fh, executable="/bin/bash")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--aug", required=True)
    ap.add_argument("--aln", required=True)
    ap.add_argument("--release", required=True)
    ap.add_argument("--code", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--catalogues", default="/mnt/sargasso1/CxSV/datasets/external_catalogues")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--bundle", default="/mnt/sargasso1/CxSV/bundles/cxsv-bundle-20261105-8fc5bef725bd",
                    help="host path of the bundle the run used: /bundle/ paths in its scoring_params.json map here")
    a = ap.parse_args()
    R, A, W = a.run, a.aug, f"{a.code}/workflow/scripts"
    C = f"{a.code}/config"
    pat = yaml.safe_load(open(f"{C}/patient.yaml"))
    sc = yaml.safe_load(open(f"{C}/scoring.yaml"))
    cal = yaml.safe_load(open(f"{C}/calling.yaml"))
    m = pat.get("panel_match", {})
    for d in ("private", "scoring", "logs"):
        os.makedirs(f"{A}/{d}", exist_ok=True)
    log = f"{A}/logs/rescore.log"
    open(log, "w").close()
    Rel = a.release
    sh(f"python {W}/classify_patient_cxsv.py --master-table {A}/cxsv/cxsv_master_table.tsv "
       f"--junctions {A}/junctions/junctions.tsv --panel-carriers {Rel}/panel_junction_carriers.tsv "
       f"--panel-sites {Rel}/panel_junction_sites.tsv --panel-individuals {Rel}/panel_individuals.tsv "
       f"--panel-calls {Rel}/panel_cxsv_calls.tsv --tolerance {m.get('tolerance', 500)} "
       f"--ins-len-ratio {m.get('ins_len_ratio', 0.5)} --min-junction-fraction {m.get('min_junction_fraction', 0.75)} "
       f"--min-carriers {m.get('min_carriers', 1)} --window {m.get('hotspot_window', 1000)} "
       f"--release-dir {Rel} "
       f"--out-all {A}/private/cxsv_panel_status.tsv --out-private {A}/private/private_cxsv_master_table.tsv "
       f"--out-bed {A}/private/private_cxsv_loci.bed --out-junctions {A}/private/junction_panel_status.tsv", log)
    cats = " ".join(f"--catalogue {t['name']}:{t['version']}:{a.catalogues}/{t['name']}/{t['version']}/{t['name']}.{t['version']}.junction_sites.tsv.gz:"
                    f"{a.catalogues}/{t['name']}/{t['version']}/{t['name']}.{t['version']}.manifest.tsv"
                    for t in sc["external_catalogues"]["tables"])
    cm = sc["external_catalogues"].get("match", {})
    sh(f"python {W}/match_catalogues.py --junctions {A}/junctions/junctions.tsv --loci {A}/private/cxsv_panel_status.tsv "
       f"{cats} --tolerance {cm.get('tolerance', 500)} --ins-len-ratio {cm.get('ins_len_ratio', 0.5)} "
       f"--min-junction-fraction {cm.get('min_junction_fraction', 0.75)} --out-loci {A}/private/cxsv_catalogue_hits.tsv "
       f"--out-junctions {A}/private/junction_catalogue_hits.tsv --out-check {A}/private/catalogue_check.tsv", log)
    cn = sc["copy_number"]
    jp, cp = cal["junction_params"], cal["cxsv_params"]
    sh(f"python {W}/locus_segments.py --junctions {A}/junctions/junctions.tsv --loci {A}/private/cxsv_panel_status.tsv "
       f"--fai {a.ref}.fai --collapse-distance {jp['collapse_distance']} --composite-fraction {cp['composite_fraction']} "
       f"--flank {cn['flank_bp']} --gap {cn['gap_bp']} --margin {cn['margin_bp']} --out {A}/scoring/locus_segments.bed", log)
    old = {}
    if os.path.exists(f"{R}/scoring/locus_depth.tsv"):
        for line in open(f"{R}/scoring/locus_depth.tsv"):
            f = line.rstrip("\n").split("\t")
            old.setdefault(f[3], []).append(line)
    need, reuse = [], []
    for line in open(f"{A}/scoring/locus_segments.bed"):
        f = line.rstrip("\n").split("\t")
        if f[3] in old:
            continue
        need.append(line)
    seen = set()
    for line in open(f"{A}/scoring/locus_segments.bed"):
        lid = line.split("\t")[3]
        if lid in old and lid not in seen:
            reuse += old[lid]
            seen.add(lid)
    with open(f"{A}/scoring/segments_new.bed", "w") as fh:
        fh.writelines(need)
    if need:
        sh(f"samtools bedcov -Q {cn['min_mapq']} -j --reference {a.ref} {A}/scoring/segments_new.bed {a.aln} "
           f"> {A}/scoring/depth_new.tsv", log)
    else:
        open(f"{A}/scoring/depth_new.tsv", "w").close()
    with open(f"{A}/scoring/locus_depth.tsv", "w") as fh:
        fh.writelines(reuse)
        fh.write(open(f"{A}/scoring/depth_new.tsv").read())
    # the run's scoring parameters; paths inside the container (/bundle/...) -> host paths
    with open(f"{R}/scoring/scoring_params.json") as fh:
        ptxt = fh.read().replace("/bundle/", a.bundle.rstrip("/") + "/")
    with open(f"{A}/scoring/scoring_params.json", "w") as fh:
        fh.write(ptxt)
    lib = f"--library-qc {R}/qc/library_qc.tsv" if os.path.exists(f"{R}/qc/library_qc.tsv") else ""
    sh(f"python {W}/score_cxsv.py --status {A}/private/cxsv_panel_status.tsv --catalogue-hits {A}/private/cxsv_catalogue_hits.tsv "
       f"--catalogue-check {A}/private/catalogue_check.tsv --junctions {A}/junctions/junctions.tsv "
       f"--junction-status {A}/private/junction_panel_status.tsv --vcf {A}/cxsv/cxsv_only.vcf "
       f"--caller-vcf Sniffles2={A}/callers/sniffles_filtered.vcf --caller-vcf cuteSV={A}/callers/cutesv_filtered.vcf "
       f"--caller-vcf SVIM={A}/callers/svim_filtered.vcf --depth {A}/scoring/locus_depth.tsv "
       f"--params {A}/scoring/scoring_params.json --genes {sc_ann(cal)['genes_bed']} --exons {sc_ann(cal)['exons_bed']} "
       f"--repeats {sc_ann(cal)['repeats_bed']} --segdups {sc_ann(cal)['segdup_bed']} --fai {a.ref}.fai {lib} "
       f"--sample {a.sample} --out-tsv {A}/scoring/cxsv_scores.tsv --out-features {A}/scoring/cxsv_score_features.tsv "
       f"--out-report {A}/scoring/cxsv_score_report.md", log)
    print(json.dumps(dict(reused_loci=len(seen), new_segments=len(need))))


RES = os.environ.get("CXSV_RESOURCES", "/mnt/sargasso1/CxSV/resources")


def sc_ann(cal):
    """Annotation BEDs of config/calling.yaml (relative names) under the resources directory."""
    return {k: (v if str(v).startswith("/") else f"{RES}/{os.path.basename(str(v))}")
            for k, v in cal["annotation"].items() if k.endswith("_bed")}


if __name__ == "__main__":
    main()
