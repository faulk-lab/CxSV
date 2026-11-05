#!/usr/bin/env python3
"""Re-run CxSV discovery (the per-sample rules after merge_sv_callers) for panel
genomes under a candidate discovery rule, into a separate tree.

R10-excess diagnosis (panel-v1.1 hold). Reads each genome's merged_sv.vcf from
the built panel (never writes there) and runs, with the CODE checkout's scripts
and the same parameters as the rules in workflow/rules/per_sample.smk and
junctions.smk (config/calling.yaml of CODE):

  vcf_to_bed.py (+ --single-caller-bnd-rule RULE) -> bedtools sort | cluster
  -> classify_cxsv.py -> extract_junctions.py -> extract_cxsv_vcf.py -> annotate_cxsv.py

Outputs per genome under OUT/samples/<cohort>/<sample>/{cxsv,junctions}; the three
filtered caller VCFs are hard-linked from the panel (read by the scorer). Writes
OUT/sample_manifest.tsv (the panel's, with this tree's paths). Run inside the
pipeline image (python deps, bedtools).

    rerun_discovery.py --panel PANEL_RESULTS --code CODE --out OUT --rule RULE
        [--svim-min-qual 10] [--samples cohort:sample,...] [--threads 16]
"""
import argparse
import csv
import os
import subprocess
import sys
from multiprocessing import Pool

import yaml


def sh(cmd, log):
    with open(log, "a") as fh:
        fh.write("$ " + cmd + "\n")
        fh.flush()
        subprocess.run(cmd, shell=True, check=True, stdout=fh, stderr=fh, executable="/bin/bash")


def one(t):
    a, cfg, cohort, sample = t
    src = f"{a.panel}/samples/{cohort}/{sample}"
    d = f"{a.out}/samples/{cohort}/{sample}"
    ok = f"{d}/rerun.ok"
    if os.path.exists(ok):
        return cohort, sample, "done"
    for sub in ("cxsv", "junctions", "callers", "logs"):
        os.makedirs(f"{d}/{sub}", exist_ok=True)
    for f in ("sniffles_filtered.vcf", "cutesv_filtered.vcf", "svim_filtered.vcf"):
        if not os.path.exists(f"{d}/callers/{f}"):
            os.link(f"{src}/callers/{f}", f"{d}/callers/{f}")
    S = f"{a.code}/workflow/scripts"
    c, j, an = cfg["cxsv_params"], cfg["junction_params"], cfg["annotation"]
    vcf, log = f"{src}/merged_sv.vcf", f"{d}/logs/rerun.log"
    fai = a.fai
    try:
        sh(f"python {S}/vcf_to_bed.py --vcf {vcf} --out {d}/cxsv/merged_sv.bed --window {c['bnd_window']} "
           f"--span-link-max {c['span_link_max']} --long-span-rule {c['long_span_rule']} "
           f"--long-out {d}/cxsv/long_span_records.tsv --single-caller-bnd-rule {a.rule} "
           f"--svim-min-qual {a.svim_min_qual} --filtered-out {d}/cxsv/single_caller_filtered.tsv", log)
        sh(f"bedtools sort -i {d}/cxsv/merged_sv.bed > {d}/cxsv/merged_sv_sorted.bed && "
           f"bedtools cluster -d {c['cluster_distance']} -i {d}/cxsv/merged_sv_sorted.bed > {d}/cxsv/clusters.bed", log)
        sh(f"python {S}/classify_cxsv.py --bed {d}/cxsv/clusters.bed --vcf {vcf} --fai {fai} "
           f"--output {d}/cxsv/cxsv_summary.tsv --collapse-distance {j['collapse_distance']} "
           f"--ins-len-ratio {j['ins_len_ratio']} --composite-fraction {c['composite_fraction']} "
           f"--min-bp {c['min_breakpoints']} --max-window {c['cluster_distance']} "
           f"--chromoanag-min {c['chromoanag_min_window']} --chromoanag-max {c['chromoanag_max_window']} "
           f"--chromoanag-min-sv {c['chromoanag_min_svs']} --extreme-min-breakends {c['extreme_min_breakends']} "
           f"--min-interleaved {c['min_interleaved']} --interleave-alt-min-intra {c['interleave_alt_min_intra']} "
           f"--interleave-alt-min-inter {c['interleave_alt_min_inter']} --max-cluster-size {c['max_cluster_size']}", log)
        sh(f"python {S}/extract_junctions.py --vcf {vcf} --clusters {d}/cxsv/clusters.bed "
           f"--summary {d}/cxsv/cxsv_summary.tsv --fai {fai} --collapse-distance {j['collapse_distance']} "
           f"--ins-len-ratio {j['ins_len_ratio']} --out {d}/junctions/junctions.tsv", log)
        sh(f"python {S}/extract_cxsv_vcf.py --vcf {vcf} --summary {d}/cxsv/cxsv_summary.tsv "
           f"--span-link-max {c['span_link_max']} --long-span-rule {c['long_span_rule']} --output {d}/cxsv/cxsv_only.vcf", log)
        sh(f"python {S}/annotate_cxsv.py --summary {d}/cxsv/cxsv_summary.tsv --vcf {d}/cxsv/cxsv_only.vcf "
           f"--genes {an['genes_bed']} --exons {an['exons_bed']} --repeats {an['repeats_bed']} "
           f"--segdups {an['segdup_bed']} --out-table {d}/cxsv/cxsv_master_table.tsv --out-bed {d}/cxsv/cxsv_loci.bed", log)
    except subprocess.CalledProcessError as e:
        return cohort, sample, f"FAILED {e.cmd[:80]}"
    open(ok, "w").close()
    return cohort, sample, "ok"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--panel", required=True, help="panel results root (…/results/panel)")
    ap.add_argument("--code", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--rule", required=True)
    ap.add_argument("--svim-min-qual", type=float, default=10)
    ap.add_argument("--aggregate-yaml", required=True, help="the panel build's aggregate.yaml (annotation paths)")
    ap.add_argument("--fai", default="/mnt/sargasso1/CxSV/refs/hs1.fa.fai")
    ap.add_argument("--samples", default=None)
    ap.add_argument("--threads", type=int, default=16)
    a = ap.parse_args()
    cfg = yaml.safe_load(open(f"{a.code}/config/calling.yaml"))
    cfg["annotation"] = yaml.safe_load(open(a.aggregate_yaml))["annotation"]
    man = list(csv.DictReader(open(f"{a.panel}/panel/sample_manifest.tsv"), delimiter="\t"))
    if a.samples:
        want = set(tuple(x.split(":")) for x in a.samples.split(","))
        man = [m for m in man if (m["cohort"], m["sample_id"]) in want]
    os.makedirs(a.out, exist_ok=True)
    tasks = [(a, cfg, m["cohort"], m["sample_id"]) for m in man]
    bad = 0
    with Pool(a.threads) as p:
        for i, (c, s, st) in enumerate(p.imap_unordered(one, tasks)):
            if st.startswith("FAILED"):
                bad += 1
                print(c, s, st, file=sys.stderr)
            if i % 200 == 0:
                print(i, file=sys.stderr)
    with open(f"{a.out}/sample_manifest.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(man[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        for m in man:
            m = dict(m)
            m["master_table"] = f"{a.out}/samples/{m['cohort']}/{m['sample_id']}/cxsv/cxsv_master_table.tsv"
            m["junctions"] = f"{a.out}/samples/{m['cohort']}/{m['sample_id']}/junctions/junctions.tsv"
            w.writerow(m)
    print(f"{len(tasks)} genomes, {bad} failed")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
