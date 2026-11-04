#!/bin/bash
# Roadmap item 13, steps 3-4: re-run ONLY the downstream per-sample rules (vcf_to_bed,
# cluster, classify, junctions, extract, annotate, report) for a subset of panel
# genomes with the branch's code, on copies of their panel-v1.0 per-sample outputs,
# and show with Snakemake dry runs which jobs a rebuild would plan.
#
#   downstream_rerun.sh SUBSET.tsv WORKDIR CODE SIF [CORES]
#
# SUBSET.tsv: cohort, sample_id, ... (header). The panel-v1.0 build is only read:
# each sample directory is copied (rsync -a, mtimes kept; caller work directories
# left out) to WORKDIR/results/panel/samples/<cohort>/<sample>. Dry runs:
#   dryrun_plain.txt     no forced rule: with no per-sample Snakemake metadata in the
#                        work dir (as in BUILD/aggregate) nothing is planned
#   dryrun_forcerun.txt  --forcerun sample_vcf_to_bed: the downstream rules only
# then the forced run executes; benchmarks/<rule>.tsv give the per-sample cost.
# LONG_SPAN_RULE=exclude|breakends|paired sets cxsv_params.long_span_rule in eval.yaml
# (default: CODE's config/calling.yaml).
set -euo pipefail
SUB=$1; W=$2; CODE=$3; SIF=$4; CORES=${5:-16}
BUILD=/mnt/sargasso1/CxSV/build/panel-v1.0
SRC=$BUILD/results/panel/samples
R=$W/results/panel/samples
mkdir -p "$R" "$W/sheets" "$W/tmp"
tail -n +2 "$SUB" | while IFS=$'\t' read -r cohort sample rest; do
  mkdir -p "$R/$cohort"
  rsync -a --exclude 'callers/cutesv_work' --exclude 'callers/svim/signatures' --exclude 'callers/svim/candidates' \
        --exclude 'benchmarks' --exclude 'panel_build' "$SRC/$cohort/$sample" "$R/$cohort/"
done
python3 - "$SUB" "$W" <<'PY'
import csv, json, sys
sub, W = sys.argv[1], sys.argv[2]
want = {}
for r in csv.DictReader(open(sub), delimiter="\t"):
    want.setdefault(r["cohort"], set()).add(r["sample_id"])
cfg = json.load(open("/mnt/sargasso1/CxSV/build/panel-v1.0/aggregate.yaml"))
cfg["panel_results_root"] = f"{W}/results/panel"
for c, spec in list(cfg["panel_cohorts"].items()):
    if c not in want:
        del cfg["panel_cohorts"][c]; continue
    rows = list(csv.DictReader(open(spec["samples"]), delimiter="\t"))
    out = f"{W}/sheets/{c}.tsv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, rows[0].keys(), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(r for r in rows if r["sample_id"] in want[c])
    spec["samples"] = out
cfg["release"]["test_release"] = True
import os
if os.environ.get("LONG_SPAN_RULE"):
    cfg.setdefault("cxsv_params", {})["long_span_rule"] = os.environ["LONG_SPAN_RULE"]
json.dump(cfg, open(f"{W}/eval.yaml", "w"), indent=1)
with open(f"{W}/targets.txt", "w") as fh:
    for c, ss in want.items():
        for s in sorted(ss):
            for o in ("cxsv/cxsv_master_table.tsv", "junctions/junctions.tsv", "cxsv/cxsv_only.vcf", "sv_summary.txt"):
                fh.write(f"{W}/results/panel/samples/{c}/{s}/{o}\n")
PY
X="singularity exec -B /mnt/sargasso1 -B /home/minknow/scratch -B $W/tmp:/tmp $SIF"
SM="snakemake -s $CODE/workflow/panel_build/Snakefile --directory $W --configfile $W/eval.yaml --cores $CORES --printshellcmds"
$X $SM -n $(cat "$W/targets.txt") > "$W/dryrun_plain.txt" 2>&1 || true
$X $SM -n --forcerun sample_vcf_to_bed $(cat "$W/targets.txt") > "$W/dryrun_forcerun.txt" 2>&1 || true
for f in dryrun_plain dryrun_forcerun; do echo "== $f"; awk '/^Job stats/{f=1} f{print} f && /^$/{exit}' "$W/$f.txt"; grep -m3 -E "Nothing to be done|Error" "$W/$f.txt" || true; done
# keep the panel-v1.0 copies of the outputs that will be re-made, for comparison
tail -n +2 "$SUB" | while IFS=$'\t' read -r cohort sample rest; do
  d=$R/$cohort/$sample; mkdir -p "$d/v1.0"
  cp -p "$d/cxsv/cxsv_summary.tsv" "$d/cxsv/cxsv_master_table.tsv" "$d/cxsv/cxsv_only.vcf" "$d/junctions/junctions.tsv" "$d/v1.0/"
  sha256sum "$d"/callers/*.vcf "$d"/callers/svim/variants.vcf "$d/merged_sv.vcf" > "$d/v1.0/upstream.sha256"
done
t0=$(date +%s)
rc=0; $X $SM --forcerun sample_vcf_to_bed --keep-going $(cat "$W/targets.txt") > "$W/run.log" 2>&1 || rc=$?
echo "forced downstream run rc $rc wall_s $(( $(date +%s) - t0 ))"
tail -n +2 "$SUB" | while IFS=$'\t' read -r cohort sample rest; do
  d=$R/$cohort/$sample; sha256sum -c --quiet "$d/v1.0/upstream.sha256" || echo "UPSTREAM CHANGED $cohort/$sample"
done
echo "upstream check done"
