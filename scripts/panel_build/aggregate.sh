#!/usr/bin/env bash
# =============================================================================
# aggregate.sh — the panel aggregate step over all shards, on freya.
#
#   bash scripts/panel_build/aggregate.sh BUILD [CORES] [--dry-run] [-- snakemake args]
#
# Requires every planned sample to be done (build_status.py --require-complete).
# Then runs workflow/panel/Snakefile in the image, in BUILD/aggregate (its own
# .snakemake), with BUILD/aggregate.yaml: all three cohort sheets in full and
# the shared results tree BUILD/results/panel. Every per-sample output already
# exists, so only the aggregate rules run: cohort call sets and QC plots,
# panel loci and junction sites with carriers per individual, callability
# masks and per-site callability, Sniffles2 population genotyping, release.
# The dry run is checked first: it must schedule no per-sample job of the shard
# stage. Since panel-v1.2 (copy-number track, docs/CN_LINKAGE.md) the per-sample
# rules after the aggregate panel-of-normals run here: AFTER_PON below.
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD="$1"; shift; CORES=24; DRY=false; EXTRA=()
while [ $# -gt 0 ]; do
    case "$1" in --dry-run) DRY=true; shift ;; --) shift; EXTRA=("$@"); break ;; *) CORES="$1"; shift ;; esac
done
export CXSV_SIF="$(python3 -c "import json; print(json.load(open('$BUILD/plan/plan.json'))['fleet']['sif'])")"
# shellcheck source=lib.sh
source "$HERE/lib.sh"
python3 "$BUILD/code/scripts/panel_build/build_status.py" "$BUILD" --write --require-complete
# Sniffles2 population genotyping (panel_sniffles_population) opens every sample's .snf in each
# worker: with the usual soft limit of 1024 open files, workers fail task after task with
# "OSError: [Errno 24] Too many open files" while sniffles carries on and exits 0 with whole
# regions missing. Raise the soft limit to the hard limit and require headroom over the panel size.
ulimit -Sn "$(ulimit -Hn)"
NSAMP=$(python3 -c "import csv; print(sum(1 for _ in csv.DictReader(open('$BUILD/plan/shard_plan.tsv'), delimiter='\t')))")
[ "$(ulimit -Sn)" -ge $((4 * NSAMP + 1024)) ] || { echo "open-file limit $(ulimit -Sn) too low for $NSAMP samples (need >= $((4 * NSAMP + 1024))); raise the hard limit (nofile)"; exit 1; }
echo "open-file limit: $(ulimit -Sn) ($NSAMP samples)"
TMP="$HOME/scratch/cxsv_panel_build/tmp/aggregate"
run() { cxsv_container "$BUILD/aggregate" "$TMP" snakemake -s "$BUILD/code/workflow/panel/Snakefile" \
        --directory "$BUILD/aggregate" --configfile "$BUILD/aggregate.yaml" --cores "$CORES" \
        --keep-going --printshellcmds "$@"; }
mkdir -p "$BUILD/aggregate"
run -n > "$BUILD/aggregate/dryrun.txt" 2>&1
# per-sample rules that need the aggregate PoN (rules/cn.smk and everything after the linkage)
AFTER_PON='cn_segment|sample_cn_linkage|sample_annotate_cxsv|sample_junctions|sample_sv_summary'
per_sample=$(awk -v ap="^($AFTER_PON)\$" '/^Job stats/{f=1} f && $1 ~ /^(run_|filter_|merge_sv|sample_|callable_mask|check_cram|cn_depth)/ && $1 !~ ap {s+=$2} END{print s+0}' "$BUILD/aggregate/dryrun.txt")
# print the first Job stats table without a pipe: under pipefail, grep | sed '/^$/q' dies of SIGPIPE
# (rc 141) when Snakemake prints a second Job stats block after sed has quit
awk '/^Job stats/{f=1} f{print} f && /^$/{exit}' "$BUILD/aggregate/dryrun.txt"
[ "$per_sample" = 0 ] || { echo "dry run schedules $per_sample per-sample jobs: a shard is incomplete or stale; stop"; exit 1; }
$DRY && exit 0
run "${EXTRA[@]}" > "$BUILD/aggregate/aggregate.log" 2>&1
# a tool that logs worker failures but exits 0 must not pass: fail on any such line in the aggregate logs
R="$(python3 -c "import json; print(json.load(open('$BUILD/aggregate.yaml')).get('panel_results_root', '$BUILD/results/panel'))")"
bad=$(grep -lE 'Too many open files|Errno 24|ERROR sniffles\.worker' "$R"/panel/logs/*.log "$R"/cohorts/*/logs/*.log 2>/dev/null || true)
[ -z "$bad" ] || { echo "worker errors in: $bad"; exit 1; }
echo "aggregate done: $BUILD/results/panel/release/"
