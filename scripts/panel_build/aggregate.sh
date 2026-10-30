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
# The dry run is checked first: it must schedule no per-sample job.
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
TMP="$HOME/scratch/cxsv_panel_build/tmp/aggregate"
run() { cxsv_container "$BUILD/aggregate" "$TMP" snakemake -s "$BUILD/code/workflow/panel/Snakefile" \
        --directory "$BUILD/aggregate" --configfile "$BUILD/aggregate.yaml" --cores "$CORES" \
        --keep-going --printshellcmds "$@"; }
mkdir -p "$BUILD/aggregate"
run -n > "$BUILD/aggregate/dryrun.txt" 2>&1
per_sample=$(awk '/^Job stats/{f=1} f && $1 ~ /^(run_|filter_|merge_sv|sample_|callable_mask|check_cram)/ {s+=$2} END{print s+0}' "$BUILD/aggregate/dryrun.txt")
grep -A40 '^Job stats' "$BUILD/aggregate/dryrun.txt" | sed '/^$/q'
[ "$per_sample" = 0 ] || { echo "dry run schedules $per_sample per-sample jobs: a shard is incomplete or stale; stop"; exit 1; }
$DRY && exit 0
run "${EXTRA[@]}" > "$BUILD/aggregate/aggregate.log" 2>&1
echo "aggregate done: $BUILD/results/panel/release/"
