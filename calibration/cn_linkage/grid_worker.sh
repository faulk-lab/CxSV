#!/usr/bin/env bash
# =============================================================================
# grid_worker.sh — test B runs (grid.py run) across the fleet (roadmap item 15).
#
#   grid_worker.sh TASKS GRID SEG LOCAL [PAR] [THREADS]
#
# TASKS: lines "<cell> <cov>". GRID: grid.py --out (cells/, cells.tsv). SEG: segment_all
# output (pon/, segments_all.tsv). LOCAL: node-local work directory (the patient
# workflow runs there; CIFS timestamp skew breaks snakemake's mtime checks) and the
# node's copy of the panel release (LOCAL/release; copied once, file count and manifest
# checked; marker LOCAL/release.copy.ok).
# Claims: GRID/claims/<cell>__<cov> (atomic mkdir on the share); done = run.ok in the
# cell's cov directory; failures in GRID/failed/. Several nodes may run at once.
# Runs inside the pipeline image (python, samtools, snakemake, mosdepth).
# =============================================================================
set -uo pipefail
TASKS=$1; GRID=$2; SEG=$3; LOCAL=$4; PAR=${5:-4}; THREADS=${6:-4}
HERE="$(cd "$(dirname "$0")" && pwd)"; CODE=$(cd "$HERE/../.." && pwd)
SIF=${SIF:-/mnt/sargasso1/CxSV/bundles/cxsv-bundle-20261105-8fc5bef725bd/container/cxsv-pipeline.sif}
REL_SRC=${REL_SRC:-/mnt/sargasso1/CxSV/build/panel-v1.1-r2/results/panel/release}
mkdir -p "$GRID/claims" "$GRID/failed" "$GRID/logs" "$LOCAL/work"
# the marker sits beside the copy: the release verifies against its manifest and
# refuses any extra file
rm -f "$LOCAL/release/.copy.ok"
verify_release() {
    [ -d "$LOCAL/release" ] || return 1
    [ "$(find "$REL_SRC" -type f | wc -l)" = "$(find "$LOCAL/release" -type f | wc -l)" ] || return 1
    cmp -s "$REL_SRC/release_manifest.tsv" "$LOCAL/release/release_manifest.tsv"
}
if [ ! -f "$LOCAL/release.copy.ok" ]; then
    verify_release || { rm -rf "$LOCAL/release"; cp -r "$REL_SRC" "$LOCAL/release"; }
    verify_release || { echo "release copy incomplete or differs"; exit 1; }
    touch "$LOCAL/release.copy.ok"
fi
export GRID SEG LOCAL THREADS CODE SIF HERE
run_one() {
    local cell=$1 cov=$2 id="$1__cov$2"
    [ -f "$GRID/cells/$cell/cov$cov/run.ok" ] && return 0
    # a systematic error stops the fleet instead of burning through the task list
    [ "$(ls "$GRID/failed" | wc -l)" -ge "${MAX_FAIL:-10}" ] && return 0
    mkdir "$GRID/claims/$id" 2>/dev/null || return 0
    echo "$(hostname -s) $(date -Is)" > "$GRID/claims/$id/host"
    local t0=$(date +%s)
    if singularity exec -B /mnt/sargasso1 -B "$LOCAL" "$SIF" python3 "$HERE/grid.py" run --out "$GRID" --cell "$cell" \
         --cov "$cov" --code "$CODE" --seg "$SEG" --release "$LOCAL/release" --threads "$THREADS" \
         --local "$LOCAL/work" > "$GRID/logs/$id.log" 2>&1; then
        echo -e "ok\t$id\t$(hostname -s)\t$(( $(date +%s) - t0 ))"
    else
        echo -e "FAILED\t$id\t$(hostname -s)\t$(( $(date +%s) - t0 ))"; rm -rf "$GRID/claims/$id"
        hostname -s > "$GRID/failed/$id"
    fi
}
export -f run_one
grep -v '^#' "$TASKS" | xargs -P "$PAR" -L 1 bash -c 'run_one "$0" "$1"'
