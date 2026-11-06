#!/usr/bin/env bash
# depth_worker.sh TASKS ROOT PAR — run depth_pass.sh for every task line
# (key<TAB>alignment) not yet done or claimed, PAR at a time on this node.
# Claims are directories ROOT/claims/<key> (mkdir is atomic on the CIFS share);
# a finished task has ROOT/<key>/depth.ok. Several nodes may run this at once.
set -uo pipefail
TASKS=$1; ROOT=$2; PAR=${3:-2}; HERE="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$ROOT/claims" "$ROOT/logs"
run_one() {
    local key=$1 aln=$2
    [ -f "$ROOT/$key/depth.ok" ] && return 0
    mkdir -p "$ROOT/claims/$(dirname "$key")"
    mkdir "$ROOT/claims/$key" 2>/dev/null || return 0
    echo "$(hostname -s) $(date -Is)" > "$ROOT/claims/$key/host"
    if bash "$HERE/depth_pass.sh" "$aln" "$ROOT/$key" "${THREADS:-4}" > "$ROOT/logs/${key//\//__}.log" 2>&1; then
        echo "ok $key $(hostname -s)"
    else
        echo "FAILED $key $(hostname -s)"; rm -rf "$ROOT/claims/$key"
        mkdir -p "$ROOT/failed"; echo "$(hostname -s)" > "$ROOT/failed/${key//\//__}"
    fi
}
export -f run_one; export ROOT HERE
grep -v '^#' "$TASKS" | xargs -P "$PAR" -L 1 bash -c 'run_one "$0" "$1"'
