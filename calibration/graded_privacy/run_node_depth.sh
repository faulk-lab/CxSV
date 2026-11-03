#!/bin/bash
# Run sample_depth_catalogues.sh over a shard (cohort<TAB>sample<TAB>alignment per
# line) with N parallel workers on this node. Skips samples already done.
# Usage: run_node_depth.sh SHARD N
set -uo pipefail
here=$(cd "$(dirname "$0")" && pwd)
shard=$1; n=$2
awk -F'\t' '{print $1" "$2" "$3}' "$shard" | xargs -P "$n" -L 1 bash "$here/sample_depth_catalogues.sh"
echo "node $(hostname) shard $shard finished $(date -Is)"
