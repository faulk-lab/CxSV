#!/bin/bash
# Stage C: score every panel sample (sample_score.sh) whose depth/catalogue stage
# is done, N in parallel on this node. Sex from families.tsv (registry).
# Usage: run_scores.sh N [MODES]     (env SCORE_TAG / WEIGHTS_OVERRIDE pass through;
#        SHARD=i/k scores only every k-th task starting at i, to split over nodes)
set -uo pipefail
here=$(cd "$(dirname "$0")" && pwd)
source "$here/config.sh"
n=${1:-16}; modes=${2:-"LOIO LOFO"}
tasks=$CAL_ROOT/score_tasks.$(hostname).txt
python3 - "$SAMPLE_MANIFEST" "$CAL_ROOT/families.tsv" "$CAL_ROOT" "${SHARD:-0/1}" > "$tasks" <<'PY'
import csv, os, sys
man, fam, root, shard = sys.argv[1:5]
i, k = map(int, shard.split("/"))
sex = {r["individual_id"]: r["sex"] for r in csv.DictReader(open(fam), delimiter="\t")}
for t, r in enumerate(csv.DictReader(open(man), delimiter="\t")):
    if t % k == i and os.path.exists(f"{root}/persample/{r['cohort']}/{r['sample_id']}/depth_catalogues.ok"):
        print(r["cohort"], r["sample_id"], sex.get(r["individual_id"]) or "unknown")
PY
echo "scoring $(wc -l < "$tasks") samples with $n workers"
xargs -P "$n" -L 1 -a "$tasks" bash -c 'bash '"$here"'/sample_score.sh "$0" "$1" "$2" "'"$modes"'"'
echo "finished $(date -Is)"
