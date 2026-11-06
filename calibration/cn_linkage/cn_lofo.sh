#!/bin/bash
# =============================================================================
# cn_lofo.sh — test D: healthy-genome Tier 1 of the subset with the copy-number-linked
# loci added (roadmap item 15, Phase 1). The item-10 / item-13 LOFO calibration
# (calibration/graded_privacy, calibration/oversized_clusters/eval_depth.sh), unchanged,
# on the augmented discovery tree written by link_all.py:
#   1 locus_callability.py   leave-out callability of every locus of the subset genomes
#   2 leaveout_status.py     LOFO panel statuses against panel-v1.1 (family removed)
#   3 eval_depth.sh          read depth: rows of loci with an unchanged id are reused from
#                            the panel-v1.1 LOFO run (REUSE_CAL); only new loci are measured
#   4 sample_score.sh LOFO   deployed scorer, no weight changes
#
#   cn_lofo.sh AUG_ROOT CAL_ROOT [PAR]
# AUG_ROOT holds samples/ and subset_manifest.tsv (panel sample_manifest rows of the subset).
# =============================================================================
set -uo pipefail
AUG=$1; export CAL_ROOT=$2; PAR=${3:-8}
HERE="$(cd "$(dirname "$0")" && pwd)"
export CODE=${CODE:-$(cd "$HERE/../.." && pwd)}
export OLD_CAL=${OLD_CAL:-/mnt/sargasso1/CxSV/calibration/item10}
export REUSE_CAL=${REUSE_CAL:-/mnt/sargasso1/CxSV/calibration/panel_v11_r2/lofo}
NEWB=/mnt/sargasso1/CxSV/build/panel-v1.1-r2
export BUILD=$NEWB RELEASE=$NEWB/results/panel/release
export SAMPLES_ROOT=$AUG/samples SAMPLE_MANIFEST=$AUG/subset_manifest.tsv
export SIF=${SIF:-/mnt/sargasso1/CxSV/bundles/cxsv-bundle-20261105-95126e14f036/container/cxsv-pipeline.sif}
source "$CODE/calibration/graded_privacy/config.sh"
MAN=$NEWB/results/panel/panel/sample_manifest.tsv
mkdir -p "$CAL_ROOT/logs"
cp "$OLD_CAL/families.tsv" "$CAL_ROOT/families.tsv"
tail -n +2 "$SAMPLE_MANIFEST" | cut -f1,2 > "$CAL_ROOT/pairs.tsv"
flank=$(python3 -c "import json;print(json.load(open('$RELEASE/release.json'))['callability']['flank_bp'])")
t0=$(date +%s)
[ -s "$CAL_ROOT/locus_callable.tsv.gz" ] || $X python "$CODE/calibration/graded_privacy/locus_callability.py" \
   --sample-manifest "$SAMPLE_MANIFEST" --samples-root "$SAMPLES_ROOT" --callability-dir "$RELEASE/callability" \
   --flank "$flank" --families "$CAL_ROOT/families.tsv" --threads 16 --out "$CAL_ROOT/locus_callable.tsv.gz" \
   2> "$CAL_ROOT/logs/locus_callability.log" || exit 1
echo "locus callability wall_s $(( $(date +%s) - t0 ))"; t0=$(date +%s)
[ -f "$CAL_ROOT/leaveout.ok" ] || { $X python "$CODE/calibration/graded_privacy/leaveout_status.py" --release-dir "$RELEASE" \
   --sample-manifest "$MAN" --samples-root "$SAMPLES_ROOT" --families "$CAL_ROOT/families.tsv" \
   --locus-callable "$CAL_ROOT/locus_callable.tsv.gz" --modes LOFO --threads 8 --samples "$CAL_ROOT/pairs.tsv" \
   --out "$CAL_ROOT/leaveout" 2> "$CAL_ROOT/logs/leaveout_status.log" && touch "$CAL_ROOT/leaveout.ok"; } || exit 1
echo "leaveout wall_s $(( $(date +%s) - t0 ))"; t0=$(date +%s)
python3 - "$SAMPLE_MANIFEST" "$SAMPLE_SHEETS" > "$CAL_ROOT/depth_tasks.txt" <<'PY'
import csv, sys
man, sheets = sys.argv[1:3]
for r in csv.DictReader(open(man), delimiter="\t"):
    sh = {x["sample_id"]: x for x in csv.DictReader(open(f"{sheets}/{r['cohort']}.tsv"), delimiter="\t")}
    print(r["cohort"], r["sample_id"], sh[r["sample_id"]]["alignment"])
PY
xargs -P "$PAR" -L 1 bash -c 'bash '"$CODE"'/calibration/oversized_clusters/eval_depth.sh "$0" "$1" "$2"' < "$CAL_ROOT/depth_tasks.txt"
echo "depth: $(ls "$CAL_ROOT"/persample/*/*/depth_catalogues.ok | wc -l) wall_s $(( $(date +%s) - t0 ))"; t0=$(date +%s)
python3 - "$SAMPLE_MANIFEST" "$CAL_ROOT/families.tsv" > "$CAL_ROOT/score_tasks.txt" <<'PY'
import csv, sys
sex = {r["individual_id"]: r["sex"] for r in csv.DictReader(open(sys.argv[2]), delimiter="\t")}
for r in csv.DictReader(open(sys.argv[1]), delimiter="\t"):
    print(r["cohort"], r["sample_id"], sex.get(r["individual_id"]) or "unknown")
PY
xargs -P "$PAR" -L 1 bash -c 'bash '"$CODE"'/calibration/graded_privacy/sample_score.sh "$0" "$1" "$2" LOFO' < "$CAL_ROOT/score_tasks.txt"
echo "scored: $(ls "$CAL_ROOT"/persample/*/*/score.ok | wc -l) wall_s $(( $(date +%s) - t0 ))"
