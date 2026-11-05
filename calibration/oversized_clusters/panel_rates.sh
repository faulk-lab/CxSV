#!/bin/bash
# =============================================================================
# panel_rates.sh — healthy-genome Tier-1/Tier-2 rates of a rebuilt panel
# (panel-v1.1, roadmap item 13): every panel genome's own CxSV loci scored by
# the deployed patient scorer (workflow/scripts/score_cxsv.py, config/scoring.yaml
# with the item-10 graded privacy, no re-tuning) against the NEW release with
# the genome's family removed (LOFO), exactly as the item-10 calibration
# (calibration/graded_privacy) did for panel-v1.0.
#
#   panel_rates.sh NEW_BUILD CAL_ROOT PAR
#
# Steps (calibration/graded_privacy scripts of this checkout):
#   1 locus_callability.py   per-locus leave-out callability of the new loci
#   2 leaveout_status.py     LOFO panel statuses against the new release
#   3 eval_depth.sh          read depth of the loci that are not public; a locus
#                            whose id is unchanged has the same junctions, so its
#                            depth rows are reused from the item-10 calibration
#                            (OLD_CAL) or the item-13 subset evaluation (REUSE_CAL);
#                            only new loci are measured; catalogue hits for all
#   4 sample_score.sh LOFO   scores_LOFO.tsv per genome
# Families and sexes: the item-10 families.tsv (registry). Results are collected
# by panel_rates_collect.py.
# =============================================================================
set -uo pipefail
NEWB=$1; export CAL_ROOT=$2; PAR=${3:-12}
HERE="$(cd "$(dirname "$0")" && pwd)"
export CODE=${CODE:-$NEWB/code}
export OLD_CAL=${OLD_CAL:-/mnt/sargasso1/CxSV/calibration/item10}
# RELEASE / SAMPLES_ROOT / SAMPLE_MANIFEST may be set to score another tree of per-sample
# CxSV outputs against the same release (e.g. a candidate discovery rule, rerun_discovery.py)
export BUILD=$NEWB RELEASE=${RELEASE:-$NEWB/results/panel/release}
export SAMPLES_ROOT=${SAMPLES_ROOT:-$NEWB/results/panel/samples}
export SAMPLE_MANIFEST=${SAMPLE_MANIFEST:-$NEWB/results/panel/panel/sample_manifest.tsv}
export SIF=${SIF:?set SIF}
source "$CODE/calibration/graded_privacy/config.sh"
mkdir -p "$CAL_ROOT/logs"
cp "$OLD_CAL/families.tsv" "$CAL_ROOT/families.tsv"
tail -n +2 "$SAMPLE_MANIFEST" | cut -f1,2 > "$CAL_ROOT/pairs.tsv"
flank=$(python3 -c "import json;print(json.load(open('$RELEASE/release.json'))['callability']['flank_bp'])")
t0=$(date +%s)
[ -s "$CAL_ROOT/locus_callable.tsv.gz" ] || $X python "$CODE/calibration/graded_privacy/locus_callability.py" \
   --sample-manifest "$SAMPLE_MANIFEST" --samples-root "$SAMPLES_ROOT" --callability-dir "$RELEASE/callability" \
   --flank "$flank" --families "$CAL_ROOT/families.tsv" --threads 16 --out "$CAL_ROOT/locus_callable.tsv.gz" \
   ${LOCUS_CALLABILITY_EXTRA:-} \
   2> "$CAL_ROOT/logs/locus_callability.log" || exit 1
echo "locus callability wall_s $(( $(date +%s) - t0 ))"; t0=$(date +%s)
[ -f "$CAL_ROOT/leaveout.ok" ] || { $X python "$CODE/calibration/graded_privacy/leaveout_status.py" --release-dir "$RELEASE" \
   --sample-manifest "$SAMPLE_MANIFEST" --samples-root "$SAMPLES_ROOT" --families "$CAL_ROOT/families.tsv" \
   --locus-callable "$CAL_ROOT/locus_callable.tsv.gz" --modes LOFO --threads 8 --samples "$CAL_ROOT/pairs.tsv" \
   --out "$CAL_ROOT/leaveout" 2> "$CAL_ROOT/logs/leaveout_status.log" && touch "$CAL_ROOT/leaveout.ok"; } || exit 1
echo "leaveout wall_s $(( $(date +%s) - t0 ))"; t0=$(date +%s)
# depth tasks: cohort sample alignment (from the registry sample sheets), largest-first so the long ones start early
python3 - "$SAMPLE_MANIFEST" "$SAMPLE_SHEETS" > "$CAL_ROOT/depth_tasks.txt" <<'PY'
import csv, os, sys
man, sheets = sys.argv[1:3]
rows = []
for r in csv.DictReader(open(man), delimiter="\t"):
    sh = {x["sample_id"]: x for x in csv.DictReader(open(f"{sheets}/{r['cohort']}.tsv"), delimiter="\t")}
    aln = sh[r["sample_id"]]["alignment"]
    rows.append((os.path.getsize(aln), r["cohort"], r["sample_id"], aln))
for _, c, s, a in sorted(rows, reverse=True):
    print(c, s, a)
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
