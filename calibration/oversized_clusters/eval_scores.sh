#!/bin/bash
# Roadmap item 13, step 4b: healthy-genome Tier-1/Tier-2 rates before and after the
# item-13 clustering on a subset of panel genomes, with the item-10 LOFO approach
# (calibration/graded_privacy: each genome's own loci scored by the unchanged patient
# scorer against the released panel-v1.0 junction sites with its family removed).
# The released junction sites do not depend on clustering (build_panel_sites.py reads
# every junction, whatever its locus), so they are valid for both arms; region-overlap
# carriers (informational, not scored) still come from the panel-v1.0 loci.
#
#   eval_scores.sh SUBSET.tsv EVALROOT NEW_SAMPLES_ROOT [PAR]
#
# before: SAMPLES_ROOT = panel-v1.0 samples; leave-out statuses, depth and catalogue
#         hits of the item-10 calibration (copied), re-scored with this checkout
# after:  SAMPLES_ROOT = NEW_SAMPLES_ROOT (downstream_rerun.sh); locus callability and
#         LOFO statuses recomputed; depth only for loci item 10 did not measure
#         (AFTER_NAME names the arm's directory; REUSE_CAL, if set, is a second
#         calibration root whose depth rows are re-used for loci with the same id)
set -uo pipefail
SUB=$1; E=$2; NEWROOT=$3; PAR=${4:-8}
HERE="$(cd "$(dirname "$0")" && pwd)"
export OLD_CAL=/mnt/sargasso1/CxSV/calibration/item10
export SIF=${SIF:-/mnt/sargasso1/CxSV/bundles/cxsv-bundle-20261103-33c075a22071/container/cxsv-pipeline.sif}
BUILD=/mnt/sargasso1/CxSV/build/panel-v1.0
MAN=$BUILD/results/panel/panel/sample_manifest.tsv
mkdir -p "$E"
awk -F'\t' 'NR==FNR{if(FNR>1) w[$1"\t"$2]=1; next} FNR==1 || (($1"\t"$2) in w)' "$SUB" "$MAN" > "$E/subset_manifest.tsv"
tail -n +2 "$SUB" | cut -f1,2 > "$E/subset_pairs.tsv"
sexof() { awk -F'\t' -v i="$1" 'NR==1{for(k=1;k<=NF;k++) if($k=="sex") c=k; next} $1==i{print $c}' "$OLD_CAL/families.tsv"; }
alnof() { awk -F'\t' -v c="$1" -v s="$2" '$1==c && $2==s {print $NF}' "$SUB"; }

# ---- before ----
export CAL_ROOT=$E/cal_before SAMPLES_ROOT=$BUILD/results/panel/samples SAMPLE_MANIFEST=$E/subset_manifest.tsv
mkdir -p "$CAL_ROOT/persample"; ln -sfn "$OLD_CAL/leaveout" "$CAL_ROOT/leaveout"
while IFS=$'\t' read -r c s; do
  mkdir -p "$CAL_ROOT/persample/$c/$s"
  cp "$OLD_CAL/persample/$c/$s/"{depth.tsv,cats.tsv,cat_check.tsv} "$CAL_ROOT/persample/$c/$s/"
done < "$E/subset_pairs.tsv"
while IFS=$'\t' read -r c s; do
  i=$(awk -F'\t' -v c="$c" -v s="$s" '$1==c && $2==s{print $3}' "$E/subset_manifest.tsv")
  echo "$c $s $(sexof "$i")"
done < "$E/subset_pairs.tsv" > "$E/score_tasks.txt"
t0=$(date +%s)
xargs -P "$PAR" -L 1 bash -c 'bash '"$HERE"'/../graded_privacy/sample_score.sh "$0" "$1" "$2" LOFO' < "$E/score_tasks.txt"
echo "before scored: $(ls "$CAL_ROOT"/persample/*/*/score.ok | wc -l) wall_s $(( $(date +%s) - t0 ))"

# ---- after ----
export CAL_ROOT=$E/${AFTER_NAME:-cal_after} SAMPLES_ROOT=$NEWROOT
mkdir -p "$CAL_ROOT/logs"
source "$HERE/../graded_privacy/config.sh"
flank=$(python3 -c "import json;print(json.load(open('$RELEASE/release.json'))['callability']['flank_bp'])")
t0=$(date +%s)
$X python "$CODE/calibration/graded_privacy/locus_callability.py" --sample-manifest "$E/subset_manifest.tsv" \
   --samples-root "$SAMPLES_ROOT" --callability-dir "$RELEASE/callability" --flank "$flank" \
   --families "$OLD_CAL/families.tsv" --threads 16 --out "$CAL_ROOT/locus_callable.tsv.gz" 2> "$CAL_ROOT/logs/locus_callability.log" || exit 1
echo "after locus callability wall_s $(( $(date +%s) - t0 ))"; t0=$(date +%s)
$X python "$CODE/calibration/graded_privacy/leaveout_status.py" --release-dir "$RELEASE" \
   --sample-manifest "$MAN" --samples-root "$SAMPLES_ROOT" --families "$OLD_CAL/families.tsv" \
   --locus-callable "$CAL_ROOT/locus_callable.tsv.gz" --modes LOFO --threads 8 --samples "$E/subset_pairs.tsv" \
   --out "$CAL_ROOT/leaveout" 2> "$CAL_ROOT/logs/leaveout_status.log" || exit 1
echo "after leaveout wall_s $(( $(date +%s) - t0 ))"; t0=$(date +%s)
while IFS=$'\t' read -r c s; do echo "$c $s $(alnof "$c" "$s")"; done < "$E/subset_pairs.tsv" > "$E/depth_tasks.txt"
xargs -P "$PAR" -L 1 bash -c 'bash '"$HERE"'/eval_depth.sh "$0" "$1" "/mnt/sargasso1/CxSV/datasets/$2"' < "$E/depth_tasks.txt"
echo "after depth: $(ls "$CAL_ROOT"/persample/*/*/depth_catalogues.ok | wc -l) wall_s $(( $(date +%s) - t0 ))"; t0=$(date +%s)
xargs -P "$PAR" -L 1 bash -c 'bash '"$HERE"'/../graded_privacy/sample_score.sh "$0" "$1" "$2" LOFO' < "$E/score_tasks.txt"
echo "after scored: $(ls "$CAL_ROOT"/persample/*/*/score.ok | wc -l) wall_s $(( $(date +%s) - t0 ))"
