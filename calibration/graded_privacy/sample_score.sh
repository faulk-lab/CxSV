#!/bin/bash
# Score one panel sample's own CxSV loci with the patient scorer
# (workflow/scripts/score_cxsv.py, unchanged) against the panel with the sample's
# individual (LOIO) or family (LOFO) removed. Inputs: leaveout_status.py
# (status_<MODE>.tsv, jstatus_<MODE>.tsv) and sample_depth_catalogues.sh
# (depth.tsv, cats.tsv). Output: scores_<MODE>.tsv (+ params.json).
# Usage: sample_score.sh COHORT SAMPLE SEX [MODES]
set -uo pipefail
source "$(dirname "$0")/config.sh"
cohort=$1; sample=$2; sex=$3; modes=${4:-LOIO LOFO}
S=$SAMPLES_ROOT/$cohort/$sample
O=$CAL_ROOT/persample/$cohort/$sample
L=$CAL_ROOT/leaveout/$cohort/$sample
tag=${SCORE_TAG:-}
[ -f "$O/score${tag}.ok" ] && exit 0
{
  echo "host $(hostname) start $(date -Is)"
  $X python "$CODE/calibration/graded_privacy/make_params.py" --release-json "$RELEASE/release.json" --sex "$sex" \
      ${WEIGHTS_OVERRIDE:+--weights-override "$WEIGHTS_OVERRIDE"} --out "$O/params${tag}.json" || exit 1
  for m in $modes; do
    t0=$(date +%s)
    $X python "$CODE/workflow/scripts/score_cxsv.py" --status "$L/status_$m.tsv" --catalogue-hits "$O/cats.tsv" \
        --catalogue-check "$O/cat_check.tsv" --junctions "$S/junctions/junctions.tsv" --junction-status "$L/jstatus_$m.tsv" \
        --vcf "$S/cxsv/cxsv_only.vcf" --caller-vcf "Sniffles2=$S/callers/sniffles_filtered.vcf" \
        --caller-vcf "cuteSV=$S/callers/cutesv_filtered.vcf" --caller-vcf "SVIM=$S/callers/svim_filtered.vcf" \
        --depth "$O/depth.tsv" --params "$O/params${tag}.json" --genes "$RES/genes.final.bed" \
        --exons "$RES/exons.filtered.bed" --repeats "$RES/repeats.final.bed" --segdups "$RES/segdup.final.bed" \
        --fai "$REF.fai" --sample "$sample" --out-tsv "$O/scores${tag}_$m.tsv.tmp" \
        --out-features "$O/features.tmp.$m" --out-report "$O/report.tmp.$m" || exit 1
    mv "$O/scores${tag}_$m.tsv.tmp" "$O/scores${tag}_$m.tsv"; rm -f "$O/features.tmp.$m" "$O/report.tmp.$m"
    echo "$m wall_s $(( $(date +%s) - t0 ))"
  done
  echo "done $(date -Is)"
} > "$O/score${tag}.log" 2>&1 && touch "$O/score${tag}.ok" || { echo "FAILED $cohort $sample" >&2; exit 1; }
