#!/bin/bash
# Run the named rows of rescore_tasks.tsv (name src label sample alignment sheet code) with
# rerun_patient_downstream.sh, two at a time, in $HOME/v13dc_rescore/<name> on this node.
#   rescore_tasks.sh NAME [NAME ...]
set -uo pipefail
D5=/mnt/sargasso1/CxSV/calibration/panel_v13_design/d5
OUT=/mnt/sargasso1/CxSV/calibration/panel_v13_design/rescore
BUN=/mnt/sargasso1/CxSV/bundles/cxsv-bundle-20261109-84d380e08cc4
mkdir -p $OUT "$HOME/v13dc_rescore"
one() {
  read -r name src label sample aln sheet code < <(awk -v n="$1" '$1==n' $OUT/rescore_tasks.tsv)
  W=$HOME/v13dc_rescore/$name; rm -rf "$W"
  bash $D5/tools/rerun_patient_downstream.sh $BUN $src $W $code $label $sample $aln $sheet 12 > $OUT/$name.log 2>&1
  rc=$?
  mkdir -p $OUT/$name
  S=$W/run/results/patient/samples/$label/$sample
  cp -r $S/scoring $S/cn $S/private $OUT/$name/ 2>/dev/null; rm -rf $OUT/$name/cn/depth.regions.bed.gz* $OUT/$name/cn/work
  cp $W/run/results/patient/panel_compatibility.txt $W/dryrun.txt $W/rerun.log $OUT/$name/ 2>/dev/null
  echo "$name rc=$rc host=$(hostname) $(date -Is)" | tee -a $OUT/status.txt
}
export -f one; export OUT D5 BUN
printf "%s\n" "$@" | xargs -P 2 -I{} bash -c 'one {}'
