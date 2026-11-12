#!/bin/bash
# D5 pass 1: segment one library at every swept nested_min_z (integration-branch cn_segment.py,
# panel-v1.2-r2 release PoN). Output NAS $D5/seg/z<Z>/<cohort>__<sample>.{tsv,json}
#   d5_seg.sh KEY STRATUM DEPTH
set -uo pipefail
KEY=$1; STRATUM=$2; DEPTH=$3
D5=/mnt/sargasso1/CxSV/calibration/panel_v13_design/d5; CODE=$D5/code
REL=/mnt/sargasso1/CxSV/build/panel-v1.2-r2/results/panel/release
SIF=/mnt/sargasso1/CxSV/bundles/cxsv-bundle-20261109-84d380e08cc4/container/cxsv-pipeline.sif
TAG=${KEY//\//__}
LOCAL=$(mktemp -d /tmp/d5seg.XXXXXX); trap 'rm -rf $LOCAL' EXIT
cp $DEPTH $LOCAL/depth.regions.bed.gz
for Z in ${ZS:-5.0 5.25 5.5 5.75 6.0}; do
  O=$D5/seg/z$Z; mkdir -p $O
  [ -s $O/$TAG.json ] && continue
  singularity exec --containall -B /mnt/sargasso1 -B $LOCAL $SIF python $CODE/workflow/scripts/cn_segment.py \
     --depth $LOCAL/depth.regions.bed.gz --pon $REL/cn/pon/$STRATUM.npz --params $D5/cn_params_branch.json \
     --stratum $STRATUM --default-stratum $STRATUM --set segmentation.nested_min_z=$Z \
     --out $LOCAL/s.tsv --stats $LOCAL/s.json 2> $LOCAL/s.log \
   && cp $LOCAL/s.tsv $O/$TAG.tsv && cp $LOCAL/s.json $O/$TAG.json || { echo "FAILED $KEY z$Z"; cat $LOCAL/s.log | tail -3; }
done
