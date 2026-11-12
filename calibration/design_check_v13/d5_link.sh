#!/bin/bash
# D5 pass 2: cn_linkage.py (integration branch) of one genome at one z against the panel table of
# that z (d5_panel_table.py) and the release pon_members.tsv. Output $D5/link/z<Z>/<tag>/
#   d5_link.sh KEY INDIVIDUAL SAMPLE_DIR Z
set -uo pipefail
KEY=$1; IND=$2; S=$3; Z=$4
D5=/mnt/sargasso1/CxSV/calibration/panel_v13_design/d5; CODE=$D5/code
REL=/mnt/sargasso1/CxSV/build/panel-v1.2-r2/results/panel/release
SIF=/mnt/sargasso1/CxSV/bundles/cxsv-bundle-20261109-84d380e08cc4/container/cxsv-pipeline.sif
TAG=${KEY//\//__}; O=$D5/link/z$Z/$TAG; mkdir -p $O
[ -s $O/cn_loci.tsv ] && exit 0
LOCAL=$(mktemp -d /tmp/d5link.XXXXXX); trap 'rm -rf $LOCAL' EXIT
singularity exec --containall -B /mnt/sargasso1 -B /home/minknow/scratch -B $LOCAL $SIF python $CODE/workflow/scripts/cn_linkage.py \
   --vcf $S/merged_sv.vcf --clusters $S/cxsv/junction_clusters.bed --summary $S/cxsv/junction_summary.tsv \
   --cxsv-vcf $S/cxsv/junction_only.vcf --single $S/cxsv/single_caller_filtered.tsv \
   --segments $D5/seg/z$Z/$TAG.tsv --panel-segments $D5/panel_tables/panel_cn_segments.z$Z.tsv \
   --pon-members $REL/cn/pon/pon_members.tsv --individual "$IND" --stats $D5/seg/z$Z/$TAG.json \
   --fai /mnt/sargasso1/CxSV/refs/hs1.fa.fai --calling $D5/cn_params_branch.json --scripts $CODE/workflow/scripts \
   --work $LOCAL/work --out-clusters $LOCAL/clusters.bed --out-summary $LOCAL/cxsv_summary.tsv --out-vcf $LOCAL/cxsv_only.vcf \
   --out-candidates $LOCAL/cn_candidates.tsv --out-loci $LOCAL/cn_loci.tsv --out-cn-only $LOCAL/cn_only.tsv 2> $LOCAL/cn_linkage.log
rc=$?
cp $LOCAL/cn_linkage.log $O/
[ $rc -eq 0 ] && cp $LOCAL/cn_candidates.tsv $LOCAL/cn_only.tsv $O/ && cp $LOCAL/cn_loci.tsv $O/ || echo "FAILED link $KEY z$Z"
