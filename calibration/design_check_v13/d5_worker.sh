#!/bin/bash
# panel-v1.3 design check D5: one genome, nested_min_z sweep (cn_segment.py --set + cn_linkage.py
# of the integration branch) against the panel-v1.2-r2 release PoN / panel segments / members.
#   d5_worker.sh KEY INDIVIDUAL STRATUM SAMPLE_DIR "Z1 Z2 ..."
set -euo pipefail
KEY=$1; IND=$2; STRATUM=$3; S=$4; ZS=$5
DC=/home/minknow/scratch/v13/design_check; CODE=$DC/repo
REL=/mnt/sargasso1/CxSV/build/panel-v1.2-r2/results/panel/release
SIF=/mnt/sargasso1/CxSV/bundles/cxsv-bundle-20261109-84d380e08cc4/container/cxsv-pipeline.sif
X="singularity exec --containall -B /mnt/sargasso1 -B /home/minknow/scratch $SIF"
TAG=${KEY//\//__}
for Z in $ZS; do
  O=$DC/d5/z$Z/$TAG; mkdir -p $O/work
  [ -s $O/done.ok ] && continue
  t0=$(date +%s.%N)
  $X python $CODE/workflow/scripts/cn_segment.py --depth $S/cn/depth.regions.bed.gz --pon $REL/cn/pon/$STRATUM.npz \
     --params $DC/d5/cn_params_branch.json --stratum $STRATUM --default-stratum $STRATUM \
     --set segmentation.nested_min_z=$Z --out $O/segments.tsv --stats $O/stats.json 2> $O/cn_segment.log
  t1=$(date +%s.%N)
  $X python $CODE/workflow/scripts/cn_linkage.py --vcf $S/merged_sv.vcf --clusters $S/cxsv/junction_clusters.bed \
     --summary $S/cxsv/junction_summary.tsv --cxsv-vcf $S/cxsv/junction_only.vcf --single $S/cxsv/single_caller_filtered.tsv \
     --segments $O/segments.tsv --panel-segments $REL/cn/panel_cn_segments.tsv --pon-members $REL/cn/pon/pon_members.tsv \
     --individual "$IND" --stats $O/stats.json --fai /mnt/sargasso1/CxSV/refs/hs1.fa.fai --calling $DC/d5/cn_params_branch.json \
     --scripts $CODE/workflow/scripts --work $O/work --out-clusters $O/clusters.bed --out-summary $O/cxsv_summary.tsv \
     --out-vcf $O/cxsv_only.vcf --out-candidates $O/cn_candidates.tsv --out-loci $O/cn_loci.tsv --out-cn-only $O/cn_only.tsv \
     2> $O/cn_linkage.log
  t2=$(date +%s.%N)
  rm -rf $O/work $O/cxsv_only.vcf
  printf "%s\t%s\t%.1f\t%.1f\n" "$KEY" "$Z" "$(echo "$t1-$t0"|bc)" "$(echo "$t2-$t1"|bc)" > $O/done.ok
done
