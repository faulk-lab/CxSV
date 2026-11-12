#!/bin/bash
# panel-v1.3 design check, step 4 (subset): re-make one panel genome's CxSV loci with a code
# tree's copy-number linkage + annotation + junction extraction, from the panel-v1.2-r2
# per-sample upstream outputs. Inputs are COPIED (never linked) into NEW; the release tree is
# only read.   subset_tree.sh CODE NEW KEY INDIVIDUAL [MODE]
#   MODE branch (default): cn_linkage.py with --pon-members (v1.3 D2 + item g)
#   MODE main:              cn_linkage.py of origin/main (fidelity check: must reproduce the release)
# CN segments: the release's (nested_min_z 6.0) unless SEG_DIR (D5 sweep output for another z)
# and PANEL_SEGMENTS are set.
set -euo pipefail
CODE=$1; NEW=$2; KEY=$3; IND=$4; MODE=${5:-branch}
REL=/mnt/sargasso1/CxSV/build/panel-v1.2-r2/results/panel/release
SRC=/mnt/sargasso1/CxSV/build/panel-v1.2-r2/results/panel/samples/$KEY
N=$NEW/$KEY
SIF=/mnt/sargasso1/CxSV/bundles/cxsv-bundle-20261109-84d380e08cc4/container/cxsv-pipeline.sif
RES=/mnt/sargasso1/CxSV/resources
PANEL_SEGMENTS=${PANEL_SEGMENTS:-$REL/cn/panel_cn_segments.tsv}
[ -s "$N/done.ok" ] && exit 0
mkdir -p $N/callers $N/cxsv $N/cn $N/junctions $N/logs
cp $SRC/merged_sv.vcf $N/
for c in sniffles cutesv svim; do cp $SRC/callers/${c}_filtered.vcf $N/callers/; done
for f in junction_clusters.bed junction_summary.tsv junction_only.vcf single_caller_filtered.tsv long_span_records.tsv; do
  cp $SRC/cxsv/$f $N/cxsv/; done
TAG=${KEY//\//__}
if [ -n "${SEG_DIR:-}" ]; then cp $SEG_DIR/$TAG.tsv $N/cn/segments.tsv; cp $SEG_DIR/$TAG.json $N/cn/stats.json
else cp $SRC/cn/segments.tsv $SRC/cn/stats.json $N/cn/; fi
X="singularity exec --containall -B /mnt/sargasso1 -B /home/minknow/scratch $SIF"
$X python /home/minknow/scratch/v13/design_check/tools/cn_params_of.py $CODE $N/cn_params.json > $N/junction_params.txt
read CD IR < $N/junction_params.txt
PM=""; [ "$MODE" = branch ] && PM="--pon-members $REL/cn/pon/pon_members.tsv"
$X python $CODE/workflow/scripts/cn_linkage.py --vcf $N/merged_sv.vcf --clusters $N/cxsv/junction_clusters.bed \
   --summary $N/cxsv/junction_summary.tsv --cxsv-vcf $N/cxsv/junction_only.vcf --single $N/cxsv/single_caller_filtered.tsv \
   --segments $N/cn/segments.tsv --panel-segments $PANEL_SEGMENTS $PM --individual "$IND" --stats $N/cn/stats.json \
   --fai /mnt/sargasso1/CxSV/refs/hs1.fa.fai --calling $N/cn_params.json --scripts $CODE/workflow/scripts --work $N/cn/work \
   --out-clusters $N/cxsv/clusters.bed --out-summary $N/cxsv/cxsv_summary.tsv --out-vcf $N/cxsv/cxsv_only.vcf \
   --out-candidates $N/cn/cn_candidates.tsv --out-loci $N/cn/cn_loci.tsv --out-cn-only $N/cn/cn_only.tsv 2> $N/logs/cn_linkage.log
rm -rf $N/cn/work
$X python $CODE/workflow/scripts/annotate_cxsv.py --summary $N/cxsv/cxsv_summary.tsv --vcf $N/cxsv/cxsv_only.vcf \
   --genes $RES/genes.final.bed --exons $RES/exons.filtered.bed --repeats $RES/repeats.final.bed --segdups $RES/segdup.final.bed \
   --out-table $N/cxsv/cxsv_master_table.tsv --out-bed $N/cxsv/cxsv_loci.bed > $N/logs/annotate_cxsv.log 2>&1
$X python $CODE/workflow/scripts/extract_junctions.py --vcf $N/merged_sv.vcf --clusters $N/cxsv/clusters.bed \
   --summary $N/cxsv/cxsv_summary.tsv --fai /mnt/sargasso1/CxSV/refs/hs1.fa.fai --collapse-distance $CD \
   --ins-len-ratio $IR --out $N/junctions/junctions.tsv 2> $N/logs/extract_junctions.log
date -Is > $N/done.ok
