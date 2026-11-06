#!/usr/bin/env bash
# Library QC deployment, step 2b: full-genome positive control (E01-E13 + L01-L04 add-ons) from the new bundle.
# Reuses the byte-identical spliced input of the tagged-bundle run (same release CRAM, same add-ons);
# the release CRAM md5 was verified by that run, so SKIP_MD5=1.
set -euo pipefail
B=/mnt/sargasso1/CxSV/bundles/cxsv-bundle-20261105-95126e14f036
W=/home/minknow/scratch/libqc2/deploy/pc_full_L; BASE=/home/minknow/scratch/item11r2_full
A=/mnt/sargasso1/CxSV/datasets/pos_control/addons
mkdir -p $W/input; T=E13+L01+L02+L03+L04
for f in pos_control.$T.cram pos_control.$T.cram.crai pos_control.$T.cram.splice.ok $T.replace.bed $T.dropped_names.txt $T.dropped_records.bam; do cp -a $BASE/input/$f $W/input/; done
cmp $BASE/input/pos_control.$T.cram $W/input/pos_control.$T.cram && echo "spliced input copied (identical)"
t0=$(date +%s)
SKIP_MD5=1 ADDONS="$A/E13 $A/L01 $A/L02 $A/L03 $A/L04" EVENTS=$A/events_E01-E13+L01-L04.yaml \
  bash $B/repo/tests/lab/pos_control_full_genome.sh $B $W 10 && rc=0 || rc=$?
echo "pc exit $rc wall $(( $(date +%s)-t0 ))"; exit $rc
