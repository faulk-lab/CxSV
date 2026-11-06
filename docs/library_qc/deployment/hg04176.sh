#!/usr/bin/env bash
# Library QC deployment, step 2c: flagged panel genome HG04176 end to end through the new bundle's patient workflow.
# Run A: the deployment configuration (bundled panel-v1.1; HG04176 is itself a panel member).
# Run B: same caller outputs (copied workdir), panel_release = the library-QC exclusion test release, which
#        leaves out all 24 flagged libraries (HG04176 included), so HG04176 is scored as an unseen patient.
set -euo pipefail
B=/mnt/sargasso1/CxSV/bundles/cxsv-bundle-20261105-95126e14f036
D=/home/minknow/scratch/libqc2/hg04176; XR=/mnt/sargasso1/CxSV/scratch/libqc_exclusion/results/panel/release
cat > $D/runA.yaml <<YAML
patient_runs:
  HG04176_flagged: {samples: "$D/hg04176.tsv"}
threads: {samtools: 4, sv_callers: 8}
YAML
t0=$(date +%s)
$B/run_cxsv.sh patient --runtime singularity --workdir $D/runA --config $D/runA.yaml --bind $D --cores 8 > $D/runA.log 2>&1 && rc=0 || rc=$?
echo "run A exit $rc wall $(( $(date +%s)-t0 )) s"; [ $rc = 0 ] || { tail -30 $D/runA.log; exit $rc; }
rm -rf $D/runB; cp -a $D/runA $D/runB
cat > $D/runB.yaml <<YAML
patient_runs:
  HG04176_flagged: {samples: "$D/hg04176.tsv"}
threads: {samtools: 4, sv_callers: 8}
panel_release: "$XR"
YAML
t0=$(date +%s)
$B/run_cxsv.sh patient --runtime singularity --workdir $D/runB --config $D/runB.yaml --bind $D --bind $XR --cores 8 > $D/runB.log 2>&1 && rc=0 || rc=$?
echo "run B exit $rc wall $(( $(date +%s)-t0 )) s"; [ $rc = 0 ] || { tail -30 $D/runB.log; exit $rc; }
for r in runA runB; do f=$(find $D/$r/results -name cxsv_score_report.md | head -1); echo "== $r $f"; sed -n 1,30p "$f"; done
