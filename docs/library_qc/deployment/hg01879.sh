#!/usr/bin/env bash
# Library QC deployment, step 2d: flagged panel library HG01879 (ONT 1kGP-1019; HG04176 had no Tier-1 candidate) end to end through the new bundle's patient workflow.
# Run A: the deployment configuration (bundled panel-v1.1; this library is itself a panel member).
# Run B: same caller outputs (copied workdir), panel_release = the library-QC exclusion test release, which
#        leaves out all 24 flagged libraries (this ONT library included; HG01879's PacBio library stays in),
#        so the noise junctions of this library are private, as for an unseen patient.
set -euo pipefail
B=/mnt/sargasso1/CxSV/bundles/cxsv-bundle-20261105-95126e14f036
D=/home/minknow/scratch/libqc2/hg01879; XR=/mnt/sargasso1/CxSV/scratch/libqc_exclusion/results/panel/release
S=/mnt/sargasso1/CxSV/datasets/1kGP_ONT_1019_CHM13/HG01879.ERR15140163.hs1.minimap2.cram
mkdir -p $D; cp $S $S.crai $D/; [ "$(stat -c %s $D/HG01879.ERR15140163.hs1.minimap2.cram)" = "$(stat -c %s $S)" ]; cmp $S $D/HG01879.ERR15140163.hs1.minimap2.cram && echo "staged (identical)"
touch -r $D/HG01879.ERR15140163.hs1.minimap2.cram $D/HG01879.ERR15140163.hs1.minimap2.cram.crai
printf "sample_id\talignment\tindividual_id\tplatform\nHG01879\t$D/HG01879.ERR15140163.hs1.minimap2.cram\tHG01879\tONT\n" > $D/hg01879.tsv
cat > $D/runA.yaml <<YAML
patient_runs:
  HG01879_flagged: {samples: "$D/hg01879.tsv"}
threads: {samtools: 4, sv_callers: 12}
YAML
t0=$(date +%s)
$B/run_cxsv.sh patient --runtime singularity --workdir $D/runA --config $D/runA.yaml --bind $D --cores 12 > $D/runA.log 2>&1 && rc=0 || rc=$?
echo "run A exit $rc wall $(( $(date +%s)-t0 )) s"; [ $rc = 0 ] || { tail -30 $D/runA.log; exit $rc; }
rm -rf $D/runB; cp -a $D/runA $D/runB
cat > $D/runB.yaml <<YAML
patient_runs:
  HG01879_flagged: {samples: "$D/hg01879.tsv"}
threads: {samtools: 4, sv_callers: 12}
panel_release: "$XR"
YAML
t0=$(date +%s)
$B/run_cxsv.sh patient --runtime singularity --workdir $D/runB --config $D/runB.yaml --bind $D --bind $XR --cores 12 > $D/runB.log 2>&1 && rc=0 || rc=$?
echo "run B exit $rc wall $(( $(date +%s)-t0 )) s"; [ $rc = 0 ] || { tail -30 $D/runB.log; exit $rc; }
for r in runA runB; do f=$(find $D/$r/results -name cxsv_score_report.md | head -1); echo "== $r $f"; sed -n 1,30p "$f"; done
