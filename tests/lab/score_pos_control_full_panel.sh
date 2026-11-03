#!/bin/bash
# Re-score the positive control's EXISTING call set against a full panel release,
# with two scorer checkouts (e.g. main = before, a feature branch = after), and
# score each against the truth set. No calling is rerun: the calls, junction
# table, caller VCFs and locus read depth come from a finished regression
# workdir (tests/lab/pos_control_regression.sh). Per checkout:
#   classify_patient_cxsv.py (full release, verified, callability) -> private/
#   match_catalogues.py                                             -> private/
#   score_cxsv.py (that checkout's config/scoring.yaml)            -> scoring/
#   score_pos_control.py                                            -> <tag>.score.*
# Usage: score_pos_control_full_panel.sh REGRESSION_WORKDIR RELEASE_DIR OUTDIR \
#            TAG1=CHECKOUT1 [TAG2=CHECKOUT2 ...]
# Runs inside the pinned SIF when SIF is set (python + samtools of the bundle).
set -euo pipefail
REG=$1; RELEASE=$2; OUT=$3; shift 3
S0=$REG/results/patient/samples/pos_control/pos_control
PC=${PC:-/mnt/sargasso1/CxSV/datasets/pos_control}
REF=${REF:-/mnt/sargasso1/CxSV/refs/hs1.fa}
RES=${RES:-/mnt/sargasso1/CxSV/resources}
CATALOGUES=${CATALOGUES:-/mnt/sargasso1/CxSV/datasets/external_catalogues}
X=${SIF:+singularity exec -B /mnt/sargasso1 -B $HOME $SIF}
mkdir -p "$OUT"
CATS=""
for nv in CoLoRSdb:v1.2.0 HGSVC3:HGSVC2024v1.0 HPRC:release2-mc-v2.1 gnomAD-SV:v4.1; do n=${nv%%:*}; v=${nv#*:}; d=$CATALOGUES/$n/$v
  CATS="$CATS --catalogue $n:$v:$d/$n.$v.junction_sites.tsv.gz:$d/$n.$v.manifest.tsv"; done
for spec in "$@"; do
  tag=${spec%%=*}; code=${spec#*=}
  D=$OUT/$tag; mkdir -p "$D/private" "$D/scoring" "$D/cxsv" "$D/junctions" "$D/callers"
  cp -f "$S0/merged_sv.vcf" "$D/merged_sv.vcf"   # copies: the NAS (CIFS) has no symlinks
  for f in "$S0"/cxsv/*; do cp -f "$f" "$D/cxsv/"; done
  cp -f "$S0/junctions/junctions.tsv" "$D/junctions/junctions.tsv"
  for c in sniffles cutesv svim; do cp -f "$S0/callers/${c}_filtered.vcf" "$D/callers/${c}_filtered.vcf"; done
  cp -f "$S0/scoring/locus_depth.tsv" "$D/scoring/locus_depth.tsv"
  $X python "$code/workflow/scripts/classify_patient_cxsv.py" --master-table "$D/cxsv/cxsv_master_table.tsv" \
     --junctions "$D/junctions/junctions.tsv" --panel-carriers "$RELEASE/panel_junction_carriers.tsv" \
     --panel-sites "$RELEASE/panel_junction_sites.tsv" --panel-individuals "$RELEASE/panel_individuals.tsv" \
     --panel-calls "$RELEASE/panel_cxsv_calls.tsv" --tolerance 500 --ins-len-ratio 0.5 --min-junction-fraction 0.75 \
     --min-carriers 1 --window 1000 --release-dir "$RELEASE" --out-all "$D/private/cxsv_panel_status.tsv" \
     --out-private "$D/private/private_cxsv_master_table.tsv" --out-bed "$D/private/private_cxsv_loci.bed" \
     --out-junctions "$D/private/junction_panel_status.tsv" 2> "$D/classify.log"
  $X python "$code/workflow/scripts/match_catalogues.py" --junctions "$D/junctions/junctions.tsv" \
     --loci "$D/private/cxsv_panel_status.tsv" $CATS --tolerance 500 --ins-len-ratio 0.5 --min-junction-fraction 0.75 \
     --out-loci "$D/private/cxsv_catalogue_hits.tsv" --out-junctions "$D/private/junction_catalogue_hits.tsv" \
     --out-check "$D/private/catalogue_check.tsv" 2> "$D/catalogues.log"
  $X python "$(dirname "$0")/../../calibration/graded_privacy/make_params.py" --release-json "$RELEASE/release.json" \
     --config-dir "$code/config" --patient-from "$S0/scoring/scoring_params.json" --out "$D/scoring/scoring_params.json"
  $X python "$code/workflow/scripts/score_cxsv.py" --status "$D/private/cxsv_panel_status.tsv" \
     --catalogue-hits "$D/private/cxsv_catalogue_hits.tsv" --catalogue-check "$D/private/catalogue_check.tsv" \
     --junctions "$D/junctions/junctions.tsv" --junction-status "$D/private/junction_panel_status.tsv" \
     --vcf "$D/cxsv/cxsv_only.vcf" --caller-vcf "Sniffles2=$D/callers/sniffles_filtered.vcf" \
     --caller-vcf "cuteSV=$D/callers/cutesv_filtered.vcf" --caller-vcf "SVIM=$D/callers/svim_filtered.vcf" \
     --depth "$D/scoring/locus_depth.tsv" --params "$D/scoring/scoring_params.json" --genes "$RES/genes.final.bed" \
     --exons "$RES/exons.filtered.bed" --repeats "$RES/repeats.final.bed" --segdups "$RES/segdup.final.bed" \
     --fai "$REF.fai" --sample pos_control --out-tsv "$D/scoring/cxsv_scores.tsv" \
     --out-features "$D/scoring/cxsv_score_features.tsv" --out-report "$D/scoring/cxsv_score_report.md" 2> "$D/score.log"
  $X python "$code/tests/positive_control/score_pos_control.py" "$PC/release/truth" "$D" \
     "$code/tests/positive_control/events.yaml" "$OUT/$tag.score" --extra-truth "$PC/addons/E13/truth" \
     --known-failures "$code/tests/positive_control/known_failures.tsv" > "$OUT/$tag.score.txt" 2>&1 || true
  tail -3 "$OUT/$tag.score.txt"
done
