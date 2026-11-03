#!/bin/bash
# Check the leave-out shortcut against a TRUE panel rebuild for one family.
#   1. build_panel_sites.py (the release's own command and parameters) on the
#      sample manifest WITHOUT the family's samples
#   2. classify_patient_cxsv.py (unchanged) on each family sample against the
#      rebuilt carriers / individuals, region calls without the family and a
#      callability directory without the family's masks
#   3. score_cxsv.py (unchanged) with that status, same depth and catalogue hits
#   4. compare_rebuild.py: per locus, shortcut (leaveout LOFO) vs rebuild
# Usage: validate_rebuild.sh FAMILY_ID
set -uo pipefail
here=$(cd "$(dirname "$0")" && pwd)
source "$here/config.sh"
fam=$1
D=$CAL_ROOT/rebuild/$fam
# the callability directory is a set of symlinks to the release masks: on node-local
# disk, because the NAS (CIFS) does not support symlinks
CD=${LOCAL_SCRATCH:-$HOME/scratch}/item10_rebuild/$fam/callability
mkdir -p "$CD" "$D/panel"
members=$(awk -F'\t' -v f="$fam" '$2==f{print $1}' "$CAL_ROOT/families.tsv" | sort)
echo "family $fam: $(echo $members | tr ' ' ',')"
echo "$members" > "$D/members.txt"
awk -F'\t' 'NR==FNR{m[$1]=1;next} FNR==1 || !($3 in m)' "$D/members.txt" "$SAMPLE_MANIFEST" > "$D/panel/sample_manifest.tsv"
awk -F'\t' 'NR==FNR{m[$1]=1;next} FNR>1 && ($3 in m){print $1"\t"$2"\t"$3}' "$D/members.txt" "$SAMPLE_MANIFEST" > "$D/family_samples.tsv"
echo "manifest rows $(($(wc -l < "$D/panel/sample_manifest.tsv")-1)); family samples $(wc -l < "$D/family_samples.tsv")"
rid=$(python3 -c "import json;print(json.load(open('$RELEASE/release.json'))['release_id'])")
if [ ! -s "$D/panel/panel_junction_carriers.tsv" ]; then
  /usr/bin/time -v $X python "$CODE/workflow/scripts/build_panel_sites.py" --manifest "$D/panel/sample_manifest.tsv" \
     --individuals /mnt/sargasso1/CxSV/datasets/registry/individuals.tsv --fai "$REF.fai" --tolerance 500 --ins-len-ratio 0.5 \
     --min-svlen 50 --source cxsv_panel --version "$rid" --sites "$D/panel/panel_junction_sites.tsv.tmp" \
     --carriers "$D/panel/panel_junction_carriers.tsv.tmp" --panel-individuals "$D/panel/panel_individuals.tsv" \
     2> "$D/build_panel_sites.log" || { tail "$D/build_panel_sites.log"; exit 1; }
  mv "$D/panel/panel_junction_sites.tsv.tmp" "$D/panel/panel_junction_sites.tsv"
  mv "$D/panel/panel_junction_carriers.tsv.tmp" "$D/panel/panel_junction_carriers.tsv"
fi
grep -E "Elapsed|Maximum resident" "$D/build_panel_sites.log"
# callability without the family
awk -F'\t' 'NR==FNR{m[$1]=1;next} FNR==1 || !($4 in m)' "$D/members.txt" "$RELEASE/callability/samples.tsv" > "$CD/samples.tsv"
tail -n +2 "$CD/samples.tsv" | cut -f1 | while read -r m; do ln -sf "$RELEASE/callability/$m" "$CD/$m"; done
awk -F'\t' 'NR==FNR{m[$1]=1;next} FNR==1{for(i=1;i<=NF;i++) if($i=="Individual_ID") k=i; print; next} !($k in m)' \
    "$D/members.txt" "$RELEASE/panel_cxsv_calls.tsv" > "$D/panel/panel_cxsv_calls.tsv"
flank=$(python3 -c "import json;print(json.load(open('$RELEASE/release.json'))['callability']['flank_bp'])")
while IFS=$'\t' read -r cohort sample ind; do
  S=$SAMPLES_ROOT/$cohort/$sample; O=$D/$cohort/$sample; mkdir -p "$O"
  $X python "$CODE/workflow/scripts/classify_patient_cxsv.py" --master-table "$S/cxsv/cxsv_master_table.tsv" \
     --junctions "$S/junctions/junctions.tsv" --panel-carriers "$D/panel/panel_junction_carriers.tsv" \
     --panel-sites "$D/panel/panel_junction_sites.tsv" --panel-individuals "$D/panel/panel_individuals.tsv" \
     --panel-calls "$D/panel/panel_cxsv_calls.tsv" --tolerance 500 --ins-len-ratio 0.5 --min-junction-fraction 0.75 \
     --min-carriers 1 --window 1000 --callability-dir "$CD" --flank "$flank" \
     --out-all "$O/status.tsv" --out-private "$O/priv.tsv" --out-bed "$O/priv.bed" --out-junctions "$O/jstatus.tsv" \
     2> "$O/classify.log" || { cat "$O/classify.log"; exit 1; }
  P=$CAL_ROOT/persample/$cohort/$sample
  $X python "$CODE/workflow/scripts/score_cxsv.py" --status "$O/status.tsv" --catalogue-hits "$P/cats.tsv" \
     --catalogue-check "$P/cat_check.tsv" --junctions "$S/junctions/junctions.tsv" --junction-status "$O/jstatus.tsv" \
     --vcf "$S/cxsv/cxsv_only.vcf" --caller-vcf "Sniffles2=$S/callers/sniffles_filtered.vcf" \
     --caller-vcf "cuteSV=$S/callers/cutesv_filtered.vcf" --caller-vcf "SVIM=$S/callers/svim_filtered.vcf" \
     --depth "$P/depth.tsv" --params "$P/params.json" --genes "$RES/genes.final.bed" --exons "$RES/exons.filtered.bed" \
     --repeats "$RES/repeats.final.bed" --segdups "$RES/segdup.final.bed" --fai "$REF.fai" --sample "$sample" \
     --out-tsv "$O/scores.tsv" --out-features "$O/features.tsv" --out-report "$O/report.md" 2> "$O/score.log" \
     || { cat "$O/score.log"; exit 1; }
  $X python "$here/compare_rebuild.py" --rebuild-status "$O/status.tsv" --rebuild-scores "$O/scores.tsv" \
     --shortcut-status "$CAL_ROOT/leaveout/$cohort/$sample/status_LOFO.tsv" --shortcut-scores "$P/scores_LOFO.tsv" \
     --label "$fam/$cohort/$sample" --out "$O/compare.tsv"
done < "$D/family_samples.tsv"
