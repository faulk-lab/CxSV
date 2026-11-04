#!/bin/bash
# Depth + catalogue inputs of the scorer for one panel sample re-made with the item-13
# code (calibration/graded_privacy/sample_depth_catalogues.sh, DEPTH_LOCI=nonpublic,
# LOFO only). Read depth is measured only for loci that the item-10 calibration did
# not measure: a locus whose id (Chr:Start-End) is unchanged has the same junctions
# (item 13 leaves every cluster that was not oversized as it was), so its segments
# and depth rows are copied from OLD_CAL. Usage: eval_depth.sh COHORT SAMPLE ALIGNMENT
set -uo pipefail
source "$(dirname "$0")/../graded_privacy/config.sh"
cohort=$1; sample=$2; aln=$3
S=$SAMPLES_ROOT/$cohort/$sample
O=$CAL_ROOT/persample/$cohort/$sample
OLD=$OLD_CAL/persample/$cohort/$sample
L=$CAL_ROOT/leaveout/$cohort/$sample
mkdir -p "$O"
[ -f "$O/depth_catalogues.ok" ] && exit 0
t0=$(date +%s)
{
  echo "host $(hostname) start $(date -Is)"
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="Panel_Status") k=i; print "Chr\tStart\tEnd"; next}
              $k!="public"{print $1"\t"$2"\t"$3}' "$L/status_LOFO.tsv" > "$O/need.tsv"
  cat "$OLD/depth.tsv" > "$O/old_depth.tsv"
  R2=${REUSE_CAL:+$REUSE_CAL/persample/$cohort/$sample/depth.tsv}
  [ -n "$R2" ] && [ -f "$R2" ] && awk -F'\t' 'NR==FNR{seen[$4]=1; print; next} !($4 in seen)' "$OLD/depth.tsv" "$R2" > "$O/old_depth.tsv"
  cut -f4 "$O/old_depth.tsv" | sort -u > "$O/old_ids.txt"
  awk -F'\t' 'NR==FNR{o[$1]=1;next} FNR==1{print;next} !(($1":"$2"-"$3) in o)' "$O/old_ids.txt" "$O/need.tsv" > "$O/new_loci.tsv"
  awk -F'\t' 'NR==FNR{if(FNR>1) n[$1":"$2"-"$3]=1;next} ($4 in n)' "$O/need.tsv" "$O/old_depth.tsv" > "$O/depth.reused.tsv"
  echo "need $(( $(wc -l < "$O/need.tsv") - 1 )) new $(( $(wc -l < "$O/new_loci.tsv") - 1 )) reused_rows $(wc -l < "$O/depth.reused.tsv")"
  : > "$O/depth.new.tsv"
  if [ "$(wc -l < "$O/new_loci.tsv")" -gt 1 ]; then
    $X python "$CODE/workflow/scripts/locus_segments.py" --junctions "$S/junctions/junctions.tsv" \
        --loci "$O/new_loci.tsv" --fai "$REF.fai" --collapse-distance 300 --composite-fraction 0.25 \
        --flank 20000 --gap 1000 --margin 200 --out "$O/segs.all.bed" || exit 1
    $X samtools view -H "$aln" | awk -F'\t' '$1=="@SQ"{for(i=2;i<=NF;i++) if($i~/^SN:/) print substr($i,4)}' > "$O/sq.txt"
    awk 'NR==FNR{ok[$1]=1;next} ($1 in ok)' "$O/sq.txt" "$O/segs.all.bed" > "$O/segs.bed"
    echo "segments $(wc -l < "$O/segs.bed") bp $(awk '{s+=$3-$2} END{print s+0}' "$O/segs.bed")"
    if [ -s "$O/segs.bed" ]; then
      REF_PATH=/nonexistent $X samtools bedcov --input-fmt-option required_fields=0x1FF -Q 20 -j \
          --reference "$REF" "$O/segs.bed" "$aln" > "$O/depth.new.tsv" || exit 1
      [ "$(wc -l < "$O/depth.new.tsv")" -eq "$(wc -l < "$O/segs.bed")" ] || { echo "depth line count mismatch"; exit 1; }
    fi
  fi
  cat "$O/depth.reused.tsv" "$O/depth.new.tsv" > "$O/depth.tsv"
  CATS=""
  for nv in $CATALOGUE_SPECS; do n=${nv%%:*}; v=${nv#*:}; d=$CATALOGUES/$n/$v
    CATS="$CATS --catalogue $n:$v:$d/$n.$v.junction_sites.tsv.gz:$d/$n.$v.manifest.tsv"; done
  $X python "$CODE/workflow/scripts/match_catalogues.py" --junctions "$S/junctions/junctions.tsv" \
      --loci "$S/cxsv/cxsv_master_table.tsv" $CATS --tolerance 500 --ins-len-ratio 0.5 --min-junction-fraction 0.75 \
      --out-loci "$O/cats.tsv" --out-junctions "$O/jcats.tsv" --out-check "$O/cat_check.tsv" || exit 1
  rm -f "$O/segs.all.bed" "$O/old_depth.tsv"
  echo "done $(date -Is) wall_s $(( $(date +%s) - t0 ))"
} > "$O/depth_catalogues.log" 2>&1 && touch "$O/depth_catalogues.ok" || { echo "FAILED $cohort $sample" >&2; exit 1; }
