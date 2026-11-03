#!/bin/bash
# Per panel sample, the two scorer inputs that do not depend on leave-out:
#   depth.tsv  read depth of every CxSV locus's segments and flanks (rule locus_depth:
#              scripts/locus_segments.py + samtools bedcov -Q 20 -j)
#   cats.tsv   external catalogue hits (rule catalogue_hits: scripts/match_catalogues.py)
# DEPTH_LOCI=nonpublic (default) measures depth only for the loci that are NOT
# public in the leave-out statuses (LOIO or LOFO): the only loci whose points a
# privacy rule can change. Public loci then have no depth rows, so the scorer
# gives them no copy-number features (their totals are lower bounds). Reading
# depth for every locus touches nearly every CRAM container of a long-read
# alignment (a container spans megabases), i.e. ~80 TB for the panel, which
# the NAS cannot serve in reasonable time. DEPTH_LOCI=all measures every locus.
# Two deviations from the rules, neither changing a depth value:
#   * segment intervals on contigs absent from the alignment header (e.g. a BND
#     mate on 'chrN') are dropped before bedcov, which otherwise exits 2;
#     the count is logged (dropped_intervals)
#   * CRAM/BAM decoding skips SEQ/QUAL/AUX (--input-fmt-option required_fields=0x1FF);
#     bedcov output was byte-identical on HG00096 with and without it
# Usage: sample_depth_catalogues.sh COHORT SAMPLE ALIGNMENT
set -uo pipefail
source "$(dirname "$0")/config.sh"
cohort=$1; sample=$2; aln=$3
S=$SAMPLES_ROOT/$cohort/$sample
O=$CAL_ROOT/persample/$cohort/$sample
mkdir -p "$O"
[ -f "$O/depth_catalogues.ok" ] && exit 0
log=$O/depth_catalogues.log
t0=$(date +%s)
{
  echo "host $(hostname) start $(date -Is)"
  loci=$S/cxsv/cxsv_master_table.tsv
  if [ "${DEPTH_LOCI:-nonpublic}" = nonpublic ]; then
    L=$CAL_ROOT/leaveout/$cohort/$sample
    awk -F'\t' 'FNR==1{for(i=1;i<=NF;i++) if($i=="Panel_Status") k=i; if(NR==1) print "Chr\tStart\tEnd"; next}
                 $k!="public"{print $1"\t"$2"\t"$3}' "$L/status_LOIO.tsv" "$L/status_LOFO.tsv" \
      | awk '!seen[$0]++' > "$O/depth_loci.tsv"
    loci=$O/depth_loci.tsv
    echo "depth_loci nonpublic $(( $(wc -l < "$loci") - 1 ))"
  fi
  $X python "$CODE/workflow/scripts/locus_segments.py" --junctions "$S/junctions/junctions.tsv" \
      --loci "$loci" --fai "$REF.fai" --collapse-distance 300 --composite-fraction 0.25 \
      --flank 20000 --gap 1000 --margin 200 --out "$O/segs.all.bed" || exit 1
  $X samtools view -H "$aln" | awk -F'\t' '$1=="@SQ"{for(i=2;i<=NF;i++) if($i~/^SN:/) print substr($i,4)}' > "$O/sq.txt"
  [ -s "$O/sq.txt" ] || { echo "no @SQ from $aln"; exit 1; }
  awk 'NR==FNR{ok[$1]=1;next} ($1 in ok)' "$O/sq.txt" "$O/segs.all.bed" > "$O/segs.bed"
  echo "dropped_intervals $(( $(wc -l < "$O/segs.all.bed") - $(wc -l < "$O/segs.bed") ))"
  if [ -s "$O/segs.bed" ]; then
    REF_PATH=/nonexistent $X samtools bedcov --input-fmt-option required_fields=0x1FF -Q 20 -j \
        --reference "$REF" "$O/segs.bed" "$aln" > "$O/depth.tsv.tmp" || exit 1
  else
    : > "$O/depth.tsv.tmp"
  fi
  [ "$(wc -l < "$O/depth.tsv.tmp")" -eq "$(wc -l < "$O/segs.bed")" ] || { echo "depth line count mismatch"; exit 1; }
  mv "$O/depth.tsv.tmp" "$O/depth.tsv"
  CATS=""
  for nv in $CATALOGUE_SPECS; do n=${nv%%:*}; v=${nv#*:}; d=$CATALOGUES/$n/$v
    CATS="$CATS --catalogue $n:$v:$d/$n.$v.junction_sites.tsv.gz:$d/$n.$v.manifest.tsv"; done
  $X python "$CODE/workflow/scripts/match_catalogues.py" --junctions "$S/junctions/junctions.tsv" \
      --loci "$S/cxsv/cxsv_master_table.tsv" $CATS --tolerance 500 --ins-len-ratio 0.5 --min-junction-fraction 0.75 \
      --out-loci "$O/cats.tsv" --out-junctions "$O/jcats.tsv" --out-check "$O/cat_check.tsv" || exit 1
  rm -f "$O/segs.all.bed"
  echo "done $(date -Is) wall_s $(( $(date +%s) - t0 ))"
} > "$log" 2>&1 && touch "$O/depth_catalogues.ok" || { echo "FAILED $cohort $sample" >&2; exit 1; }
