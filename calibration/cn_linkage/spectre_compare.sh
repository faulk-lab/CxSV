#!/usr/bin/env bash
# External long-read CNV caller as comparator (roadmap item 15 Phase 1): Spectre
# (Sedlazeck lab, v0.2.1) on the depth pass's 1 kb mosdepth output.
#
#   spectre_compare.sh DEPTH_DIR SAMPLE SEX OUTDIR [REF] [META] [BLACKLIST]
#
# DEPTH_DIR holds depth.regions.bed.gz(.csi) (depth_pass.sh). SEX male|female: male
# runs with --ploidy-chr chrX:1,chrY:1; female on chr1-22,X only (Spectre's default
# treats every contig as diploid).
# META: Spectre RemoveNs metadata of the reference (made once, cached next to OUTDIR).
# hs1 has no N runs, so RemoveNs masks nothing; BLACKLIST (BED, e.g. the PoN-masked
# bins from spectre_blacklist.py) is passed as --blacklist (Phase 2 sanity check).
# Output: OUTDIR/<SAMPLE>/ (Spectre VCF/BED) and OUTDIR/<SAMPLE>.cnv.bed
# (chrom start end type copies), the form the comparison in CN_LINKAGE_PHASE1.md uses.
set -euo pipefail
DEPTH=$1; S=$2; SEX=$3; OUT=$4
REF=${5:-/mnt/sargasso1/CxSV/refs/hs1.fa}
META=${6:-$(dirname "$OUT")/spectre_meta/hs1.mdr}
BL=${7:-}
mkdir -p "$OUT/$S" "$(dirname "$META")"
if [ ! -s "$META" ]; then
  spectre RemoveNs --reference "$REF" --output-dir "$(dirname "$META")" --output-file "$(basename "$META")" --save-only
fi
# Spectre reads <dir>/<prefix>.regions.bed.gz; give it a private copy dir
W="$OUT/$S/cov"; mkdir -p "$W"
cp "$DEPTH/depth.regions.bed.gz" "$W/$S.regions.bed.gz"
[ -f "$DEPTH/depth.regions.bed.gz.csi" ] && cp "$DEPTH/depth.regions.bed.gz.csi" "$W/$S.regions.bed.gz.csi"
CHR=$(printf 'chr%s,' $(seq 1 22))chrX
if [ "$SEX" = male ]; then PL=(--ploidy-chr chrX:1,chrY:1 --only-chr "${CHR},chrY"); else PL=(--only-chr "$CHR"); fi
t0=$(date +%s)
spectre CNVCaller --coverage "$W/$S.regions.bed.gz" --sample-id "$S" --output-dir "$OUT/$S" --reference "$REF" \
  --metadata "$META" "${PL[@]}" ${BL:+--blacklist "$BL"} > "$OUT/$S/spectre.log" 2>&1
echo -e "$S\t$(( $(date +%s)-t0 ))" > "$OUT/$S/wall_s.tsv"
# VCF -> simple BED
V=$(ls "$OUT/$S"/*.vcf.gz "$OUT/$S"/*.vcf 2>/dev/null | head -1 || true)
if [ -n "$V" ]; then
  (zcat -f "$V" | awk -F'\t' -v OFS='\t' '!/^#/{
      end=""; t=""; cn="."
      n=split($8,a,";"); for(i=1;i<=n;i++){split(a[i],b,"="); if(b[1]=="END") end=b[2]; if(b[1]=="SVTYPE") t=b[2]}
      for(i=1;i<=n;i++){split(a[i],b,"="); if(b[1]=="CN") cn=b[2]}      # Spectre: CN is in INFO
      if(cn=="."){ split($9,f,":"); split($10,v,":"); for(i=1;i in f;i++) if(f[i]=="CN") cn=v[i] }
      print $1,$2,end,t,cn}') > "$OUT/$S.cnv.bed"
fi
