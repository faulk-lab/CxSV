#!/usr/bin/env bash
# Relocate the 100 held 1kGP-LRSC BAMs from sargasso2 to sargasso1.
#
# Two passes, deliberately separate:
#   copy   - copy + verify every BAM/BAI pair, leaving the source untouched
#   delete - remove source files, but ONLY those whose destination passed
#
# Deleting 11 TB is not reversible, so the delete pass refuses to act on any
# file that is not recorded as verified in the ledger written by the copy pass.
#
# Verification per file:
#   1. destination size == source size
#   2. samtools quickcheck on the destination BAM (validates header + BGZF EOF,
#      which is what catches a truncated transfer)
#   3. the .bai index copied alongside and non-empty
#
# Resumable: a file already verified in the ledger is skipped, so the script can
# be re-run after an interruption without re-copying.

set -uo pipefail

SRC="/mnt/sargasso2/CxSV/datasets/1KGP_haplotagged_bams"
DST="/mnt/sargasso1/CxSV/datasets/1kGP_LRSC_ONT_CHM13"
LOG_DIR="/mnt/sargasso1/CxSV/datasets/logs"
LEDGER="${LOG_DIR}/relocate_ledger.tsv"
JOBS="${JOBS:-4}"

usage() { echo "usage: $0 {copy|delete|status} [--jobs N]"; exit 2; }
[ $# -ge 1 ] || usage
MODE="$1"; shift
while [ $# -gt 0 ]; do
  case "$1" in
    --jobs) JOBS="$2"; shift 2;;
    *) usage;;
  esac
done

mkdir -p "$DST" "$LOG_DIR"
[ -f "$LEDGER" ] || printf 'bam\tsrc_bytes\tdst_bytes\tquickcheck\tbai_bytes\tstatus\tts\n' > "$LEDGER"

is_verified() {  # $1 = bam basename
  awk -F'\t' -v b="$1" '$1==b && $6=="verified" {found=1} END{exit !found}' "$LEDGER"
}

copy_one() {
  local bam="$1"
  local base; base="$(basename "$bam")"

  if is_verified "$base"; then
    echo "SKIP (already verified): $base"
    return 0
  fi

  local src_sz dst_sz bai_sz qc
  src_sz="$(stat -c %s "$SRC/$base")"

  cp -f "$SRC/$base" "$DST/$base" || { echo "COPY FAILED: $base"; \
    printf '%s\t%s\t\t\t\tcopy_failed\t%s\n' "$base" "$src_sz" "$(date -Is)" >> "$LEDGER"; return 1; }
  cp -f "$SRC/$base.bai" "$DST/$base.bai" || { echo "INDEX COPY FAILED: $base"; \
    printf '%s\t%s\t\t\t\tindex_copy_failed\t%s\n' "$base" "$src_sz" "$(date -Is)" >> "$LEDGER"; return 1; }

  dst_sz="$(stat -c %s "$DST/$base")"
  bai_sz="$(stat -c %s "$DST/$base.bai" 2>/dev/null || echo 0)"
  if samtools quickcheck -v "$DST/$base" 2>/dev/null; then qc="pass"; else qc="FAIL"; fi

  local status="verified"
  [ "$src_sz" = "$dst_sz" ] || status="size_mismatch"
  [ "$qc" = "pass" ]        || status="quickcheck_failed"
  [ "$bai_sz" -gt 0 ]       || status="index_missing"

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$base" "$src_sz" "$dst_sz" "$qc" "$bai_sz" "$status" "$(date -Is)" >> "$LEDGER"
  echo "$status: $base ($(( dst_sz / 1000000000 )) GB)"
  [ "$status" = "verified" ]
}
export -f copy_one is_verified
export SRC DST LEDGER

case "$MODE" in
  copy)
    mapfile -t BAMS < <(find "$SRC" -maxdepth 1 -name '*.bam' | sort)
    echo "copying ${#BAMS[@]} BAMs with $JOBS parallel streams"
    printf '%s\n' "${BAMS[@]}" | parallel -j "$JOBS" --line-buffer copy_one {}
    echo "--- summary ---"
    awk -F'\t' 'NR>1{c[$6]++} END{for(s in c) print s": "c[s]}' "$LEDGER"
    ;;
  delete)
    # Only verified files, and re-confirm the destination is still intact
    # immediately before unlinking the source.
    n=0
    while IFS=$'\t' read -r bam src_sz dst_sz qc bai_sz status ts; do
      [ "$status" = "verified" ] || continue
      [ -f "$DST/$bam" ] || { echo "REFUSE (dest missing): $bam"; continue; }
      cur="$(stat -c %s "$DST/$bam")"
      [ "$cur" = "$src_sz" ] || { echo "REFUSE (dest size drifted): $bam"; continue; }
      samtools quickcheck "$DST/$bam" || { echo "REFUSE (dest quickcheck): $bam"; continue; }
      rm -f "$SRC/$bam" "$SRC/$bam.bai" && { echo "removed source: $bam"; n=$((n+1)); }
    done < <(tail -n +2 "$LEDGER")
    echo "removed $n source BAM(s)"
    ;;
  status)
    awk -F'\t' 'NR>1{c[$6]++} END{for(s in c) print s": "c[s]}' "$LEDGER"
    echo "dest BAMs: $(find "$DST" -maxdepth 1 -name '*.bam' | wc -l)"
    echo "src  BAMs: $(find "$SRC" -maxdepth 1 -name '*.bam' 2>/dev/null | wc -l)"
    ;;
  *) usage;;
esac
