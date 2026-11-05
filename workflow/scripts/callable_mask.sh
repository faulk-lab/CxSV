#!/usr/bin/env bash
# =============================================================================
# callable_mask.sh — CALLABLE intervals of one sample (rule callable_mask,
# workflow/rules/panel_callability.smk; definitions in lib/cxsv_callability.py).
#
#   callable_mask.sh ALN REF PREFIX MASK MIN_DEPTH MAPQ THREADS [CRAM31]
#
# mosdepth quantizes the depth of reads with MAPQ >= MAPQ into LOW / CALLABLE
# at MIN_DEPTH; the CALLABLE intervals are written to MASK (gzip -n, BED3).
# Writes PREFIX.mosdepth.summary.txt; the quantized BED is removed.
#
# CRAM 3.1 input: the pinned mosdepth 0.3.10 is statically linked with an
# htslib that cannot decode CRAM 3.1. It logs "Slice decode failure" for every
# contig, exits 0 and leaves an EMPTY mask (208 of 1,017 1kGP ONT CRAMs in the
# panel-v1.0 build; the build transcoded them at staging). CRAM31 decides:
#   transcode (default)  rewrite as CRAM 3.0 into PREFIX.v30.cram (samtools,
#                        lossless; record counts must match), use it, delete
#                        it afterwards; the input is never modified
#   fail                 stop with an error naming the fix
# Whatever the input: an htslib/mosdepth decode error in mosdepth's output, or
# a mask with 0 callable bp, fails the step (exit 1) instead of shipping an
# empty mask.
# For CRAMs, REF_PATH / REF_CACHE must point at the MD5 cache (HTS_ENV).
# =============================================================================
set -euo pipefail
[ $# -ge 7 ] || { sed -n '4,6p' "$0" >&2; exit 2; }
ALN="$1" REF="$2" PREFIX="$3" MASK="$4" MIN_DEPTH="$5" MAPQ="$6" THREADS="$7" CRAM31="${8:-transcode}"
say() { echo "callable_mask: $*"; }

cram_version() {        # "3.1", "3.0", ... for a CRAM; "" otherwise
    local b; b=$(head -c 6 "$1" | od -An -tu1)
    set -- $b
    if [ "${1:-}" = 67 ] && [ "${2:-}" = 82 ] && [ "${3:-}" = 65 ] && [ "${4:-}" = 77 ]; then echo "$5.$6"; fi
}

IN="$ALN"; TMP=""
cleanup() { [ -n "$TMP" ] && rm -f "$TMP" "$TMP.crai"; rm -f "$MASK.tmp"; }
trap cleanup EXIT
V="$(cram_version "$ALN")"
if [ -n "$V" ]; then say "input $ALN (CRAM $V)"; else say "input $ALN (not CRAM)"; fi
mkdir -p "$(dirname "$PREFIX")"
if [ "$V" = 3.1 ]; then
    case "$CRAM31" in
        fail)
            say "ERROR: $ALN is CRAM 3.1; mosdepth 0.3.10 cannot decode CRAM 3.1 (it writes an empty mask" \
                "and exits 0). Rewrite the file as CRAM 3.0 (samtools view -C --output-fmt-option version=3.0)" \
                "or set panel_callability.cram31: transcode." >&2
            exit 1 ;;
        transcode)
            TMP="$PREFIX.v30.cram"; rm -f "$TMP" "$TMP.crai"
            samtools view -@ "$THREADS" -C --output-fmt-option version=3.0 -T "$REF" \
                --write-index -o "$TMP##idx##$TMP.crai" "$ALN"
            n_in=$(samtools view -@ "$THREADS" -c -T "$REF" "$ALN")
            n_out=$(samtools view -@ "$THREADS" -c -T "$REF" "$TMP")
            if [ -z "$n_in" ] || [ "$n_in" != "$n_out" ]; then
                say "ERROR: CRAM 3.1 -> 3.0 record counts differ ($n_in vs $n_out)" >&2; exit 1
            fi
            say "transcoded CRAM 3.1 -> 3.0 for mosdepth: $n_in records ($TMP, removed afterwards)"
            IN="$TMP" ;;
        *) say "ERROR: panel_callability.cram31 must be transcode or fail, not '$CRAM31'" >&2; exit 2 ;;
    esac
fi

set +e
MOSDEPTH_Q0=LOW MOSDEPTH_Q1=CALLABLE mosdepth -t "$THREADS" -n -x --mapq "$MAPQ" \
    --quantize "0:$MIN_DEPTH:" -f "$REF" "$PREFIX" "$IN" > "$PREFIX.mosdepth.out" 2>&1
rc=$?
set -e
cat "$PREFIX.mosdepth.out"
if [ $rc -ne 0 ]; then say "ERROR: mosdepth exited $rc" >&2; rm -f "$PREFIX.mosdepth.out"; exit 1; fi
if grep -qiE "decode failure|failed to decode|\[E::|truncated file|failed to populate reference" "$PREFIX.mosdepth.out"; then
    say "ERROR: mosdepth reported decode errors (above); the mask would be incomplete" >&2
    rm -f "$PREFIX.mosdepth.out"; exit 1
fi
rm -f "$PREFIX.mosdepth.out"

zcat "$PREFIX.quantized.bed.gz" \
    | awk -F'\t' -v OFS='\t' '$4 == "CALLABLE" {print $1, $2, $3}' \
    | gzip -n > "$MASK.tmp"
rm -f "$PREFIX.quantized.bed.gz" "$PREFIX.quantized.bed.gz.csi"
# %.0f, not %d: mawk's %d stops at 2^31-1 (2147483647), and a genome mask is ~3 Gb
bp=$(zcat "$MASK.tmp" | awk -F'\t' '{s += $3 - $2} END {printf "%.0f\n", s}')
if [ "$bp" -le 0 ]; then
    say "ERROR: the mask has 0 callable bp (depth >= $MIN_DEPTH at MAPQ >= $MAPQ nowhere);" \
        "an empty mask would make every junction uncallable in this sample" >&2
    exit 1
fi
mv "$MASK.tmp" "$MASK"
say "callable_bp $bp"
