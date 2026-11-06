#!/usr/bin/env bash
# =============================================================================
# cn_depth.sh — 1 kb binned read depth of one sample for the copy-number track
# (rule cn_depth, workflow/rules/cn.smk; docs/CN_LINKAGE.md).
#
#   cn_depth.sh ALN REF OUTDIR BIN MAPQ THREADS [PRECOMPUTED_DIR]
#
# mosdepth (-n, fast mode -x as callable_mask.sh, --mapq MAPQ, --by BIN) writes
# OUTDIR/depth.regions.bed.gz (+ .csi) and OUTDIR/depth.mosdepth.summary.txt.
# CRAM 3.1 is rewritten as CRAM 3.0 first (the pinned mosdepth cannot decode 3.1;
# record counts checked; the copy is deleted). Checks: no decode error in mosdepth's
# log; one bin per BIN bp of every contig of REF.fai, or the step fails.
# OUTDIR/depth.sha256 records the sha256 of depth.regions.bed.gz.
#
# PRECOMPUTED_DIR (optional; the full-panel pass of cn_depth_v1, docs/CN_LINKAGE.md):
# a directory with the same files from the same pass (depth.ok, depth.sha256). Its
# files are hard-linked (copied if the filesystem refuses) after the sha256 and the
# bin count are checked, instead of reading the alignment again.
# For CRAMs, REF_PATH / REF_CACHE must point at the MD5 cache (HTS_ENV).
# =============================================================================
set -euo pipefail
[ $# -ge 6 ] || { sed -n '4,6p' "$0" >&2; exit 2; }
ALN="$1" REF="$2" OUT="$3" BIN="$4" MAPQ="$5" THREADS="$6" PRE="${7:-}"
say() { echo "cn_depth: $*"; }
mkdir -p "$OUT"
want=$(awk -v b="$BIN" '{n+=int(($2+b-1)/b)} END{print n}' "$REF.fai")
check() {   # bin count of $1
    local got; got=$(zcat "$1" | wc -l)
    [ "$got" = "$want" ] || { say "ERROR: bin count $got != expected $want ($1)" >&2; exit 1; }
}
if [ -n "$PRE" ] && [ "$PRE" != "none" ]; then
    [ -f "$PRE/depth.ok" ] && [ -s "$PRE/depth.sha256" ] || { say "ERROR: $PRE has no finished depth (depth.ok, depth.sha256)" >&2; exit 1; }
    sha=$(sha256sum "$PRE/depth.regions.bed.gz" | awk '{print $1}')
    [ "$sha" = "$(cat "$PRE/depth.sha256")" ] || { say "ERROR: sha256 of $PRE/depth.regions.bed.gz differs from depth.sha256" >&2; exit 1; }
    check "$PRE/depth.regions.bed.gz"
    for f in depth.regions.bed.gz depth.regions.bed.gz.csi depth.mosdepth.summary.txt; do
        rm -f "$OUT/$f"; ln "$PRE/$f" "$OUT/$f" 2>/dev/null || cp "$PRE/$f" "$OUT/$f"
    done
    echo "$sha" > "$OUT/depth.sha256"
    say "precomputed depth from $PRE (sha256 $sha)"
    exit 0
fi
cram_version() { local b; b=$(head -c 6 "$1" | od -An -tu1); set -- $b
    if [ "${1:-}" = 67 ] && [ "${2:-}" = 82 ] && [ "${3:-}" = 65 ] && [ "${4:-}" = 77 ]; then echo "$5.$6"; fi; }
IN="$ALN"; TMP=""
cleanup() { [ -n "$TMP" ] && rm -f "$TMP" "$TMP.crai"; return 0; }
trap cleanup EXIT
V="$(cram_version "$ALN")"
if [ "$V" = 3.1 ]; then
    TMP="$OUT/v30.cram"; rm -f "$TMP" "$TMP.crai"
    samtools view -@ "$THREADS" -C --output-fmt-option version=3.0 -T "$REF" --write-index -o "$TMP##idx##$TMP.crai" "$ALN"
    a=$(samtools view -@ "$THREADS" -c -T "$REF" "$ALN"); b=$(samtools view -@ "$THREADS" -c -T "$REF" "$TMP")
    [ -n "$a" ] && [ "$a" = "$b" ] || { say "ERROR: CRAM 3.1 -> 3.0 record counts differ ($a vs $b)" >&2; exit 1; }
    say "transcoded CRAM 3.1 -> 3.0 for mosdepth: $a records"
    IN="$TMP"
fi
mosdepth -t "$THREADS" -n -x --mapq "$MAPQ" --by "$BIN" -f "$REF" "$OUT/depth" "$IN" > "$OUT/mosdepth.log" 2>&1
if grep -qiE "decode failure|failed to decode|\[E::|truncated file|failed to populate reference" "$OUT/mosdepth.log"; then
    say "ERROR: mosdepth decode errors" >&2; cat "$OUT/mosdepth.log" >&2; exit 1
fi
rm -f "$OUT/depth.mosdepth.global.dist.txt" "$OUT/depth.mosdepth.region.dist.txt"
check "$OUT/depth.regions.bed.gz"
sha256sum "$OUT/depth.regions.bed.gz" | awk '{print $1}' > "$OUT/depth.sha256"
say "$(cat "$OUT/depth.sha256") $want bins"
