#!/usr/bin/env bash
# =============================================================================
# depth_pass.sh — binned read depth of one genome (roadmap item 15, Phase 1;
# docs/CN_LINKAGE_PHASE1.md). Prototype outside the pipeline.
#
#   depth_pass.sh ALN OUTDIR [THREADS] [BIN]
#
# mosdepth 0.3.10 from the pipeline image (SIF, env SIF), reads with MAPQ >= 20,
# no per-base output (-n), fast mode (-x, as callable_mask.sh), --by BIN (default
# 1000; the 10 kb bins of the plan are sums of ten 1 kb bins, so one pass gives both).
# CRAM 3.1 is rewritten as CRAM 3.0 into OUTDIR/tmp first (the transcode path of
# workflow/scripts/callable_mask.sh, record counts checked), then deleted.
# Outputs in OUTDIR: depth.regions.bed.gz (+ .csi), depth.mosdepth.summary.txt,
# timing.tsv (host, bytes, seconds per step, MB/s = input bytes / total seconds), depth.ok.
# Checks: no decode error in mosdepth's log; one bin per BIN bp of every contig
# in the .fai (count compared), or the step fails.
# =============================================================================
set -euo pipefail
ALN=$1; OUT=$2; T=${3:-4}; BIN=${4:-1000}
SIF=${SIF:-/mnt/sargasso1/CxSV/bundles/cxsv-bundle-20261105-8fc5bef725bd/container/cxsv-pipeline.sif}
REF=${REF:-/mnt/sargasso1/CxSV/refs/hs1.fa}
CACHE=${CACHE:-/mnt/sargasso1/CxSV/refs/ref_cache}
RT=$(command -v singularity || command -v apptainer)
mkdir -p "$OUT/tmp"; rm -f "$OUT/depth.ok"
X() { "$RT" exec -B /mnt/sargasso1 -B "$(dirname "$ALN")" -B "$OUT" \
      --env "REF_PATH=$CACHE/%2s/%2s/%s,REF_CACHE=$CACHE/%2s/%2s/%s" "$SIF" "$@"; }
cram_version() { local b; b=$(head -c 6 "$1" | od -An -tu1); set -- $b
                 [ "${1:-}" = 67 ] && [ "${2:-}" = 82 ] && [ "${3:-}" = 65 ] && [ "${4:-}" = 77 ] && echo "$5.$6" || true; }
bytes=$(stat -L -c %s "$ALN"); V=$(cram_version "$ALN"); IN=$ALN; t_tr=0; TMP=""
t0=$(date +%s)
if [ "$V" = 3.1 ]; then
    TMP="$OUT/tmp/v30.cram"
    X samtools view -@ "$T" -C --output-fmt-option version=3.0 -T "$REF" --write-index -o "$TMP##idx##$TMP.crai" "$ALN"
    a=$(X samtools view -@ "$T" -c -T "$REF" "$ALN"); b=$(X samtools view -@ "$T" -c -T "$REF" "$TMP")
    [ -n "$a" ] && [ "$a" = "$b" ] || { echo "transcode record counts differ: $a vs $b" >&2; exit 1; }
    IN=$TMP; t_tr=$(( $(date +%s) - t0 ))
fi
t1=$(date +%s)
X mosdepth -t "$T" -n -x --mapq 20 --by "$BIN" -f "$REF" "$OUT/depth" "$IN" > "$OUT/mosdepth.log" 2>&1
t_md=$(( $(date +%s) - t1 ))
if grep -qiE "decode failure|failed to decode|\[E::|truncated file|failed to populate reference" "$OUT/mosdepth.log"; then
    echo "mosdepth decode errors" >&2; cat "$OUT/mosdepth.log" >&2; exit 1; fi
[ -n "$TMP" ] && rm -f "$TMP" "$TMP.crai"
want=$(awk -v b="$BIN" '{n+=int(($2+b-1)/b)} END{print n}' "$REF.fai")
got=$(zcat "$OUT/depth.regions.bed.gz" | wc -l)
[ "$want" = "$got" ] || { echo "bin count $got != expected $want" >&2; exit 1; }
rm -f "$OUT/depth.mosdepth.global.dist.txt" "$OUT/depth.mosdepth.region.dist.txt"
tot=$(( $(date +%s) - t0 ))
printf 'host\taln\tbytes\tcram_version\ttranscode_s\tmosdepth_s\ttotal_s\tMB_per_s\tbins\tthreads\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$(hostname -s)" "$ALN" "$bytes" "${V:-none}" "$t_tr" "$t_md" "$tot" \
    "$(awk -v b="$bytes" -v t="$tot" 'BEGIN{printf "%.1f", b/1e6/(t>0?t:1)}')" "$got" "$T" > "$OUT/timing.tsv"
touch "$OUT/depth.ok"
