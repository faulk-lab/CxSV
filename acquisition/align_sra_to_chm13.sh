#!/usr/bin/env bash
# Align one SRA run to T2T-CHM13v2 (hs1), reproducing the recipe used for the
# 100 held 1kGP-LRSC BAMs so that calls are comparable across the background set.
#
#   bash align_sra_to_chm13.sh ERR15139516.sra NA07019
#
# The held BAMs were produced by the CARD cardEndToEndVcfMethyl pipeline with:
#   minimap2 -t 10 -y -x map-ont -a --eqx -k 17 -K 10g T2T-CHM13v2.fasta
# We match -x map-ont, -a, --eqx, -k 17 (baked into the prebuilt index) and -K 10g.
#
# Deliberate difference: **-y is omitted.** It copies FASTQ comments into SAM
# tags, which is how the original carried MM/ML methylation tags over from an
# unaligned modBAM. Our input is SRA-derived FASTQ with no modification tags,
# where -y would only copy fasterq-dump's "length=NNN" comment into every record
# as noise. Methylation is not used for CxSV calling, so nothing is lost.
#
# --eqx is not cosmetic: it emits =/X instead of M in the CIGAR, and the SV
# callers in this pipeline were run against BAMs that have it. Dropping it would
# make these samples differ from the held 100 in a way that affects calling.

set -euo pipefail

SRA_FILE="${1:?usage: align_sra_to_chm13.sh <file.sra> <sample_id>}"
SAMPLE="${2:?usage: align_sra_to_chm13.sh <file.sra> <sample_id>}"

REF_MMI="${REF_MMI:-${REF_DIR:-/home/minknow/refs/hs1}/hs1.map-ont.k17.mmi}"
REF_FA="${REF_FA:-${REF_DIR:-/home/minknow/refs/hs1}/hs1.fa}"
OUT_DIR="${OUT_DIR:-/mnt/sargasso1/CxSV/datasets/1kGP_ONT_1019_CHM13}"
LOG_DIR="${LOG_DIR:-/mnt/sargasso1/CxSV/datasets/logs/align}"
# Scratch must be local: sorting through a CIFS mount is far slower and the
# temp files are written and deleted in the same minute.
SCRATCH="${SCRATCH:-/home/minknow/scratch/cxsv_align}"
THREADS_MM="${THREADS_MM:-24}"
THREADS_SORT="${THREADS_SORT:-8}"
SORT_MEM="${SORT_MEM:-2G}"
THREADS_FQ="${THREADS_FQ:-6}"
# Which conda env supplies fasterq-dump/minimap2/samtools. freya, shadowfax and
# athena carry minimap2 2.24 system-wide and only need sra-tools; betsy shipped
# 2.27, so it gets a pinned env instead. All four must align with the SAME
# minimap2 version as the held 100 BAMs (2.24-r1122) or the background set
# acquires a version batch effect that looks like private variation.
CONDA_ENV="${CONDA_ENV:-sratools}"
# minimap2 minibatch size. The held BAMs used -K 10g; shadowfax has 61 GB RAM
# rather than 125, so it runs with a smaller batch. This changes throughput, not
# the alignments themselves.
MM_BATCH="${MM_BATCH:-4g}"

ACC="$(basename "$SRA_FILE" .sra)"
OUT_CRAM="${OUT_DIR}/${SAMPLE}.${ACC}.hs1.minimap2.cram"
TMP="${SCRATCH}/${ACC}"
LOG="${LOG_DIR}/${SAMPLE}.${ACC}.log"

mkdir -p "$OUT_DIR" "$LOG_DIR" "$TMP"
trap 'rm -rf "$TMP"' EXIT

if [ -s "${OUT_CRAM}.crai" ] && samtools quickcheck "$OUT_CRAM" 2>/dev/null; then
  echo "SKIP (already aligned and indexed): $OUT_CRAM"
  exit 0
fi

source ~/miniforge3/etc/profile.d/conda.sh
conda activate "$CONDA_ENV"

MM_VERSION="$(minimap2 --version)"
if [ "$MM_VERSION" != "2.24-r1122" ]; then
  echo "REFUSING: minimap2 is $MM_VERSION, expected 2.24-r1122 (held BAMs' version)" >&2
  exit 4
fi

echo "=== $SAMPLE / $ACC  started $(date -Is)" | tee "$LOG"
START=$(date +%s)

# fasterq-dump needs --concatenate-reads to stream to stdout at all: its default
# split-3 mode would emit multiple files and it refuses stdout in that mode.
set -o pipefail
fasterq-dump --concatenate-reads --stdout --threads "$THREADS_FQ" --temp "$TMP" "$SRA_FILE" \
  2>> "$LOG" \
  | minimap2 -t "$THREADS_MM" -x map-ont -a --eqx -K "$MM_BATCH" \
      -R "@RG\tID:${ACC}\tSM:${SAMPLE}\tPL:ONT\tPU:${ACC}\tDS:PRJEB89727_SRA_mirror" \
      "$REF_MMI" - 2>> "$LOG" \
  | samtools sort -@ "$THREADS_SORT" -m "$SORT_MEM" -T "${TMP}/sort" \
      --reference "$REF_FA" -O cram -o "${OUT_CRAM}.tmp" 2>> "$LOG"

# Only name it as the final CRAM once it passes a structural check, so an
# interrupted run never leaves a truncated file that looks finished.
samtools quickcheck -v "${OUT_CRAM}.tmp"
mv "${OUT_CRAM}.tmp" "$OUT_CRAM"
samtools index -@ "$THREADS_SORT" "$OUT_CRAM"

END=$(date +%s)
BYTES=$(stat -c %s "$OUT_CRAM")
echo "=== done in $((END-START))s, $((BYTES/1000000000)) GB" | tee -a "$LOG"

# Per-sample QC: mapped fraction and achieved depth, written next to the BAM so
# the registry can be populated without re-reading the alignment.
samtools flagstat -@ 4 "$OUT_CRAM" > "${OUT_DIR}/${SAMPLE}.${ACC}.flagstat.txt"
if command -v mosdepth >/dev/null; then
  ( cd "$TMP" && mosdepth -n --fast-mode -t 4 -f "$REF_FA" "${SAMPLE}" "$OUT_CRAM" \
      && cp "${SAMPLE}.mosdepth.summary.txt" "${OUT_DIR}/${SAMPLE}.${ACC}.mosdepth.summary.txt" )
fi

# Integrity gate for the later .sra cleanup pass: the number of primary
# alignments must equal the read count the archive itself reports. Recorded here
# so the delete pass can check it without re-reading either file.
SRA_READS=$(sra-stat --quick "$SRA_FILE" 2>/dev/null | awk -F'|' '{print $3}' | awk -F: '{print $1}')
CRAM_READS=$(samtools view -c -F 0x900 --reference "$REF_FA" "$OUT_CRAM")
# One file per sample, written to a temp name and renamed into place.
#
# Do NOT append to a shared ledger here. Every node runs this script, and
# O_APPEND is not atomic across SMB clients: a client holding a stale cached
# file size writes at the old end-of-file and overwrites another node's row.
# That corrupted the shared align_integrity.tsv (NUL-filled rows, rows lost)
# once six nodes were appending to it. rename() is atomic on SMB, so a
# per-sample file can never be clobbered by a concurrent writer.
INTEG_DIR="${LOG_DIR}/integrity"
mkdir -p "$INTEG_DIR"
VERDICT=$([ -n "$SRA_READS" ] && [ "$SRA_READS" = "$CRAM_READS" ] && echo read_count_match || echo READ_COUNT_MISMATCH)
TMPREC="${INTEG_DIR}/.${ACC}.$$.tmp"
printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$SAMPLE" "$ACC" "$SRA_READS" "$CRAM_READS" "$VERDICT" "$BYTES" > "$TMPREC"
mv -f "$TMPREC" "${INTEG_DIR}/${ACC}.tsv"
echo "OK $SAMPLE $ACC $((END-START))s ${BYTES} sra_reads=${SRA_READS} cram_reads=${CRAM_READS}"
