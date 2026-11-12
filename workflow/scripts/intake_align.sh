#!/usr/bin/env bash
# =============================================================================
# intake_align.sh — align one patient sample to hs1 with the panel recipe
# (roadmap item 16k; docs/INTAKE.md). Called by workflow/intake/Snakefile; also
# runs on its own:
#
#   REF_FA=hs1.fa REF_MMI=hs1.map-ont.k17.mmi \
#     bash workflow/scripts/intake_align.sh <sample_id> <out_dir> <input> [<input> ...]
#
# Inputs: unaligned BAM / dorado modBAM, FASTQ (.fastq/.fq, optionally .gz), or a
# BAM/CRAM from another aligner or reference naming. Several inputs (flow cells)
# are concatenated into one alignment. Only primary records of BAM/CRAM inputs
# are used (samtools view -F 0x900), so each read enters exactly once.
#
# Recipe (the panel's, docs/POSITIVE_CONTROL.md step 1; on BOC_4433 it reproduced the
# panel-recipe CRAM record for record, docs/INTAKE.md):
#   samtools fastq -T MM,ML,MN,qs | minimap2 2.24-r1122 -y -ax map-ont --eqx -k 17 -K <batch>
#     | samtools sort | CRAM (reference hs1.fa), indexed
# -y copies the MM/ML/MN modified-base tags and dorado's qs into the alignments.
# Only those four tags are copied: '-T *' would also copy dorado's per-read RG:Z,
# which would clash with the read group written here. FASTQ comments are reduced
# to the same four tags (intake_tool.py fastq-tags), because -y copies the comment
# verbatim and non-SAM tokens (runid=..., ch=...) would make invalid SAM.
# -K (minibatch) changes memory and speed only, not the alignments.
# The index must be a k=17 w=10 map-ont index of the same hs1 (checked from its header).
# The panel CRAMs' @PG has no '-k 17' (k came from the prebuilt index); here -k 17 is also
# written on the command line, which minimap2 overrides with the index's k, so the
# alignments are the same and the @PG records k explicitly.
#
# Read group: ID = SM = <sample_id> (the sample sheet's ID, so caller VCF sample
# columns match the patient sheet), PL, DS = intake description + recipe.
#
# Checks (any failure exits non-zero and removes the CRAM):
#   minimap2 is 2.24-r1122; index header k=17 w=10 and one sequence per .fai line;
#   input primary reads == output primary mapped + unmapped; output header has one
#   @RG with SM == sample_id and a minimap2 @PG with the recipe; if the first 2,000
#   input records carry MM/ML, the first 2,000 output records do too.
# Outputs in <out_dir>: <sample>.hs1.cram(.crai), <sample>.flagstat.txt,
#   <sample>.idxstats.tsv, input_header.txt, provenance.tsv, intake.json, logs/.
#
# Environment (defaults in brackets):
#   REF_FA, REF_MMI   hs1 FASTA (with .fai) and minimap2 index        [required]
#   INPUT_REF         FASTA to decode CRAM inputs made against another reference
#                     (otherwise htslib uses REF_PATH, i.e. the hs1 MD5 cache)
#   THREADS [16]  MM_BATCH [4g]  SORT_THREADS [8]  SORT_MEM [2G]
#   PLATFORM [ONT]  DESCRIPTION [""]  TMP_DIR [<out_dir>/.tmp]
# =============================================================================
set -euo pipefail
SAMPLE="${1:?usage: intake_align.sh <sample_id> <out_dir> <input> [<input> ...]}"
OUT="${2:?usage: intake_align.sh <sample_id> <out_dir> <input> [<input> ...]}"
shift 2
[ $# -ge 1 ] || { echo "no input given" >&2; exit 1; }
INPUTS=("$@")
: "${REF_FA:?REF_FA (hs1 FASTA) is required}"; : "${REF_MMI:?REF_MMI (minimap2 index) is required}"
THREADS="${THREADS:-16}"; MM_BATCH="${MM_BATCH:-4g}"; SORT_THREADS="${SORT_THREADS:-8}"; SORT_MEM="${SORT_MEM:-2G}"
PLATFORM="${PLATFORM:-ONT}"; DESCRIPTION="${DESCRIPTION:-}"; INPUT_REF="${INPUT_REF:-}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOL="${TOOL:-python3 $HERE/intake_tool.py}"
RECIPE_MM="2.24-r1122"

mkdir -p "$OUT/logs"
TMP_DIR="${TMP_DIR:-$OUT/.tmp}"; T="$TMP_DIR/intake_${SAMPLE}.$$"; mkdir -p "$T"
CRAM="$OUT/$SAMPLE.hs1.cram"
cleanup() {
  rc=$?; rm -rf "$T"
  if [ $rc -ne 0 ]; then
    rm -f "$CRAM" "$CRAM.crai"
    for l in "$OUT"/logs/*.stderr; do [ -s "$l" ] && { echo "== $l"; tail -n 5 "$l"; } >&2; done
    echo "intake $SAMPLE FAILED (exit $rc)" >&2
  fi
}
trap cleanup EXIT

MM_VERSION="$(minimap2 --version)"
[ "$MM_VERSION" = "$RECIPE_MM" ] || { echo "REFUSING: minimap2 is $MM_VERSION, the panel recipe needs $RECIPE_MM" >&2; exit 4; }
$TOOL mmi-check "$REF_MMI" --fai "$REF_FA.fai" > "$OUT/logs/mmi_check.json"

fmt_of() { case "${1,,}" in *.fastq|*.fq|*.fastq.gz|*.fq.gz) echo fastq ;; *.bam) echo bam ;; *.cram) echo cram ;;
           *) echo "unsupported input $1" >&2; return 1 ;; esac; }
VIEW_REF=(); [ -n "$INPUT_REF" ] && VIEW_REF=(--reference "$INPUT_REF")

{
  echo -e "sample\t$SAMPLE"; echo -e "inputs\t$(IFS=,; echo "${INPUTS[*]}")"
  for f in "${INPUTS[@]}"; do echo -e "input_bytes\t$(basename "$f")\t$(stat -c %s "$f")"; done
  echo -e "node\t$(hostname)"; echo -e "start\t$(date -Is)"
  echo -e "minimap2\t$MM_VERSION"; echo -e "samtools\t$(samtools --version | head -1)"
  echo -e "ref_fa\t$REF_FA"; echo -e "ref_mmi\t$REF_MMI"; echo -e "input_reference\t${INPUT_REF:-none}"
  echo -e "threads\t$THREADS"; echo -e "mm_batch\t$MM_BATCH"; echo -e "sort\t-@ $SORT_THREADS -m $SORT_MEM"
} > "$OUT/provenance.tsv"

# headers of the inputs (alignment inputs), kept for the record
: > "$OUT/input_header.txt"
for f in "${INPUTS[@]}"; do
  if [ "$(fmt_of "$f")" != fastq ]; then
    { echo "# $f"; samtools view -H "${VIEW_REF[@]}" "$f" | grep -v '^@SQ'; } >> "$OUT/input_header.txt"
  fi
done

# preflight: every BAM/CRAM input must decode (a CRAM made against another
# reference needs input_reference); the first 2,000 primary records of the first
# input also give the modified-base tag sample
i=0
for f in "${INPUTS[@]}"; do
  i=$((i + 1))
  [ "$(fmt_of "$f")" = fastq ] && continue
  { samtools view -F 0x900 "${VIEW_REF[@]}" "$f" 2> "$T/peek.$i.err" || true; } | head -n 2000 > "$T/peek.$i.sam" || true
  if grep -q -E '^\[E::|[Ff]ail' "$T/peek.$i.err"; then
    echo "cannot decode $f: $(grep -m1 -E '^\[E::|[Ff]ail' "$T/peek.$i.err")" >&2
    echo "  (a CRAM aligned to another reference needs input_reference = that FASTA in the intake sheet)" >&2
    exit 3
  fi
done
first="${INPUTS[0]}"
case "$(fmt_of "$first")" in
  fastq) MOD_IN=$( { if [[ "${first,,}" == *.gz ]]; then gzip -dc "$first"; else cat "$first"; fi; } 2>/dev/null \
                   | head -n 8000 | awk 'NR % 4 == 1' | grep -c 'MM:Z:' || true) ;;
  *)     MOD_IN=$(grep -c $'\tMM:Z:' "$T/peek.1.sam" || true) ;;
esac

RG_DS="cxsv-intake: ${DESCRIPTION:+$DESCRIPTION; }minimap2 $RECIPE_MM -y -ax map-ont --eqx -k 17 to hs1"
RG="@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:${PLATFORM}\tDS:${RG_DS}"

# reads written by samtools fastq ("[M::bam2fq_mainloop] processed N reads")
count_reads() { python3 - "$1" <<'PY'
import re, sys
m = re.findall(r"processed (\d+) reads", open(sys.argv[1]).read())
if not m:
    sys.exit("samtools fastq reported no read count")
print(m[-1])
PY
}

# one stream of reads from every input; input md5s from the same read pass
stream_reads() {
  local i=0
  for f in "${INPUTS[@]}"; do
    i=$((i + 1)); mkfifo "$T/in.$i.fifo"
    md5sum < "$T/in.$i.fifo" | awk '{print $1}' > "$T/in.$i.md5" &
    case "$(fmt_of "$f")" in
      fastq)
        if [[ "${f,,}" == *.gz ]]; then
          tee "$T/in.$i.fifo" < "$f" | gzip -dc | $TOOL fastq-tags --count-out "$T/count.$i" - 2>> "$OUT/logs/fastq.stderr"
        else
          tee "$T/in.$i.fifo" < "$f" | $TOOL fastq-tags --count-out "$T/count.$i" - 2>> "$OUT/logs/fastq.stderr"
        fi ;;
      bam|cram)
        tee "$T/in.$i.fifo" < "$f" | samtools view -u -F 0x900 "${VIEW_REF[@]}" - 2>> "$OUT/logs/view.stderr" \
          | samtools fastq -T MM,ML,MN,qs - 2> "$T/fastq.$i.stderr"
        cat "$T/fastq.$i.stderr" >> "$OUT/logs/fastq.stderr"
        count_reads "$T/fastq.$i.stderr" > "$T/count.$i" ;;
    esac
  done
  wait
}

START=$(date +%s)
stream_reads \
  | minimap2 -t "$THREADS" -y -ax map-ont --eqx -k 17 -K "$MM_BATCH" -R "$RG" "$REF_MMI" - 2> "$OUT/logs/minimap2.stderr" \
  | samtools sort -@ "$SORT_THREADS" -m "$SORT_MEM" -T "$T/sort" --reference "$REF_FA" -O cram -o "$CRAM" - 2> "$OUT/logs/sort.stderr"
ALIGN_END=$(date +%s)

samtools quickcheck -v "$CRAM"
samtools index -@ "$SORT_THREADS" "$CRAM"
samtools flagstat -@ "$SORT_THREADS" --input-fmt-option reference="$REF_FA" "$CRAM" > "$OUT/$SAMPLE.flagstat.txt"
samtools idxstats "$CRAM" > "$OUT/$SAMPLE.idxstats.tsv"
samtools view -H "$CRAM" > "$T/out_header.txt"
MOD_OUT=$(samtools view -F 0x900 --reference "$REF_FA" "$CRAM" 2>/dev/null | head -n 2000 | grep -c $'\tMM:Z:' || true)

READS_IN=0; i=0
for f in "${INPUTS[@]}"; do
  i=$((i + 1)); n=$(cat "$T/count.$i"); READS_IN=$((READS_IN + n))
  echo -e "input_reads\t$(basename "$f")\t$n\t$(cat "$T/in.$i.md5")" >> "$OUT/provenance.tsv"
done
{
  echo -e "align_wall_s\t$((ALIGN_END - START))"; echo -e "end\t$(date -Is)"
  echo -e "reads_in\t$READS_IN"; echo -e "cram_md5\t$(md5sum "$CRAM" | awk '{print $1}')"
  echo -e "cram_bytes\t$(stat -c %s "$CRAM")"; echo -e "minimap2_peak\t$(grep -o 'Peak RSS: [0-9.]* GB' "$OUT/logs/minimap2.stderr" | tail -1)"
} >> "$OUT/provenance.tsv"

$TOOL verify --sample "$SAMPLE" --reads-in "$READS_IN" --flagstat "$OUT/$SAMPLE.flagstat.txt" \
  --header "$T/out_header.txt" --mod-in "$MOD_IN" --mod-out "$MOD_OUT" \
  --provenance "$OUT/provenance.tsv" --out "$OUT/intake.json"
echo "intake $SAMPLE OK: $READS_IN reads, $((ALIGN_END - START)) s alignment"
