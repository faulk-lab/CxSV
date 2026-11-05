#!/usr/bin/env bash
# Library QC survey (docs/LIBRARY_QC.md): workflow/scripts/library_qc.py over a task
# list (sample, alignment, merged VCF), PAR at a time, one TSV per task.
#   survey.sh CODE TASKS OUT SIF PAR
# TASKS: tab-separated label, alignment, merged_sv.vcf (no header).
set -uo pipefail
CODE=$1 TASKS=$2 OUT=$3 SIF=$4 PAR=${5:-12}
REF=/mnt/sargasso1/CxSV/refs/hs1.fa RC=/mnt/sargasso1/CxSV/refs/ref_cache
mkdir -p "$OUT"; EXTRA_BIND=""; [ -d /home/minknow/scratch ] && EXTRA_BIND="-B /home/minknow/scratch"
one() {
    IFS=$'\t' read -r lab aln vcf <<<"$1"
    [ -s "$OUT/$lab.tsv" ] && return 0
    singularity exec --containall --cleanenv -B /mnt/sargasso1 $EXTRA_BIND \
        --env REF_PATH="$RC/%2s/%2s/%s" --env REF_CACHE="$RC/%2s/%2s/%s" "$SIF" \
        python3 "$CODE/workflow/scripts/library_qc.py" --aln "$aln" --ref "$REF" \
        --windows "$CODE/config/library_qc_windows.hs1.bed" --merged-vcf "$vcf" --sample "$lab" \
        --out "$OUT/$lab.tsv.tmp" 2> "$OUT/$lab.log" && mv "$OUT/$lab.tsv.tmp" "$OUT/$lab.tsv" || echo "FAILED $lab" >&2
}
export -f one; export OUT SIF CODE REF RC EXTRA_BIND
xargs -a "$TASKS" -d '\n' -P "$PAR" -I{} bash -c 'one "$@"' _ {}
