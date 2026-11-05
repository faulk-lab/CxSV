#!/usr/bin/env bash
# =============================================================================
# check_callable_identity.sh - is a held panel-v1.0 callability mask what the
# rebuild's code would produce? (panel-v1.1 rebuild, roadmap item 13)
#
#   check_callable_identity.sh CODE OLD_RESULTS OUT SIF THREADS PAR COHORT:SAMPLE:ALN ...
#
# Runs CODE/workflow/scripts/callable_mask.sh (the rule callable_mask script of
# the rebuild commit) in the pinned image on each listed alignment, with the
# rule's parameters (config/panel.yaml panel_callability, calling.yaml
# sv_params.min_mapq), and compares callable.bed.gz and
# mosdepth.mosdepth.summary.txt byte for byte with OLD_RESULTS/samples/<cohort>/
# <sample>/callability/. Writes OUT/identity.tsv (one row per genome) and exits
# 1 unless every genome is identical.
# =============================================================================
set -euo pipefail
CODE="$1" OLD="$2" OUT="$3" SIF="$4" THREADS="$5" PAR="$6"; shift 6
REF=/mnt/sargasso1/CxSV/refs/hs1.fa; RC=/mnt/sargasso1/CxSV/refs/ref_cache
read -r MIN_DEPTH CRAM31 MAPQ < <(singularity exec --containall --bind /mnt/sargasso1 "$SIF" python3 - "$CODE" <<'PY'
import sys, yaml
c = sys.argv[1]
p = yaml.safe_load(open(f"{c}/config/panel.yaml"))["panel_callability"]
s = yaml.safe_load(open(f"{c}/config/calling.yaml"))["sv_params"]
print(p["min_depth"], p.get("cram31", "transcode"), s["min_mapq"])
PY
)
mkdir -p "$OUT"
one() {
    IFS=: read -r cohort sample aln <<<"$1"
    local d="$OUT/$cohort/$sample" t0 t1 v
    mkdir -p "$d/tmp"; t0=$(date +%s)
    singularity exec --containall --cleanenv --bind /mnt/sargasso1 --bind "$OUT" --bind "$d/tmp:/tmp" \
        --env REF_PATH="$RC/%2s/%2s/%s" --env REF_CACHE="$RC/%2s/%2s/%s" "$SIF" \
        bash "$CODE/workflow/scripts/callable_mask.sh" "$aln" "$REF" "$d/mosdepth" "$d/callable.bed.gz" \
        "$MIN_DEPTH" "$MAPQ" "$THREADS" "$CRAM31" > "$d/callable_mask.log" 2>&1 || { echo -e "$cohort\t$sample\tFAILED\t-\t-\t-" > "$d/row.tsv"; return 0; }
    t1=$(date +%s)
    local o="$OLD/samples/$cohort/$sample/callability"
    local m s; cmp -s "$d/callable.bed.gz" "$o/callable.bed.gz" && m=identical || m=DIFFERENT
    cmp -s "$d/mosdepth.mosdepth.summary.txt" "$o/mosdepth.mosdepth.summary.txt" && s=identical || s=DIFFERENT
    v=$(head -c 6 "$aln" | od -An -tu1 | awk '$1==67{print "CRAM" $5 "." $6; next}{print "BAM"}')
    echo -e "$cohort\t$sample\t$v\t$m\t$s\t$((t1 - t0))\t$(sha256sum "$d/callable.bed.gz" | cut -c1-16)" > "$d/row.tsv"
}
export -f one; export OUT OLD CODE SIF THREADS REF RC MIN_DEPTH MAPQ CRAM31
printf '%s\n' "$@" | xargs -P "$PAR" -I{} bash -c 'one "$@"' _ {}
{ echo -e "cohort\tsample\tinput\tcallable_bed_gz\tmosdepth_summary\twall_s\tsha256_16"
  for a in "$@"; do IFS=: read -r c s _ <<<"$a"; cat "$OUT/$c/$s/row.tsv"; done; } > "$OUT/identity.tsv"
cat "$OUT/identity.tsv"
! grep -qE "DIFFERENT|FAILED" "$OUT/identity.tsv"
