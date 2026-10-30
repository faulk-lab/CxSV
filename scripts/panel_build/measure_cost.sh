#!/usr/bin/env bash
# =============================================================================
# measure_cost.sh — per-rule cost of the panel per-sample path on test genomes
#
#   bash scripts/panel_build/measure_cost.sh COST_ROOT CODE CORES THREADS STAGE KEY [KEY...]
#
#   COST_ROOT  NAS directory for this measurement (results/, work/, logs/, sheets/)
#   CODE       repository checkout to run (workflow/panel_build/Snakefile)
#   CORES      snakemake --cores per sample instance
#   THREADS    threads.sv_callers (Sniffles2, cuteSV); SVIM is single-threaded
#   STAGE      node-local directory for alignment copies, or "none" to read
#              the alignments from the NAS in place
#   KEY        COHORT/SAMPLE of the site sheets (config/sites/faulk_lab_panel.yaml)
#
# Per KEY, in order: copy the alignment + index to STAGE (timed; the build's
# node runner does the same), write a one-sample overlay, then start one
# snakemake instance (workflow/panel_build/Snakefile) in the background; the
# next KEY is staged while earlier ones run. Every per-sample rule gets a
# benchmark file and the sample a cost record ({S}/panel_build/cost.tsv);
# /usr/bin/time -v around each instance gives its total CPU time as a
# cross-check. The NIC receive counter is logged every 60 s (logs/nic_rx.tsv).
# Outputs: COST_ROOT/logs/<sample>.{stage.tsv,snakemake.log,time.txt,exit}
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib.sh
source "$HERE/lib.sh"
ROOT="$1"; CODE="$2"; CORES="$3"; THREADS="$4"; STAGE="$5"; shift 5
SITE="$CODE/config/sites/faulk_lab_panel.yaml"
TMPBASE="${CXSV_LOCAL_TMP:-$HOME/scratch/cxsv_tmp}"
mkdir -p "$ROOT/logs" "$ROOT/work" "$ROOT/overlays"
[ "$STAGE" != none ] && export CXSV_EXTRA_BINDS="$STAGE" && mkdir -p "$STAGE"

NIC=$(ip route get 192.168.0.169 | grep -oP 'dev \K\S+')
( while true; do echo -e "$(date +%s)\t$(cat /sys/class/net/$NIC/statistics/rx_bytes)"; sleep 60; done ) \
    > "$ROOT/logs/nic_rx.tsv" & NICPID=$!
trap 'kill $NICPID 2>/dev/null || true' EXIT

sheet_col() {  # cohort sample column -> value from the site sheet
    awk -F'\t' -v s="$2" -v col="$3" 'NR==1{for(i=1;i<=NF;i++) if($i==col) c=i; next} $1==s{print $c}' \
        "/mnt/sargasso1/CxSV/datasets/registry/sample_sheets/$1.tsv"
}

pids=()
for key in "$@"; do
    cohort="${key%%/*}"; s="${key#*/}"
    stage_arg=()
    if [ "$STAGE" != none ]; then
        aln=$(sheet_col "$cohort" "$s" alignment); idx=$(sheet_col "$cohort" "$s" index)
        ta=$(cxsv_stage "$aln" "$STAGE/$cohort/$(basename "$aln")")
        ti=$(cxsv_stage "$idx" "$STAGE/$cohort/$(basename "$idx")")
        echo -e "key\tbytes\tstage_s\n$key\t$(stat -c %s "$aln")\t$ta" > "$ROOT/logs/$s.stage.tsv"
        echo "$(date -Is) staged $key: $(stat -c %s "$aln") bytes in ${ta}s (+ index ${ti}s)"
        stage_arg=(--stage-dir "$STAGE")
    fi
    cxsv_container "$ROOT/work" "$TMPBASE/$s" python3 "$HERE/cost_overlay.py" --site "$SITE" \
        --root "$ROOT" --threads "$THREADS" --tag "$s" "${stage_arg[@]}" "$key" > "$ROOT/overlays/$s.yaml"
    ( set +e
      /usr/bin/time -v -o "$ROOT/logs/$s.time.txt" \
          bash -c "source '$HERE/lib.sh'; cxsv_snakemake '$ROOT/work' '$TMPBASE/$s' '$CODE' '$CORES' '$ROOT/overlays/$s.yaml'" \
          > "$ROOT/logs/$s.snakemake.log" 2>&1
      echo "$key exit=$? end=$(date -Is)" > "$ROOT/logs/$s.exit" ) &
    pids+=($!)
    echo "$(date -Is) started $key (pid $!)"
done
for p in "${pids[@]}"; do wait "$p" || true; done
echo "$(date -Is) all instances finished"; cat "$ROOT"/logs/*.exit
