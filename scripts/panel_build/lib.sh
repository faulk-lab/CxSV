# shellcheck shell=bash
# =============================================================================
# scripts/panel_build/lib.sh — shared helpers for the panel build scripts
# (sourced, not executed). docs/PANEL_BUILD.md describes the whole procedure.
#
#   CXSV_SIF       pinned pipeline image (default: the validated bundle's SIF);
#                  "none" runs the tools from the current environment instead
#                  (the synthetic CI test, tests/synthetic/run_panel_build_synthetic.sh)
#   cxsv_container WORKDIR LOCAL_TMP CMD...   run CMD in the image:
#       --containall --cleanenv (no host $HOME or environment), the NAS bound at
#       the same path (/mnt/sargasso1), LOCAL_TMP (on the node's own disk, never
#       a RAM tmpfs) bound as /tmp, working directory WORKDIR.
#   CXSV_EXTRA_BINDS  space-separated host directories also bound at the same
#       path (e.g. the node-local staging directory holding alignment copies)
#   cxsv_snakemake WORKDIR LOCAL_TMP CODE CORES CONFIG... -- [snakemake args]
#       snakemake -s CODE/workflow/panel_build/Snakefile in the image.
#   cxsv_stage SRC DEST   copy one file to node-local disk (temp + rename,
#       size checked); prints the seconds it took. CXSV_STAGE_BWLIMIT_MB > 0
#       throttles the copy (rsync --bwlimit) to leave NAS bandwidth to others.
# =============================================================================
: "${CXSV_SIF:=/mnt/sargasso1/CxSV/bundles/cxsv-bundle-20261029-1dcabd9770f7/container/cxsv-pipeline.sif}"
: "${CXSV_NAS:=/mnt/sargasso1}"

cxsv_runtime() {
    command -v singularity 2>/dev/null || command -v apptainer 2>/dev/null \
        || { echo "no singularity/apptainer on $(hostname)" >&2; return 1; }
}

cxsv_container() {
    local wd="$1" tmp="$2"; shift 2
    mkdir -p "$tmp" "$wd"
    if [ "$CXSV_SIF" = none ]; then
        ( cd "$wd" && TMPDIR="$tmp" "$@" ); return
    fi
    local rt; rt="$(cxsv_runtime)" || return 1
    local binds=(--bind "$CXSV_NAS" --bind "$tmp:/tmp") b
    for b in ${CXSV_EXTRA_BINDS:-}; do binds+=(--bind "$b"); done
    "$rt" exec --containall --cleanenv "${binds[@]}" --pwd "$wd" \
        --env TMPDIR=/tmp --env XDG_CACHE_HOME="$wd/.cache" \
        "$CXSV_SIF" "$@"
}

# cxsv_snakemake WORKDIR LOCAL_TMP CODE CORES CONFIG [CONFIG...] [-- EXTRA...]
cxsv_snakemake() {
    local wd="$1" tmp="$2" code="$3" cores="$4"; shift 4
    local cfgs=() extra=()
    while [ $# -gt 0 ]; do
        if [ "$1" = "--" ]; then shift; extra=("$@"); break; fi
        cfgs+=("$1"); shift
    done
    cxsv_container "$wd" "$tmp" snakemake -s "$code/workflow/panel_build/Snakefile" \
        --directory "$wd" --cores "$cores" --configfile "${cfgs[@]}" \
        --nolock --keep-going --rerun-incomplete --printshellcmds "${extra[@]}"
}

# cxsv_stage SRC DEST: one NAS -> local copy, atomically; echoes elapsed seconds
cxsv_stage() {
    local src="$1" dst="$2" t0 t1
    mkdir -p "$(dirname "$dst")"
    if [ -f "$dst" ] && [ "$(stat -c %s "$dst")" = "$(stat -c %s "$src")" ]; then echo 0; return 0; fi
    t0=$(date +%s.%N)
    if [ "${CXSV_STAGE_BWLIMIT_MB:-0}" -gt 0 ]; then      # throttled copy, KiB/s
        rsync --bwlimit=$((CXSV_STAGE_BWLIMIT_MB * 1000000 / 1024)) "$src" "$dst.part" \
            && mv "$dst.part" "$dst" || { rm -f "$dst.part"; return 1; }
    else
        cp "$src" "$dst.part" && mv "$dst.part" "$dst" || { rm -f "$dst.part"; return 1; }
    fi
    t1=$(date +%s.%N)
    [ "$(stat -c %s "$dst")" = "$(stat -c %s "$src")" ] || { echo "size mismatch: $dst" >&2; rm -f "$dst"; return 1; }
    awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.1f\n", b-a}'
}
