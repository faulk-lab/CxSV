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
#       size checked, source mtime kept); prints the seconds it took.
#       CXSV_STAGE_BWLIMIT_MB > 0 throttles the copy (rsync --bwlimit) to leave
#       NAS bandwidth to others.
#   cxsv_cram30 SRC DEST DEST_INDEX   if the staged copy DEST is CRAM 3.1,
#       rewrite it as CRAM 3.0 in place (samtools in the image; record counts
#       must match) with a new index, both given SRC's mtime, and write
#       DEST.cram30.ok. The pinned mosdepth 0.3.10 is statically linked with an
#       htslib that cannot decode CRAM 3.1: it logs "Slice decode failure" for
#       every contig, exits 0 and leaves an empty callability mask. Needs
#       CXSV_REF (FASTA) and CXSV_REF_CACHE (htslib MD5 cache directory).
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
        rsync -t --bwlimit=$((CXSV_STAGE_BWLIMIT_MB * 1000000 / 1024)) "$src" "$dst.part" \
            && mv "$dst.part" "$dst" || { rm -f "$dst.part"; return 1; }
    else
        cp --preserve=timestamps "$src" "$dst.part" && mv "$dst.part" "$dst" || { rm -f "$dst.part"; return 1; }
    fi
    t1=$(date +%s.%N)
    [ "$(stat -c %s "$dst")" = "$(stat -c %s "$src")" ] || { echo "size mismatch: $dst" >&2; rm -f "$dst"; return 1; }
    awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.1f\n", b-a}'
}

# cxsv_cram_version FILE: "3.1", "3.0", ... for a CRAM, "" otherwise
cxsv_cram_version() {
    local b; b=$(head -c 6 "$1" | od -An -tu1) || return 1
    set -- $b; [ "${1:-}" = 67 ] && [ "${2:-}" = 82 ] && [ "${3:-}" = 65 ] && [ "${4:-}" = 77 ] && echo "$5.$6"
}

# cxsv_cram30 SRC DEST DEST_INDEX (see the header)
cxsv_cram30() {
    local src="$1" dst="$2" idx="$3" v n_in n_out t="${CXSV_TRANSCODE_THREADS:-4}"
    v=$(cxsv_cram_version "$dst") || return 1
    [ "$v" = "3.1" ] || return 0
    : "${CXSV_REF:?CXSV_REF not set}" "${CXSV_REF_CACHE:?CXSV_REF_CACHE not set}"
    local wd; wd=$(dirname "$dst")
    local hts="export REF_PATH='$CXSV_REF_CACHE/%2s/%2s/%s' REF_CACHE='$CXSV_REF_CACHE/%2s/%2s/%s';"
    rm -f "$dst.v30.part" "$dst.v30.part.crai"
    cxsv_container "$wd" "${CXSV_STAGE_TMP:-$wd}" bash -c "$hts samtools view -@ $t -C --output-fmt-option version=3.0 \
        -T '$CXSV_REF' --write-index -o '$dst.v30.part##idx##$dst.v30.part.crai' '$dst'" || { rm -f "$dst.v30.part"*; return 1; }
    n_in=$(cxsv_container "$wd" "${CXSV_STAGE_TMP:-$wd}" bash -c "$hts samtools view -@ $t -c '$dst'") || return 1
    n_out=$(cxsv_container "$wd" "${CXSV_STAGE_TMP:-$wd}" bash -c "$hts samtools view -@ $t -c '$dst.v30.part'") || return 1
    [ -n "$n_in" ] && [ "$n_in" = "$n_out" ] || { echo "record count mismatch: $n_in vs $n_out" >&2; rm -f "$dst.v30.part"*; return 1; }
    mv "$dst.v30.part.crai" "$idx" && mv "$dst.v30.part" "$dst" || return 1
    touch -r "$src" "$dst" "$idx"
    printf 'source\t%s\nsource_bytes\t%s\nsource_cram_version\t3.1\ncram_version\t3.0\nrecords\t%s\n' \
        "$src" "$(stat -c %s "$src")" "$n_in" > "$dst.cram30.ok"
    echo "transcoded CRAM 3.1 -> 3.0: $n_in records"
}
