#!/usr/bin/env bash
# =============================================================================
# Re-run the DOWNSTREAM part of a finished whole-genome positive-control run
# (tests/lab/pos_control_full_genome.sh) with another code checkout, keeping the
# bundle's container and the run's caller outputs (roadmap item 11, detection
# test of L01-L04 with main after item 13).
#
#   bash tests/lab/pos_control_rerun_downstream.sh BUNDLE SRC_WORKDIR WORKDIR CODE [CORES]
#
#   BUNDLE       the bundle the source run used (its SIF and bundle_site.yaml)
#   SRC_WORKDIR  WORKDIR of the finished pos_control_full_genome.sh run (only read)
#   WORKDIR      new run directory; SRC_WORKDIR/run is copied here (rsync -a, mtimes kept)
#   CODE         repository checkout whose workflow/ is run (e.g. main after item 13)
#
#   1. copy       run/ without .snakemake and caches; the copied caller VCFs are
#                 checksummed before and after
#   2. dry run    snakemake -n --forcerun sample_vcf_to_bed calling_fingerprint
#                 --rerun-triggers mtime (the fingerprint is recomputed for CODE):
#                 the plan must not contain a caller, caller-filter or merge rule
#   3. run        the same, in the bundle's SIF (--containall, network off), with
#                 allow_fingerprint_mismatch: true, because CODE's discovery code
#                 (item 13) changes the calling fingerprint and no panel release
#                 has been built with it yet. Callers are not re-run; privacy is
#                 still junction-level against the bundle's panel release, whose
#                 loci were made with the bundle's code (exploratory use only)
#   4. scoring    CODE's score_pos_control.py against the same truth (+ ADDONS)
#   1b. reference  when BUNDLE is not the bundle of the source run (SRC_BUNDLE, default:
#                 the bundle named in SRC_WORKDIR/pos_control_full_genome.log, next to
#                 BUNDLE), both ref_cache/md5_manifest.tsv must be byte-identical; then
#                 the reference check, callers, filters and merge are marked up to date
#                 with `snakemake --touch merged_sv.vcf` (tests/lab/lib/reuse_callers.sh). Their
#                 outputs are not re-made; without this, the newer cache's mtimes would
#                 re-run every caller and the step-2 guard would stop the run. A
#                 differing manifest stops the rerun: use pos_control_full_genome.sh
# ADDONS / EVENTS as for pos_control_full_genome.sh. LITERATURE=DIR (a literature
# set, acquisition/literature) also turns on the literature_match annotation.
# DRY_ONLY=1 stops after step 2.
# STRICT_FINGERPRINT=1 writes allow_fingerprint_mismatch: false (the run stops unless
# CODE's calling fingerprint equals the bundle's panel release).
# =============================================================================
set -euo pipefail
BUNDLE="$(cd "$1" && pwd)"; SRCW="$(cd "$2" && pwd)"; W="$(mkdir -p "$3" && cd "$3" && pwd)"
CODE="$(cd "$4" && pwd)"; CORES="${5:-12}"
PC="${PC:-/mnt/sargasso1/CxSV/datasets/pos_control}"
EVENTS="${EVENTS:-$CODE/tests/positive_control/events.yaml}"
ADDONS="${ADDONS:-$PC/addons/E13}"
source "$BUNDLE/bundle.env"
SIF="$BUNDLE/container/cxsv-pipeline.sif"
RT=$(command -v singularity || command -v apptainer)
S_REL=results/patient/samples/pos_control/pos_control
hdr() { echo; echo "== $* ($(date '+%F %T'))"; }
# shellcheck source=tests/lab/lib/reuse_callers.sh
source "$(dirname "${BASH_SOURCE[0]}")/lib/reuse_callers.sh"

COMMIT=$(git -C "$CODE" rev-parse HEAD 2>/dev/null || echo unknown)
echo "bundle $CXSV_BUNDLE_ID (SIF $CXSV_SIF_SHA256), code $CODE at $COMMIT"
echo "source run $SRCW"

hdr "1. copy the finished run"
mkdir -p "$W/run"
rsync -a --exclude '.snakemake' --exclude '.cache' --exclude 'tmp' "$SRCW/run/" "$W/run/"
mkdir -p "$W/run/tmp"
( cd "$W/run/$S_REL" && sha256sum callers/*.vcf merged_sv.vcf > "$W/upstream.sha256" )
LITBIND=(); LITCFG="literature_set: {root: null}"
if [ -n "${LITERATURE:-}" ]; then
    LITERATURE="$(cd "$LITERATURE" && pwd)"; LITBIND=(--bind "$LITERATURE:/literature:ro")
    LITCFG="literature_set: {root: /literature}"
fi
cat > "$W/run/config/rerun_downstream.yaml" <<YAML
# written by tests/lab/pos_control_rerun_downstream.sh
allow_fingerprint_mismatch: $([ "${STRICT_FINGERPRINT:-0}" = 1 ] && echo false || echo true)
$LITCFG
YAML

X=("$RT" exec --containall --net --network none --pwd /work
   --bind "$BUNDLE:/bundle:ro" --bind "$W/run:/work" --bind "$CODE:/code:ro"
   --bind "$SRCW/input:$SRCW/input:ro" --bind "$PC:$PC:ro" "${LITBIND[@]}"
   --env "CXSV_BUNDLE_ID=$CXSV_BUNDLE_ID,CXSV_COMMIT=$COMMIT,CXSV_IMAGE_ID=$CXSV_IMAGE_ID"
   --env "CXSV_SIF_SHA256=$CXSV_SIF_SHA256,CXSV_RUNTIME=$(basename "$RT"),XDG_CACHE_HOME=/work/.cache,TMPDIR=/work/tmp"
   "$SIF")
SM_BASE=(snakemake -s /code/workflow/patient/Snakefile --directory /work --cores "$CORES"
         --rerun-triggers mtime --configfile /bundle/config/bundle_site.yaml
         /work/config/pos_control_full.yaml /work/config/rerun_downstream.yaml)
SM=("${SM_BASE[@]}" --keep-going --printshellcmds --forcerun sample_vcf_to_bed calling_fingerprint)

hdr "1b. reference of the source run"
SRC_ID=$(awk '$1=="bundle"{print $2; exit}' "$SRCW/pos_control_full_genome.log" 2>/dev/null || true)
SRC_BUNDLE="${SRC_BUNDLE:-${SRC_ID:+$(dirname "$BUNDLE")/$SRC_ID}}"
MANIFEST_REL=reference/ref_cache/md5_manifest.tsv
if [ -z "$SRC_BUNDLE" ] || [ ! -d "$SRC_BUNDLE" ]; then
    echo "source bundle unknown (set SRC_BUNDLE); nothing touched, the step-2 guard decides"
elif [ "$(cd "$SRC_BUNDLE" && pwd)" = "$BUNDLE" ]; then
    echo "source run used this bundle; nothing to touch"
elif same_reference_manifest "$SRC_BUNDLE/$MANIFEST_REL" "$BUNDLE/$MANIFEST_REL"; then
    echo "md5_manifest.tsv identical in $SRC_BUNDLE and $BUNDLE ($(sha256sum < "$BUNDLE/$MANIFEST_REL" | cut -c1-12)):"
    echo "marking ${CXSV_UPSTREAM_RULES[*]} up to date (snakemake --touch)"
    touch_upstream "/work/$S_REL/merged_sv.vcf" "${X[@]}" "${SM_BASE[@]}" > "$W/touch_upstream.txt" 2>&1 || { tail -30 "$W/touch_upstream.txt"; exit 1; }
    echo "  touched: $(sed -nE 's/^(local)?rule ([A-Za-z0-9_]+):$/\2/p' "$W/touch_upstream.txt" | sort -u | tr '\n' ' ')"
else
    echo "md5_manifest.tsv differs between $SRC_BUNDLE and $BUNDLE: the callers must re-run;"
    echo "use tests/lab/pos_control_full_genome.sh instead"; exit 1
fi

hdr "2. dry run"
"${X[@]}" "${SM[@]}" -n > "$W/dryrun.txt" 2>&1 || { tail -30 "$W/dryrun.txt"; exit 1; }
awk '/^Job stats/{f=1} f{print} f && /^$/{exit}' "$W/dryrun.txt"
grep -E '^(localrule|rule) ' "$W/dryrun.txt" | sort | uniq -c
if grep -E '^(localrule|rule) (run_sniffles|run_cutesv|run_svim|filter_|merge_sv_callers)' "$W/dryrun.txt"; then
    echo "the plan re-runs a caller rule; stopping"; exit 1
fi
[ -n "${DRY_ONLY:-}" ] && { echo "DRY_ONLY set: stopping after the dry run"; exit 0; }

hdr "3. downstream run with $COMMIT"
t=$(date +%s)
"${X[@]}" "${SM[@]}" > "$W/rerun.log" 2>&1 || { tail -40 "$W/rerun.log"; exit 1; }
echo "downstream run finished: $(( $(date +%s)-t )) s"
( cd "$W/run/$S_REL" && sha256sum -c --quiet "$W/upstream.sha256" ) && echo "caller VCFs unchanged"
grep -E '^(this run|panel|status)' "$W/run/results/patient/panel_compatibility.txt" | sed 's/^/  /' || true

hdr "4. score against the truth set"
XT=""; for a in $ADDONS; do XT="$XT --extra-truth $a/truth"; done
"${X[@]}" python /code/tests/positive_control/score_pos_control.py "$PC/release/truth" "/work/$S_REL" \
    "$EVENTS" /work/score $XT --known-failures /code/tests/positive_control/known_failures.tsv \
    > "$W/score.txt" 2>&1 || true
mv "$W/run/score".* "$W/" 2>/dev/null || true
tail -25 "$W/score.txt"
echo "done"
