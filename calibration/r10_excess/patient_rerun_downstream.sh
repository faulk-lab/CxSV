#!/usr/bin/env bash
# =============================================================================
# patient_rerun_downstream.sh — re-run the DOWNSTREAM part of a finished patient
# run (bundle run_cxsv.sh patient, or tests/lab/pos_control_full_genome.sh) with
# another code checkout and extra config, keeping the run's caller outputs
# (R10-excess diagnosis: candidate discovery rules on the positive control and
# the stand-in HG01868). Generalises tests/lab/pos_control_rerun_downstream.sh.
#
#   bash patient_rerun_downstream.sh BUNDLE SRC_RUN WORKDIR CODE RUN_CONFIG CORES [EXTRA_YAML]
#
#   SRC_RUN     the finished run directory (holds config/ and results/; only read)
#   RUN_CONFIG  its run config inside SRC_RUN/config (e.g. pos_control_full.yaml)
#   EXTRA_YAML  file whose contents are appended to the rerun config (e.g.
#               cxsv_params: {single_caller_bnd_rule: svim_lowq})
# Env: BINDS="dir ..." read-only host directories the run's sample sheets point
# into (alignments). The plan must not contain a caller, filter or merge rule.
# allow_fingerprint_mismatch is set: the discovery code differs from the
# release's, so privacy is exploratory (junction level, as in the source run).
# =============================================================================
set -euo pipefail
BUNDLE="$(cd "$1" && pwd)"; SRC="$(cd "$2" && pwd)"; W="$(mkdir -p "$3" && cd "$3" && pwd)"
CODE="$(cd "$4" && pwd)"; RUNCFG="$5"; CORES="${6:-12}"; EXTRA="${7:-}"
source "$BUNDLE/bundle.env"
SIF="$BUNDLE/container/cxsv-pipeline.sif"; RT=$(command -v singularity || command -v apptainer)
mkdir -p "$W/run"
rsync -a --exclude '.snakemake' --exclude '.cache' --exclude 'tmp' "$SRC/" "$W/run/"
mkdir -p "$W/run/tmp"
find "$W/run/results" -path '*/callers/*.vcf' -o -name merged_sv.vcf | sort | xargs sha256sum > "$W/upstream.sha256"
{ echo "allow_fingerprint_mismatch: true"; echo "literature_set: {root: null}"; [ -n "$EXTRA" ] && cat "$EXTRA"; } \
    > "$W/run/config/rerun_downstream.yaml"
B=(); for d in ${BINDS:-}; do B+=(--bind "$d:$d:ro"); done
X=("$RT" exec --containall --net --network none --pwd /work --bind "$BUNDLE:/bundle:ro" --bind "$W/run:/work"
   --bind "$CODE:/code:ro" "${B[@]}" --env "XDG_CACHE_HOME=/work/.cache,TMPDIR=/work/tmp" "$SIF")
SM=(snakemake -s /code/workflow/patient/Snakefile --directory /work --cores "$CORES" --keep-going --printshellcmds
    --rerun-triggers mtime --forcerun sample_vcf_to_bed calling_fingerprint
    --configfile /bundle/config/bundle_site.yaml "/work/config/$RUNCFG" /work/config/rerun_downstream.yaml)
"${X[@]}" "${SM[@]}" -n > "$W/dryrun.txt" 2>&1 || { tail -30 "$W/dryrun.txt"; exit 1; }
awk '/^Job stats/{f=1} f{print} f && /^$/{exit}' "$W/dryrun.txt"
if grep -E '^(localrule|rule) (run_sniffles|run_cutesv|run_svim|filter_|merge_sv_callers)' "$W/dryrun.txt"; then
    echo "the plan re-runs a caller rule; stopping"; exit 1
fi
t=$(date +%s)
"${X[@]}" "${SM[@]}" > "$W/rerun.log" 2>&1 || { tail -40 "$W/rerun.log"; exit 1; }
echo "downstream run finished: $(( $(date +%s)-t )) s"
sha256sum -c --quiet "$W/upstream.sha256" && echo "caller VCFs and merged VCF unchanged"
