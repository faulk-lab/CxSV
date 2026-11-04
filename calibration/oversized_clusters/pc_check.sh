#!/bin/bash
# Roadmap item 13, step 4c: the whole-genome positive control as a CHECK of the new
# clustering (never an input to its design). Re-uses a finished patient run of the
# bundle (tests/lab/pos_control_full_genome.sh) by copying its run directory, with its
# .snakemake metadata, so Snakemake's own rerun triggers decide what is re-made with
# this checkout's code bound over the bundle's repo; the panel is the bundle's
# panel-v1.0 release (fingerprint differs by design: allow_fingerprint_mismatch).
#
#   pc_check.sh BUNDLE OLD_PC_WORKDIR NEW_WORKDIR CODE [CORES]
set -euo pipefail
BUNDLE=$1; OLD=$2; W=$3; CODE=$4; CORES=${5:-12}
PC=/mnt/sargasso1/CxSV/datasets/pos_control
ADDON=$PC/addons/E13
source "$BUNDLE/bundle.env"
mkdir -p "$W"
[ -d "$W/run" ] || rsync -a "$OLD/run/" "$W/run/"
printf 'allow_fingerprint_mismatch: true\n' > "$W/run/config/item13_check.yaml"
cp "$OLD/input/pos_control_full.yaml" "$W/run/config/pos_control_full.yaml"
S=(singularity exec --containall --pwd /work --bind "$BUNDLE:/bundle:ro" --bind "$CODE:/bundle/repo:ro"
   --bind "$W/run:/work" --bind "$OLD/input:$OLD/input:ro" --bind "$PC:$PC:ro" --net --network none
   --env XDG_CACHE_HOME=/work/.cache --env CXSV_BUNDLE_ID="$CXSV_BUNDLE_ID" --env CXSV_COMMIT=item13
   "$BUNDLE/container/cxsv-pipeline.sif")
SM=(snakemake -s /bundle/repo/workflow/patient/Snakefile --directory /work --cores "$CORES" --keep-going
    --printshellcmds --configfile /bundle/config/bundle_site.yaml /work/config/pos_control_full.yaml
    /work/config/item13_check.yaml)
"${S[@]}" "${SM[@]}" -n > "$W/dryrun.txt" 2>&1 || true
awk '/^Job stats/{f=1} f{print} f && /^$/{exit}' "$W/dryrun.txt"
up=$(awk '/^Job stats/{f=1} f && $1 ~ /^(run_|filter_|merge_sv|check_cram|callable_mask)/ {s+=$2} END{print s+0}' "$W/dryrun.txt")
[ "$up" = 0 ] || { echo "dry run plans $up caller/merge jobs: stop"; exit 1; }
t=$(date +%s)
"${S[@]}" "${SM[@]}" > "$W/patient_run.log" 2>&1 || { tail -30 "$W/patient_run.log"; exit 1; }
echo "patient rerun $(( $(date +%s) - t )) s"
R=/work/results/patient/samples/pos_control/pos_control
"${S[@]}" python /bundle/repo/tests/positive_control/score_pos_control.py "$PC/release/truth" "$R" \
    /bundle/repo/tests/positive_control/events.yaml /work/score --extra-truth "$ADDON/truth" \
    --known-failures /bundle/repo/tests/positive_control/known_failures.tsv > "$W/score.txt" 2>&1 || true
tail -25 "$W/score.txt"
grep -E '^(this run|panel|status)' "$W/run/results/patient/panel_compatibility.txt" || true
