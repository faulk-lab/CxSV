#!/usr/bin/env bash
# =============================================================================
# panel-v1.3 design check, step 4: re-run the DOWNSTREAM part of a finished patient run
# (any sample: test genome, stand-in, positive control) with another code checkout,
# keeping the bundle's container and the run's caller outputs. Generalises
# tests/lab/pos_control_rerun_downstream.sh to any patient run.
#
#   rerun_patient_downstream.sh BUNDLE SRC W CODE LABEL SAMPLE ALIGNMENT SHEET [CORES]
#
#   SRC        directory holding the finished run's results/ (only read; rsync -a, mtimes kept)
#   W          new node-local run directory
#   LABEL      patient_runs key of the run; SAMPLE its sample id
#   ALIGNMENT  the run's input alignment (its index next to it); the sheet SHEET is copied
#              with its alignment/index columns pointed at it
#
#   CODE's files must be older than the run's outputs (the frozen calibration copy is
#   touched to 2026-10-01), else every rule that lists its script as input re-runs.
#   1. copy    results/ without .snakemake; caller VCFs checksummed before and after
#   2. dry run snakemake -n --forcerun sample_vcf_to_bed calling_fingerprint cn_params_json
#              --rerun-triggers mtime, target = scoring/cxsv_scores.tsv (the clinician-report
#              rule is not a target: the v1.2 image has no python-docx/reportlab). The plan
#              must not contain a caller, caller-filter or merge rule
#   3. run     the same in the bundle's SIF (--containall, network off) with
#              allow_fingerprint_mismatch: true (branch calling fingerprint != release)
# =============================================================================
set -euo pipefail
BUNDLE=$1; SRC=$2; W=$3; CODE=$4; LABEL=$5; SAMPLE=$6; ALN=$7; SHEET=$8; CORES=${9:-12}
source "$BUNDLE/bundle.env"
SIF="$BUNDLE/container/cxsv-pipeline.sif"
RT=$(command -v singularity || command -v apptainer)
S_REL=results/patient/samples/$LABEL/$SAMPLE
COMMIT=$(git -C "$CODE" rev-parse HEAD 2>/dev/null || cat "$CODE/COMMIT_USED" 2>/dev/null || echo unknown)
hdr() { echo; echo "== $* ($(date '+%F %T'))"; }
echo "bundle $CXSV_BUNDLE_ID, code $CODE at $COMMIT, source $SRC, sample $LABEL/$SAMPLE"

hdr "1. copy"
mkdir -p "$W/run/config" "$W/run/tmp"
rsync -a --exclude '.snakemake' --exclude '.cache' --exclude 'tmp' "$SRC/results" "$W/run/"
( cd "$W/run/$S_REL" && sha256sum callers/*.vcf merged_sv.vcf > "$W/upstream.sha256" )
ALNDIR=$(dirname "$ALN")
python3 - "$SHEET" "$ALN" "$W/run/sheet.tsv" <<'PY'
import csv, os, sys
src, aln, out = sys.argv[1:4]
rows = list(csv.DictReader(open(src), delimiter="\t"))
assert len(rows) == 1, rows
r = rows[0]; r["alignment"] = aln
if "index" in r:
    r["index"] = next(aln + x for x in (".crai", ".bai", ".csi") if os.path.exists(aln + x))
w = csv.DictWriter(open(out, "w"), fieldnames=list(r), delimiter="\t"); w.writeheader(); w.writerow(r)
PY
cat > "$W/run/config/run.yaml" <<YAML
patient_runs:
  $LABEL: {samples: "/work/sheet.tsv"}
threads: {samtools: 4, sv_callers: $CORES}
allow_fingerprint_mismatch: true
YAML
X=("$RT" exec --containall --net --network none --pwd /work
   --bind "$BUNDLE:/bundle:ro" --bind "$W/run:/work" --bind "$CODE:/code:ro" --bind "$ALNDIR:$ALNDIR:ro"
   --env "CXSV_BUNDLE_ID=$CXSV_BUNDLE_ID,CXSV_COMMIT=$COMMIT,CXSV_IMAGE_ID=$CXSV_IMAGE_ID"
   --env "CXSV_SIF_SHA256=$CXSV_SIF_SHA256,CXSV_RUNTIME=$(basename "$RT"),XDG_CACHE_HOME=/work/.cache,TMPDIR=/work/tmp"
   "$SIF")
# Only the rules downstream of the callers may run (--allowed-rules): callers, caller filters,
# the merge, cn_depth and library_qc keep their copied outputs, whatever their mtimes (a NAS
# copy can show a caller's raw VCF as newer than its filtered VCF). Nothing is touched.
DOWNSTREAM_RULES="sample_vcf_to_bed sample_cluster_sv sample_classify_cxsv sample_extract_cxsv_vcf
  cn_params_json cn_segment sample_cn_linkage sample_annotate_cxsv sample_junctions calling_fingerprint
  verify_panel_release check_panel_compatibility classify_patient_cxsv catalogue_hits literature_hits
  locus_depth scoring_params score_cxsv"
SM=(snakemake "/work/$S_REL/scoring/cxsv_scores.tsv" "/work/results/patient/panel_compatibility.txt"
    -s /code/workflow/patient/Snakefile --directory /work --cores "$CORES"
    --rerun-triggers mtime --keep-going --printshellcmds
    --forcerun sample_vcf_to_bed calling_fingerprint cn_params_json scoring_params
    --allowed-rules $DOWNSTREAM_RULES
    --configfile /bundle/config/bundle_site.yaml /work/config/run.yaml)

hdr "2. dry run"
"${X[@]}" "${SM[@]}" -n > "$W/dryrun.txt" 2>&1 || { tail -30 "$W/dryrun.txt"; exit 1; }
grep -E '^(localrule|rule) ' "$W/dryrun.txt" | sort | uniq -c
if grep -E '^(localrule|rule) (check_cram_reference|run_sniffles|run_cutesv|run_svim|filter_|merge_sv_callers)' "$W/dryrun.txt"; then
    echo "the plan re-runs a caller rule; stopping"; exit 1
fi
[ -n "${DRY_ONLY:-}" ] && { echo "DRY_ONLY set: stopping after the dry run"; exit 0; }

hdr "3. downstream run with $COMMIT"
t=$(date +%s)
"${X[@]}" "${SM[@]}" > "$W/rerun.log" 2>&1 || { tail -40 "$W/rerun.log"; exit 1; }
echo "downstream run finished: $(( $(date +%s)-t )) s"
( cd "$W/run/$S_REL" && sha256sum -c --quiet "$W/upstream.sha256" ) && echo "caller VCFs unchanged"
grep -E '^(this run|panel|status)' "$W/run/results/patient/panel_compatibility.txt" | sed 's/^/  /' || true
echo done
