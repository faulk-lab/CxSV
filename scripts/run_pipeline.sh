#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh — run the CxSV panel builder or the patient scorer
# =============================================================================
# Usage:
#   bash scripts/run_pipeline.sh panel   [--extra-config FILE] [--resume] [--dry-run]
#   bash scripts/run_pipeline.sh patient [--extra-config FILE] [--resume] [--dry-run]
#
#   panel    build the background CxSV panel and its release
#            (workflow/panel/Snakefile; config/calling.yaml + config/panel.yaml)
#   patient  call one patient run with the panel's calling module, check the
#            calling fingerprint against the panel release, label CxSVs
#            public/private (workflow/patient/Snakefile; config/calling.yaml +
#            config/patient.yaml)
#
#   --extra-config FILE   site overlay (paths, cohorts, patient sheet); repeatable
#   --resume              re-run incomplete jobs from an interrupted run
#   --dry-run             resolve the DAG and list jobs, run nothing
#   CORES=<n>             cores for snakemake (default 32)
#
# Prerequisite: conda env create -f envs/cxsv-pipeline.yaml && conda activate cxsv-pipeline
# Run from anywhere; snakemake is started in the project root.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
CORES="${CORES:-32}"

usage() { sed -n '4,22p' "$0" | sed 's/^# \{0,1\}//'; }

[ $# -ge 1 ] || { usage; exit 1; }
WORKFLOW="$1"; shift
case "$WORKFLOW" in
    panel|patient) ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: first argument must be 'panel' or 'patient' (got '$WORKFLOW')"; usage; exit 1 ;;
esac

EXTRA=(); RESUME=false; DRY=false
while [ $# -gt 0 ]; do
    case "$1" in
        --extra-config) EXTRA+=("$(cd "$(dirname "$2")" && pwd)/$(basename "$2")"); shift 2 ;;
        --resume)       RESUME=true; shift ;;
        --dry-run|-n)   DRY=true; shift ;;
        -h|--help)      usage; exit 0 ;;
        *) echo "Unknown option: $1"; usage; exit 1 ;;
    esac
done

SNAKEFILE="$PROJECT_DIR/workflow/$WORKFLOW/Snakefile"
cd "$PROJECT_DIR"

command -v snakemake >/dev/null || { echo "ERROR: snakemake not on PATH; activate the cxsv-pipeline env"; exit 1; }

FLAGS=(-s "$SNAKEFILE" --cores "$CORES" --keep-going --printshellcmds)
[ ${#EXTRA[@]} -gt 0 ] && FLAGS+=(--configfile "${EXTRA[@]}")
$RESUME && FLAGS+=(--rerun-incomplete)

echo "CxSV $WORKFLOW workflow"
echo "  project : $PROJECT_DIR"
echo "  config  : config/calling.yaml config/$WORKFLOW.yaml ${EXTRA[*]:-}"
echo "  cores   : $CORES"

if $DRY; then
    snakemake "${FLAGS[@]}" --dry-run --quiet rules
    exit 0
fi

mkdir -p "$PROJECT_DIR/results/logs"
LOG="$PROJECT_DIR/results/logs/${WORKFLOW}_$(date +%Y%m%d_%H%M%S).log"
echo "  log     : $LOG"
set +e
snakemake "${FLAGS[@]}" 2>&1 | tee "$LOG"
RC=${PIPESTATUS[0]}
set -e
if [ "$RC" -ne 0 ]; then
    echo "Finished with errors (exit $RC). Resume with: bash scripts/run_pipeline.sh $WORKFLOW ${EXTRA[*]:+--extra-config ${EXTRA[*]}} --resume"
    exit "$RC"
fi
echo "Done."
