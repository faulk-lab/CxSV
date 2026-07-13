#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh — CxSV Pipeline Runner
# =============================================================================
# Runs the Snakemake pipeline against whatever BAMs are already sitting in
# bamfiles/ — your own, local BAMs, or ones fetched via the separate,
# optional scripts/download_1000g_bams.sh helper. This script does NOT
# download anything itself.
#
# Usage:
#   bash scripts/run_pipeline.sh              # run the full pipeline
#   bash scripts/run_pipeline.sh --resume      # resume an interrupted run
#
# Prerequisites:
#   conda activate snakemake_env
#   bamfiles/ must contain at least one *.bam file (any name works)
#     - your own BAMs, copied/symlinked in, or
#     - bash scripts/download_1000g_bams.sh   (optional 1000G-ONT download)
# =============================================================================

set -euo pipefail

# ─── Configuration ────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
DATA_DIR="$PROJECT_DIR/bamfiles"
RESULTS_DIR="$PROJECT_DIR/results"
SNAKEFILE="$PROJECT_DIR/workflow/Snakefile"
CONFIG="$PROJECT_DIR/config/config.yaml"
CORES="${CORES:-32}"

# ─── Argument parsing ────────────────────────────────────────────────────────
RESUME=false

for arg in "$@"; do
    case "$arg" in
        --resume)   RESUME=true ;;
        --help|-h)
            head -20 "$0" | grep '^#' | sed 's/^# \{0,1\}//'
            exit 0
            ;;
        *)
            echo "Unknown option: $arg"
            echo "Usage: bash scripts/run_pipeline.sh [--resume]"
            exit 1
            ;;
    esac
done

# ─── Banner ───────────────────────────────────────────────────────────────────
echo ""
echo "╔══════════════════════════════════════════════════════════════╗"
echo "║           CxSV Long-Read Detection Pipeline                  ║"
echo "╚══════════════════════════════════════════════════════════════╝"
echo ""
echo "Project directory : $PROJECT_DIR"
echo "BAM directory      : $DATA_DIR"
echo "Results directory : $RESULTS_DIR"
echo "Cores             : $CORES"
echo ""

# ─── Prerequisites check ─────────────────────────────────────────────────────
echo "Checking prerequisites..."

if [ ! -f "$CONFIG" ]; then
    echo "ERROR: Config not found: $CONFIG"
    exit 1
fi

if [ ! -f "$SNAKEFILE" ]; then
    echo "ERROR: Snakefile not found: $SNAKEFILE"
    exit 1
fi

N_BAMS=$(ls "$DATA_DIR"/*.bam 2>/dev/null | wc -l || true)
if [ "$N_BAMS" -eq 0 ]; then
    echo ""
    echo "ERROR: No BAMs found in $DATA_DIR"
    echo "Either copy your own .bam file(s) there (any filename works),"
    echo "or run the optional 1000G-ONT download first:"
    echo "  bash scripts/download_1000g_bams.sh"
    exit 1
fi
echo "BAMs found in bamfiles/: $N_BAMS"
echo ""

mkdir -p "$RESULTS_DIR"

# config.yaml and the Snakefile use paths relative to the project root
# (e.g. "bamfiles", "workflow/scripts"), so snakemake must be run from there.
cd "$PROJECT_DIR"

# ─── Dry run ──────────────────────────────────────────────────────────────────
echo "════════════════════════════════════════════════════════════════"
echo "  Dry run (checking DAG)"
echo "════════════════════════════════════════════════════════════════"
snakemake \
    -s "$SNAKEFILE" \
    --configfile "$CONFIG" \
    --cores "$CORES" \
    -n \
    --quiet \
    2>&1 | tail -8
echo ""

# ─── Full run ─────────────────────────────────────────────────────────────────
echo "════════════════════════════════════════════════════════════════"
echo "  Running Snakemake pipeline"
echo "  Samples : $N_BAMS"
echo "  Cores   : $CORES"
echo "════════════════════════════════════════════════════════════════"
echo ""

LOG="$RESULTS_DIR/snakemake_run_$(date +%Y%m%d_%H%M%S).log"
echo "Full log: $LOG"
echo ""

SNAKE_FLAGS=(--cores "$CORES" --configfile "$CONFIG" --keep-going)
if [ "$RESUME" = true ]; then
    SNAKE_FLAGS+=(--rerun-incomplete)
    echo "(Resuming incomplete run)"
fi

snakemake -s "$SNAKEFILE" "${SNAKE_FLAGS[@]}" 2>&1 | tee "$LOG"
SNAKE_EXIT=${PIPESTATUS[0]}

if [ "$SNAKE_EXIT" -ne 0 ]; then
    echo ""
    echo "════════════════════════════════════════════════════════════════"
    echo "  Pipeline finished with errors (exit code $SNAKE_EXIT)"
    echo "  Some samples may have failed — check: $LOG"
    echo "  To resume: bash scripts/run_pipeline.sh --resume"
    echo "════════════════════════════════════════════════════════════════"
    exit "$SNAKE_EXIT"
fi

echo ""
echo "════════════════════════════════════════════════════════════════"
echo "  Pipeline complete!"
echo ""

N_MERGED=$(find "$RESULTS_DIR/per_sample" -name "*_merged_sv.vcf" 2>/dev/null | wc -l)
N_CXSV=$(find "$RESULTS_DIR/per_sample" -name "*_cxsv_only.vcf"  2>/dev/null | wc -l)
echo "  Per-sample merged VCFs : $N_MERGED / $N_BAMS"
echo "  Per-sample CxSV VCFs   : $N_CXSV / $N_BAMS"
echo ""
if [ -f "$RESULTS_DIR/final/cxsv_count.txt" ]; then
    echo "  CxSV summary:"
    cat "$RESULTS_DIR/final/cxsv_count.txt" | sed 's/^/    /'
fi
echo ""
echo "  Outputs:"
echo "    results/final/cxsv_master_table.tsv"
echo "    results/final/long_read_filtered.vcf"
echo "    results/final/cxsv_only.vcf"
echo "    results/plots/*.pdf"
echo "════════════════════════════════════════════════════════════════"
