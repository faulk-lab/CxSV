#!/usr/bin/env bash
# =============================================================================
# download_1000g_bams.sh — OPTIONAL 1000G-ONT BAM downloader
# =============================================================================
# Downloads long-read BAMs from the public 1000G-ONT S3 bucket into
# bamfiles/<cohort>/, saved as {sample}.bam — the pipeline itself has no
# naming requirement, so this is just this script's own local convention.
# This download is entirely optional and separate from the pipeline itself —
# if you already have your own local BAMs, just copy/symlink them into
# bamfiles/<cohort>/ (any filename works) and skip this script, then run
# scripts/run_pipeline.sh --cohort <cohort> directly.
#
# S3 bucket settings are documented in config/remote_1000g.yaml.
#
# Usage:
#   bash scripts/download_1000g_bams.sh                     # download into bamfiles/1000g/
#   bash scripts/download_1000g_bams.sh --cohort <name>      # download into bamfiles/<name>/
#   bash scripts/download_1000g_bams.sh --verify-only        # just check what's present
#
# Prerequisites:
#   conda install -c conda-forge awscli     # for fast parallel downloads (optional)
#   config/samples.txt must exist           # generate with:
#                                            #   python scripts/fetch_sample_list.py
#
# Estimated size: ~118 GB per BAM (100 samples = ~11.8 TB total)
# =============================================================================

set -euo pipefail

# ─── Configuration ────────────────────────────────────────────────────────────
# Kept in sync with config/remote_1000g.yaml — that file is the documented
# source of truth for these values; this script hardcodes them for simplicity
# since it has no YAML parser dependency.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
SAMPLES_FILE="$PROJECT_DIR/config/samples.txt"
COHORT="1000g"

S3_BUCKET="1000g-ont"
S3_PREFIX="ALIGNMENT_AND_ASSEMBLY_DATA/FIRST_100/NAPU_PIPELINE/CHM13"
S3_BASE="https://s3.amazonaws.com/${S3_BUCKET}"
# This is how the bucket names files server-side — unrelated to the pipeline's
# (lack of a) local naming convention. Downloaded BAMs are saved locally as
# plain {sample}.bam regardless.
REMOTE_BAM_SUFFIX="_PMDV_FINAL.haplotagged.bam"

# Number of parallel downloads (adjust based on your internet connection)
PARALLEL_DOWNLOADS="${PARALLEL_DOWNLOADS:-4}"

# ─── Argument parsing ────────────────────────────────────────────────────────
VERIFY_ONLY=false
while [ $# -gt 0 ]; do
    case "$1" in
        --cohort)      COHORT="$2"; shift 2 ;;
        --verify-only) VERIFY_ONLY=true; shift ;;
        --help|-h)
            head -22 "$0" | grep '^#' | sed 's/^# \{0,1\}//'
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: bash scripts/download_1000g_bams.sh [--cohort <name>] [--verify-only]"
            exit 1
            ;;
    esac
done

DATA_DIR="$PROJECT_DIR/bamfiles/$COHORT"

echo ""
echo "╔══════════════════════════════════════════════════════════════╗"
echo "║       1000G-ONT BAM Download (optional)                      ║"
echo "╚══════════════════════════════════════════════════════════════╝"
echo ""
echo "Cohort            : $COHORT"
echo "Target directory  : $DATA_DIR"
echo ""

if [ ! -f "$SAMPLES_FILE" ]; then
    echo "ERROR: Sample list not found: $SAMPLES_FILE"
    echo "Generate it by running:"
    echo "  python scripts/fetch_sample_list.py"
    exit 1
fi

SAMPLES=( $(grep -v '^#' "$SAMPLES_FILE" | grep '\S') )
TOTAL=${#SAMPLES[@]}
echo "Samples listed    : $TOTAL"

if [ "$TOTAL" -eq 0 ]; then
    echo "ERROR: No samples found in $SAMPLES_FILE"
    exit 1
fi

mkdir -p "$DATA_DIR"

N_ALREADY=$(ls "$DATA_DIR"/*.bam 2>/dev/null | wc -l || true)
echo "Already present   : $N_ALREADY / $TOTAL"
echo ""

if [ "$VERIFY_ONLY" = true ]; then
    for SAMPLE in "${SAMPLES[@]}"; do
        if [ ! -f "$DATA_DIR/${SAMPLE}.bam" ]; then
            echo "  MISSING: ${SAMPLE}.bam"
        fi
    done
    exit 0
fi

# Check disk space available under bamfiles/<cohort>/
AVAILABLE_GB=$(df -BG "$DATA_DIR" | awk 'NR==2 {gsub("G",""); print $4}')
NEEDED_GB=$(( TOTAL * 118 ))
echo "Disk available    : ~${AVAILABLE_GB} GB"
echo "Disk needed       : ~${NEEDED_GB} GB  (${TOTAL} BAMs × 118 GB)"
echo ""

if [ "$AVAILABLE_GB" -lt "$NEEDED_GB" ]; then
    echo "WARNING: You may not have enough disk space."
    read -r -p "Continue anyway? [y/N] " CONFIRM
    if [[ ! "$CONFIRM" =~ ^[Yy]$ ]]; then
        echo "Aborted."
        exit 1
    fi
fi

echo "════════════════════════════════════════════════════════════════"
echo "  Downloading $TOTAL BAMs from S3"
echo "  Parallel downloads: $PARALLEL_DOWNLOADS"
echo "════════════════════════════════════════════════════════════════"
echo ""

if command -v aws &>/dev/null; then
    USE_AWS=true
    echo "Using: aws cli (recommended — faster, supports resume)"
else
    USE_AWS=false
    echo "Using: curl (install awscli for faster downloads)"
    echo "  conda install -c conda-forge awscli"
fi
echo ""

download_bam() {
    local SAMPLE="$1"
    local BAM_FILE="${DATA_DIR}/${SAMPLE}.bam"
    local S3_KEY="${S3_PREFIX}/${SAMPLE}/${SAMPLE}${REMOTE_BAM_SUFFIX}"
    local BAM_URL="${S3_BASE}/${S3_KEY}"

    if [ -f "$BAM_FILE" ]; then
        echo "  [skip]  $SAMPLE — already present"
        return 0
    fi

    echo "  [start] $SAMPLE"
    if [ "${USE_AWS}" = true ]; then
        aws s3 cp \
            "s3://${S3_BUCKET}/${S3_KEY}" \
            "$BAM_FILE" \
            --no-sign-request \
            --no-progress \
            2>&1 | tail -1
    else
        curl -L \
            --retry 5 \
            --retry-delay 10 \
            --retry-all-errors \
            --continue-at - \
            --silent \
            --show-error \
            "$BAM_URL" \
            -o "$BAM_FILE"
    fi
    echo "  [done]  $SAMPLE — $(du -sh "$BAM_FILE" 2>/dev/null | cut -f1)"
}

export -f download_bam
export DATA_DIR REMOTE_BAM_SUFFIX S3_BUCKET S3_PREFIX S3_BASE USE_AWS

printf '%s\n' "${SAMPLES[@]}" \
    | xargs -P "$PARALLEL_DOWNLOADS" -I{} bash -c 'download_bam "$@"' _ {}

echo ""
N_DOWNLOADED=$(ls "$DATA_DIR"/*.bam 2>/dev/null | wc -l || true)
echo "Download complete: $N_DOWNLOADED / $TOTAL BAMs in bamfiles/$COHORT/"
echo "Total disk usage: $(du -sh "$DATA_DIR" 2>/dev/null | cut -f1)"
echo ""

MISSING=0
for SAMPLE in "${SAMPLES[@]}"; do
    if [ ! -f "$DATA_DIR/${SAMPLE}.bam" ]; then
        echo "  MISSING: ${SAMPLE}.bam"
        MISSING=$(( MISSING + 1 ))
    fi
done

if [ "$MISSING" -gt 0 ]; then
    echo ""
    echo "ERROR: $MISSING BAMs failed to download."
    echo "Re-run this script to retry — downloads resume automatically."
    exit 1
fi

echo "All $TOTAL BAMs verified."
echo ""
echo "Next: bash scripts/run_pipeline.sh --cohort $COHORT"
