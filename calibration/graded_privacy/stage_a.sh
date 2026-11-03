#!/bin/bash
# Stage A (one node): families -> per-locus leave-out callability -> leave-out statuses.
# Usage: stage_a.sh [THREADS_CALLABLE] [THREADS_STATUS]
set -uo pipefail
source "$(dirname "$0")/config.sh"
tc=${1:-16}; ts=${2:-8}
mkdir -p "$CAL_ROOT/logs" "$CAL_ROOT/leaveout"
flank=$(python3 -c "import json;print(json.load(open('$RELEASE/release.json'))['callability']['flank_bp'])")
set -e
$X python "$CODE/calibration/graded_privacy/families.py" --panel-individuals "$RELEASE/panel_individuals.tsv" \
   --ped "$REGISTRY_SOURCES/igsr_3202_samples_ped_population.txt" \
   --sample-info "$REGISTRY_SOURCES/igsr_20130606_sample_info.txt" --out "$CAL_ROOT/families.tsv" \
   2> "$CAL_ROOT/logs/families.log"
cat "$CAL_ROOT/logs/families.log"
[ -s "$CAL_ROOT/locus_callable.tsv.gz" ] || /usr/bin/time -v $X python "$CODE/calibration/graded_privacy/locus_callability.py" \
   --sample-manifest "$SAMPLE_MANIFEST" --samples-root "$SAMPLES_ROOT" --callability-dir "$RELEASE/callability" \
   --flank "$flank" --families "$CAL_ROOT/families.tsv" --threads "$tc" --out "$CAL_ROOT/locus_callable.tsv.gz" \
   2> "$CAL_ROOT/logs/locus_callability.log"
tail -3 "$CAL_ROOT/logs/locus_callability.log"
/usr/bin/time -v $X python "$CODE/calibration/graded_privacy/leaveout_status.py" --release-dir "$RELEASE" \
   --sample-manifest "$SAMPLE_MANIFEST" --samples-root "$SAMPLES_ROOT" --families "$CAL_ROOT/families.tsv" \
   --locus-callable "$CAL_ROOT/locus_callable.tsv.gz" --modes "${MODES:-LOIO,LOFO}" --threads "$ts" \
   ${SUBSET:+--samples "$SUBSET"} --out "${STATUS_OUT:-$CAL_ROOT/leaveout}" 2> "$CAL_ROOT/logs/leaveout_status${SUBSET:+_subset}.log"
tail -3 "$CAL_ROOT/logs/leaveout_status${SUBSET:+_subset}.log"
