#!/usr/bin/env bash
# =============================================================================
# prepare_build.sh — freeze the code, check the fingerprint, write the plan,
# and dry-run every shard. Run on freya. Launches nothing.
#
#   bash scripts/panel_build/prepare_build.sh COMMIT COST_MODEL [make_shards args...]
#
#   COMMIT      git commit the build runs (every node runs this exact code)
#   COST_MODEL  cost_model.tsv from scripts/panel_build/cost_model.py
#
# Steps (build root: $CXSV_BUILD_ROOT, else build_root of config/sites/faulk_lab_fleet.yaml, e.g.
# /mnt/sargasso1/CxSV/build/panel-v1.0):
#   1. git archive COMMIT -> BUILD/code (read-only copy + COMMIT file); an
#      existing BUILD/code from another commit is an error
#   2. calling fingerprint of BUILD/code in the image (workflow/panel target
#      calling_fingerprint.json in BUILD/aggregate); must start with the
#      fleet's expected_fingerprint
#   3. make_shards.py --write (refuses to overwrite an existing plan)
#   4. snakemake -n per node shard (workflow/panel_build) and over the full
#      panel (workflow/panel, all sheets): job counts to BUILD/plan/dryrun_*.txt
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
COMMIT="$1"; COST="$2"; shift 2
FLEET="$REPO/config/sites/faulk_lab_fleet.yaml"
SITE="$REPO/config/sites/faulk_lab_panel.yaml"
REGISTRY=/mnt/sargasso1/CxSV/datasets/registry/sample_registry.tsv
val() { python3 -c "import sys; sys.path.insert(0, '$HERE'); import shardlib; print(shardlib.load_yaml('$FLEET')['$1'])"; }
BUILD="${CXSV_BUILD_ROOT:-$(val build_root)}"; export CXSV_SIF="$(val sif)"; EXPECT="$(val expected_fingerprint)"
# shellcheck source=lib.sh
source "$HERE/lib.sh"
TMP="$HOME/scratch/cxsv_panel_build/tmp/prepare"

# 1. freeze the code
full=$(git -C "$REPO" rev-parse "$COMMIT^{commit}")
if [ -d "$BUILD/code" ]; then
    have=$(cat "$BUILD/code/COMMIT")
    [ "$have" = "$full" ] || { echo "BUILD/code is $have, not $full: use a new build root"; exit 1; }
else
    mkdir -p "$BUILD/code.part"
    git -C "$REPO" archive "$full" | tar -x -C "$BUILD/code.part"
    echo "$full" > "$BUILD/code.part/COMMIT"
    mv "$BUILD/code.part" "$BUILD/code"
fi
echo "code: $BUILD/code at $full"
cp "$COST" "$BUILD/cost_model.tsv" 2>/dev/null || true

# 3. plan (before 2: the aggregate overlay names the results root)
python3 "$BUILD/code/scripts/panel_build/make_shards.py" --site "$BUILD/code/config/sites/faulk_lab_panel.yaml" \
    --fleet "$BUILD/code/config/sites/faulk_lab_fleet.yaml" --registry "$REGISTRY" \
    --cost-model "$BUILD/cost_model.tsv" --build "$BUILD" --write --json "$BUILD/plan/estimate.json" "$@"

# 2. fingerprint
cxsv_container "$BUILD/aggregate" "$TMP" snakemake -s "$BUILD/code/workflow/panel/Snakefile" \
    --directory "$BUILD/aggregate" --configfile "$BUILD/aggregate.yaml" --cores 2 \
    "$BUILD/results/panel/calling_fingerprint.json" > "$BUILD/plan/fingerprint.log" 2>&1
fp=$(python3 -c "import json; print(json.load(open('$BUILD/results/panel/calling_fingerprint.json'))['sha256'][:12])")
[ "$fp" = "$EXPECT" ] || { echo "calling fingerprint $fp != expected $EXPECT: stop"; exit 1; }
echo "fingerprint: $fp (expected $EXPECT)"

# 4. dry runs
for nd in "$BUILD"/nodes/*/; do
    n=$(basename "$nd")
    cxsv_container "$nd/work" "$TMP" snakemake -s "$BUILD/code/workflow/panel_build/Snakefile" \
        --directory "$nd/work" --configfile "$nd/node.yaml" --cores 8 --nolock -n \
        > "$BUILD/plan/dryrun_$n.txt" 2>&1
    echo "dry run $n: $(grep -E '^total' "$BUILD/plan/dryrun_$n.txt" || echo FAILED)"
done
cxsv_container "$BUILD/aggregate" "$TMP" snakemake -s "$BUILD/code/workflow/panel/Snakefile" \
    --directory "$BUILD/aggregate" --configfile "$BUILD/aggregate.yaml" --cores 16 -n \
    > "$BUILD/plan/dryrun_full_panel.txt" 2>&1
echo "full panel dry run: $(grep -E '^total' "$BUILD/plan/dryrun_full_panel.txt" || echo FAILED)"
