#!/usr/bin/env bash
# =============================================================================
# Sharded panel build on the synthetic dataset (docs/PANEL_BUILD.md), checked
# against the unsharded panel of run_synthetic.sh.
#
#   bash tests/synthetic/run_panel_build_synthetic.sh WORKDIR
#
# Needs WORKDIR from a finished tests/synthetic/run_synthetic.sh (data, site
# overlay, unsharded panel release). Runs, without a container (CXSV_SIF=none):
#   1. make_shards.py over a two-node synthetic fleet (nodeA, nodeB, both this
#      machine) -> plan, node sheets and overlays; the plan must be a partition
#   2. run_node.py nodeA and nodeB at the same time: staging, one snakemake
#      instance per sample, per-node ledgers, per-sample cost records
#   3. a second run_node.py nodeA must change nothing (resume = continue)
#   4. build_status.py --require-complete, then aggregate.sh over both shards
#   5. the sharded release's junction sites, carriers and callability must be
#      identical to the unsharded panel's
# =============================================================================
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
W="$(cd "$1" && pwd)"
PB="$REPO/scripts/panel_build"
B="$W/panel_build"
CORES="${CORES:-4}"
[ -f "$W/results/panel/release/panel_junction_sites.tsv" ] || { echo "run tests/synthetic/run_synthetic.sh $W first"; exit 1; }
rm -rf "$B"; mkdir -p "$B"
export CXSV_SIF=none

cat > "$B/fleet.yaml" <<YAML
build_root: "$B"
sif: "none"
expected_fingerprint: "unused"
nas_aggregate_mb_s: 500
max_concurrent_stage: 1
max_attempts: 2
nodes:
  nodeA:
    ssh: local
    threads: 64
    ram_gb: 8
    reserve_threads: 0
    reserve_ram_gb: 0
    slots: 1
    cores_per_sample: $CORES
    sv_threads: 2
    speed: 1.0
    nic_mb_s: 500
    stage_dir: "$B/local/nodeA/stage"
    tmp_dir: "$B/local/nodeA/tmp"
    stage_max_gb: 10
    active_hours: ""
    enabled: true
  nodeB:
    ssh: local
    threads: 64
    ram_gb: 8
    reserve_threads: 0
    reserve_ram_gb: 0
    slots: 2
    cores_per_sample: $CORES
    sv_threads: 2
    speed: 1.0
    nic_mb_s: 500
    stage_dir: "$B/local/nodeB/stage"
    tmp_dir: "$B/local/nodeB/tmp"
    stage_max_gb: 10
    active_hours: ""
    enabled: true
YAML
printf 'group\ttest_samples\ttest_coverage_x\tcpu_h_per_x\twall_h_per_x\tbusy_threads\tpeak_rss_gb\tkept_bytes_per_x\tstage_mb_s\tref_speed\n' > "$B/cost_model.tsv"
printf 'ONT\tsynthetic\t1\t0.01\t0.01\t1\t0.0001\t1000\t500\t1.0\n' >> "$B/cost_model.tsv"
# the unsharded run's site overlay, results under the build root instead
sed "s#^panel_results_root:.*#panel_results_root: \"$B/results/panel\"#" "$W/site.yaml" > "$B/site.yaml"
# synthetic sheets have no coverage_x: give every sample 1x for the planner
awk -F'\t' -v OFS='\t' 'NR==1{print $0, "coverage_x"; next} {print $0, 1}' "$W/panel_sheet.tsv" > "$B/panel_sheet.tsv"
sed -i.bak "s#$W/panel_sheet.tsv#$B/panel_sheet.tsv#" "$B/site.yaml"

echo "== 1. plan"
python3 "$PB/make_shards.py" --site "$B/site.yaml" --fleet "$B/fleet.yaml" --cost-model "$B/cost_model.tsv" \
    --build "$B" --write
mkdir -p "$B/code" && ln -sfn "$REPO/workflow" "$B/code/workflow" && ln -sfn "$REPO/config" "$B/code/config" \
    && ln -sfn "$REPO/scripts" "$B/code/scripts"
python3 - "$B" <<'PY'
import csv, sys
b = sys.argv[1]
plan = list(csv.DictReader(open(f"{b}/plan/shard_plan.tsv"), delimiter="\t"))
keys = [(r["cohort"], r["sample_id"]) for r in plan]
assert len(keys) == len(set(keys)) == 3, keys
assert {r["node"] for r in plan} == {"nodeA", "nodeB"}, plan
print("PASS plan is a partition of 3 samples over nodeA, nodeB")
PY

echo "== 2. both nodes at once"
python3 "$PB/run_node.py" "$B" nodeA > "$B/nodeA.runner.log" 2>&1 & PA=$!
python3 "$PB/run_node.py" "$B" nodeB > "$B/nodeB.runner.log" 2>&1 & PB_=$!
wait $PA || { tail -30 "$B/nodeA.runner.log"; exit 1; }
wait $PB_ || { tail -30 "$B/nodeB.runner.log"; exit 1; }
for n in nodeA nodeB; do
    grep -q "exit_runner" "$B/nodes/$n/ledger.tsv" || { echo "FAIL $n ledger has no exit_runner"; exit 1; }
done
ls "$B"/local/*/stage/*/* 2>/dev/null | grep -v runner.lock && { echo "FAIL staged copies left behind"; exit 1; } || true

echo "== 3. resume: a second runner on nodeA must do nothing"
before=$(grep -c started "$B/nodes/nodeA/ledger.tsv")
python3 "$PB/run_node.py" "$B" nodeA > "$B/nodeA.runner2.log" 2>&1
after=$(grep -c started "$B/nodes/nodeA/ledger.tsv")
[ "$before" = "$after" ] && echo "PASS rerun started no sample" || { echo "FAIL rerun started samples again"; exit 1; }

echo "== 4. status + aggregate"
python3 "$PB/build_status.py" "$B" --write --require-complete --costs
bash "$PB/aggregate.sh" "$B" "$CORES"

echo "== 5. sharded == unsharded"
RC=0
for f in panel_junction_sites.tsv panel_junction_carriers.tsv panel_junction_callability.tsv panel_individuals.tsv; do
    if cmp -s "$W/results/panel/release/$f" "$B/results/panel/release/$f"; then echo "PASS identical $f"
    else echo "FAIL $f differs"; diff <(head -5 "$W/results/panel/release/$f") <(head -5 "$B/results/panel/release/$f") || true; RC=1; fi
done
for s in panel1 panel2 panel3; do
    a="$W/results/panel/samples/synthetic_panel/$s/junctions/junctions.tsv"; b="$B/results/panel/samples/synthetic_panel/$s/junctions/junctions.tsv"
    cmp -s "$a" "$b" && echo "PASS identical junctions.tsv for $s" || { echo "FAIL junctions.tsv differs for $s"; RC=1; }
    [ -s "$B/results/panel/samples/synthetic_panel/$s/benchmarks/run_svim.tsv" ] || { echo "FAIL no benchmark for $s"; RC=1; }
done
exit $RC
