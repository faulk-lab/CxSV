#!/usr/bin/env bash
# =============================================================================
# launch.sh — start (or continue) the panel build on fleet nodes. Run on freya.
#
#   bash scripts/panel_build/launch.sh BUILD NODE [NODE...]     e.g. ... shadowfax keurig
#   bash scripts/panel_build/launch.sh BUILD all                 every node with a shard
#
# Per node: removes a STOP file left by an earlier stop, then starts
# BUILD/code/scripts/panel_build/run_node.py BUILD NODE detached (setsid nohup)
# on that node, logging to BUILD/nodes/NODE/logs/runner.log. Re-running is
# safe: a node whose runner is alive refuses a second one (flock on the node's
# stage_dir), and a new runner continues where the last one stopped.
# =============================================================================
set -euo pipefail
BUILD="$1"; shift
[ -f "$BUILD/plan/shard_plan.tsv" ] || { echo "no plan in $BUILD (prepare_build.sh)"; exit 1; }
NODES=("$@"); [ "${NODES[0]:-}" = all ] && mapfile -t NODES < <(ls "$BUILD/nodes")
for n in "${NODES[@]}"; do
    nd="$BUILD/nodes/$n"; [ -f "$nd/runner.json" ] || { echo "$n: no shard"; continue; }
    rm -f "$nd/STOP"; mkdir -p "$nd/logs"
    cmd="setsid nohup python3 '$BUILD/code/scripts/panel_build/run_node.py' '$BUILD' '$n' >> '$nd/logs/runner.log' 2>&1 < /dev/null &"
    ssh_alias=$(python3 -c "import json; print(json.load(open('$BUILD/plan/plan.json'))['fleet']['nodes']['$n']['ssh'])")
    if [ "$ssh_alias" = local ]; then bash -c "$cmd"; else ssh -o BatchMode=yes "$ssh_alias" "$cmd"; fi
    echo "$(date -Is) launched runner on $n"
done
