#!/bin/bash
# Stage B launcher (run on freya): shard the panel samples over the fleet by
# weight and run run_node_depth.sh on each node over ssh; waits for all nodes.
# Usage: launch_depth.sh "freya:12 shadowfax:24 keurig:24 tia@athena:10"
set -uo pipefail
here=$(cd "$(dirname "$0")" && pwd)
source "$here/config.sh"
nodes=${1:-"freya:12 shadowfax:24 keurig:24 tia@athena:10"}
mkdir -p "$CAL_ROOT/shards" "$CAL_ROOT/logs"
python3 - "$SAMPLE_MANIFEST" "$SAMPLE_SHEETS" "$CAL_ROOT/shards" $nodes <<'PY'
import csv, glob, os, sys
man, sheets, out, nodes = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
aln = {}
for f in glob.glob(f"{sheets}/*.tsv"):
    coh = os.path.basename(f)[:-4]
    for r in csv.DictReader(open(f), delimiter="\t"):
        aln[(coh, r["sample_id"])] = r["alignment"]
rows = [(r["cohort"], r["sample_id"]) for r in csv.DictReader(open(man), delimiter="\t")]
missing = [r for r in rows if r not in aln]
assert not missing, missing[:5]
w = [(n.split(":")[0], int(n.split(":")[1])) for n in nodes]
tot = sum(x for _, x in w)
# deal samples round-robin weighted (interleaves cohorts and sizes)
slots = [n for n, x in w for _ in range(x)]
shard = {n: [] for n, _ in w}
for k, r in enumerate(sorted(rows, key=lambda r: (r[1], r[0]))):
    shard[slots[k % tot]].append(r)
for n, rs in shard.items():
    with open(f"{out}/{n.split('@')[-1]}.tsv", "w") as fh:
        for c, s in rs:
            fh.write(f"{c}\t{s}\t{aln[(c, s)]}\n")
    print(n, len(rs))
PY
for nw in $nodes; do
  n=${nw%%:*}; k=${nw##*:}; h=${n##*@}
  if [ "$h" = "$(hostname)" ]; then
    bash "$here/run_node_depth.sh" "$CAL_ROOT/shards/$h.tsv" "$k" > "$CAL_ROOT/logs/depth_$h.log" 2>&1 &
  else
    ssh -o BatchMode=yes "$n" "bash $here/run_node_depth.sh $CAL_ROOT/shards/$h.tsv $k" > "$CAL_ROOT/logs/depth_$h.log" 2>&1 &
  fi
done
wait
n_ok=$(find "$CAL_ROOT/persample" -name depth_catalogues.ok | wc -l)
echo "depth+catalogues done for $n_ok samples"
