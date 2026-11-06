#!/usr/bin/env bash
# =============================================================================
# cn_depth_build.sh — full-panel depth pass for the copy-number track (roadmap item 15,
# Phase 2). Same pass as Phase 1 (depth_pass.sh: mosdepth 0.3.10, MAPQ >= 20, fast
# mode, 1 kb bins, i.e. 3,117,304 bins on hs1; the 10 kb bins are sums of ten).
#
#   cn_depth_build.sh tasks  ROOT                 write ROOT/tasks.tsv (key, alignment, bytes)
#                                                 for every panel-v1.1-r2 library, largest first
#   cn_depth_build.sh reuse  ROOT PHASE1_DEPTH    hard-link (or copy) the Phase 1 depth of the
#                                                 same alignment where it is identical by
#                                                 construction (same alignment path and size,
#                                                 same image, full bin count); sha256 recorded
#   cn_depth_build.sh worker ROOT [PAR] [THREADS] depth_worker.sh over ROOT/tasks.tsv
#   cn_depth_build.sh ledger ROOT                 ROOT/ledger.tsv: one row per library
#                                                 (status, host, bytes, seconds, MB/s, bins,
#                                                 sha256, source) and a summary on stderr
#
# Layout: ROOT/<cohort>/<sample>/{depth.regions.bed.gz,.csi,depth.mosdepth.summary.txt,
# timing.tsv,depth.sha256,depth.ok[,reused_from]}; claims in ROOT/claims/<cohort>/<sample>
# (atomic mkdir on the share); failures in ROOT/failed. Resumable: a library with
# depth.ok is skipped; a claim left by a killed worker is released by deleting its
# claims/<cohort>/<sample> directory (after checking no process still runs it).
# =============================================================================
set -euo pipefail
MODE=$1; ROOT=$2
HERE="$(cd "$(dirname "$0")" && pwd)"
R2=${R2:-/mnt/sargasso1/CxSV/build/panel-v1.1-r2/results/panel}
SHEETS=${SHEETS:-/mnt/sargasso1/CxSV/datasets/registry/sample_sheets}
FAI=${FAI:-/mnt/sargasso1/CxSV/refs/hs1.fa.fai}
mkdir -p "$ROOT"
case "$MODE" in
tasks)
    python3 - "$R2/panel/sample_manifest.tsv" "$SHEETS" "$ROOT/tasks.tsv" <<'PY'
import csv, os, sys
man, sheets, out = sys.argv[1:]
rows = []
cache = {}
for r in csv.DictReader(open(man), delimiter="\t"):
    c = r["cohort"]
    if c not in cache:
        cache[c] = {x["sample_id"]: x["alignment"] for x in csv.DictReader(open(f"{sheets}/{c}.tsv"), delimiter="\t")}
    aln = cache[c][r["sample_id"]]
    rows.append((f"{c}/{r['sample_id']}", aln, os.stat(aln).st_size))
rows.sort(key=lambda x: -x[2])
with open(out, "w") as fh:
    for k, a, b in rows:
        fh.write(f"{k}\t{a}\t{b}\n")
print(f"{len(rows)} libraries, {sum(r[2] for r in rows) / 1e12:.2f} TB -> {out}", file=sys.stderr)
PY
    ;;
reuse)
    P1=$3
    want=$(awk '{n+=int(($2+999)/1000)} END{print n}' "$FAI")
    n=0; m=0
    while IFS=$'\t' read -r key aln bytes; do
        src="$P1/$key"; dst="$ROOT/$key"
        [ -f "$src/depth.ok" ] || continue
        [ -f "$dst/depth.ok" ] && continue
        a=$(awk -F'\t' 'NR==2{print $2}' "$src/timing.tsv"); b=$(awk -F'\t' 'NR==2{print $3}' "$src/timing.tsv")
        if [ "$a" != "$aln" ] || [ "$b" != "$bytes" ]; then echo "not reused (alignment differs): $key" >&2; m=$((m+1)); continue; fi
        got=$(zcat "$src/depth.regions.bed.gz" | wc -l)
        [ "$got" = "$want" ] || { echo "not reused (bins $got): $key" >&2; m=$((m+1)); continue; }
        mkdir -p "$dst"
        for f in depth.regions.bed.gz depth.regions.bed.gz.csi depth.mosdepth.summary.txt timing.tsv; do
            ln "$src/$f" "$dst/$f" 2>/dev/null || cp "$src/$f" "$dst/$f"
        done
        cmp -s "$src/depth.regions.bed.gz" "$dst/depth.regions.bed.gz" || { echo "copy differs: $key" >&2; exit 1; }
        sha256sum "$dst/depth.regions.bed.gz" | awk '{print $1}' > "$dst/depth.sha256"
        echo "$src" > "$dst/reused_from"
        touch "$dst/depth.ok"; n=$((n+1))
    done < "$ROOT/tasks.tsv"
    echo "reused $n, not reused $m" >&2
    ;;
worker)
    PAR=${3:-1}; export THREADS=${4:-4}
    bash "$HERE/depth_worker.sh" <(cut -f1,2 "$ROOT/tasks.tsv") "$ROOT" "$PAR"
    ;;
ledger)
    python3 - "$ROOT" "$FAI" <<'PY'
import csv, os, sys
root, fai = sys.argv[1:]
want = sum(-(-int(l.split()[1]) // 1000) for l in open(fai))
out = open(f"{root}/ledger.tsv", "w")
out.write("key\tstatus\thost\tbytes\tcram_version\ttotal_s\tMB_per_s\tbins\tsha256\tsource\n")
n = {"done": 0, "running": 0, "failed": 0, "todo": 0}
tb = ts = 0
for l in open(f"{root}/tasks.tsv"):
    key, aln, bytes_ = l.rstrip("\n").split("\t")
    d = f"{root}/{key}"
    if os.path.exists(f"{d}/depth.ok"):
        t = next(csv.DictReader(open(f"{d}/timing.tsv"), delimiter="\t"))
        sha = open(f"{d}/depth.sha256").read().strip() if os.path.exists(f"{d}/depth.sha256") else "NA"
        src = open(f"{d}/reused_from").read().strip() if os.path.exists(f"{d}/reused_from") else "this_pass"
        st = "done" if int(t["bins"]) == want and sha != "NA" else "CHECK"
        out.write("\t".join([key, st, t["host"], t["bytes"], t["cram_version"], t["total_s"], t["MB_per_s"], t["bins"], sha, src]) + "\n")
        n["done"] += st == "done"
        if src == "this_pass":
            tb += int(t["bytes"]); ts += int(t["total_s"])
    elif os.path.exists(f"{root}/failed/{key.replace('/', '__')}"):
        out.write(f"{key}\tfailed\t\t{bytes_}\t\t\t\t\t\t\n"); n["failed"] += 1
    elif os.path.isdir(f"{root}/claims/{key}"):
        h = open(f"{root}/claims/{key}/host").read().split()[0] if os.path.exists(f"{root}/claims/{key}/host") else ""
        out.write(f"{key}\trunning\t{h}\t{bytes_}\t\t\t\t\t\t\n"); n["running"] += 1
    else:
        out.write(f"{key}\ttodo\t\t{bytes_}\t\t\t\t\t\t\n"); n["todo"] += 1
out.close()
print(n, f"this pass: {tb / 1e12:.2f} TB in {ts / 3600:.1f} stream-h", file=sys.stderr)
PY
    ;;
*) echo "mode: tasks|reuse|worker|ledger" >&2; exit 2 ;;
esac
