#!/usr/bin/env bash
# =============================================================================
# nas_throughput.sh — sequential read throughput from the NAS, per node alone
# and with all nodes reading at once (the panel build reads ~80 TB of
# alignments, several passes each, so the NAS is a shared bottleneck).
#
#   bash scripts/panel_build/nas_throughput.sh OUT.tsv GIB FILE_LIST NODE...
#   STREAMS=2 SKIP_ALONE=1 bash ... (phase 2 only, 2 concurrent reads per node)
#
# Run on freya. FILE_LIST: one large alignment per line, at least 2 x the
# number of nodes, all distinct (no node re-reads pages another node pulled
# into the NAS cache). Phase 1 reads one file per node, one node at a time;
# phase 2 all nodes at once on fresh files. Each read is GIB GiB with
# `dd bs=4M`, timed on the reading node (a time-limited dd stopped by SIGINT
# does not print its summary on every node, so the amount is fixed instead).
# Output columns: node, phase, file, bytes, seconds, MB_per_s.
# =============================================================================
set -uo pipefail
OUT="$1"; GIB="$2"; LIST="$3"; shift 3
mapfile -t FILES < <(grep -v '^#' "$LIST" | grep -v '^$')
STREAMS="${STREAMS:-1}"; SKIP_ALONE="${SKIP_ALONE:-0}"
need=$(( (SKIP_ALONE == 1 ? 0 : $#) + STREAMS * $# ))
[ "${#FILES[@]}" -ge "$need" ] || { echo "need >= $need files in $LIST"; exit 1; }
BYTES=$((GIB * 1073741824))
probe() {   # node file -> "seconds"
    local cmd="t0=\$(date +%s.%N); dd if='$2' of=/dev/null bs=4M count=$((GIB * 256)) status=none; t1=\$(date +%s.%N); echo \"\$t1 \$t0\" | awk '{print \$1-\$2}'"
    if [ "$1" = "$(hostname)" ]; then bash -c "$cmd"; else ssh -o BatchMode=yes "$1" "$cmd"; fi
}
row() { echo -e "$1\t$2\t$3\t$BYTES\t$4\t$(awk -v b="$BYTES" -v s="$4" 'BEGIN{printf "%.0f", b/s/1e6}')" >> "$OUT"; }
[ -f "$OUT" ] || echo -e "node\tphase\tfile\tbytes\tseconds\tMB_per_s" > "$OUT"
i=0
if [ "$SKIP_ALONE" != 1 ]; then
    for n in "$@"; do                               # phase 1: alone
        row "$n" alone "${FILES[$i]}" "$(probe "$n" "${FILES[$i]}")"; i=$((i + 1))
    done
fi
tmp=$(mktemp -d); j=$i
for n in "$@"; do                                   # phase 2: concurrent, STREAMS per node
    for k in $(seq 1 "$STREAMS"); do
        ( probe "$n" "${FILES[$j]}" > "$tmp/$n.$k" ) & j=$((j + 1))
    done
done
wait
for n in "$@"; do for k in $(seq 1 "$STREAMS"); do
    row "$n" "concurrent_x$STREAMS" "${FILES[$i]}" "$(cat "$tmp/$n.$k")"; i=$((i + 1))
done; done
rm -rf "$tmp"
awk -F'\t' -v ph="concurrent_x$STREAMS" 'NR > 1 && $2 == ph {t += $6}
            END {printf "aggregate %s: %.0f MB/s (sum of per-read rates)\n", ph, t}' "$OUT"
