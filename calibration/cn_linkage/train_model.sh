#!/usr/bin/env bash
# train_model.sh HOST_ALN PRESET OUTDIR [REGION] — badread error / qscore models and read
# statistics of one host (test B), as the positive-control model was built:
# 20,000 primary reads > 1 kb from REGION (default chr1:20,000,000-60,000,000), aligned
# with minimap2 -c (PAF with cigar), badread error_model / qscore_model, read length mean/sd
# and identity (matches / alignment length) mean, sd, p99.5; identity_max_sim = min(p99.5, 99.5).
# Env: badread, minimap2, samtools (cxsv-posctrl on freya).
set -euo pipefail
ALN=$1; PRESET=$2; OUT=$3; REGION=${4:-chr1:20000000-60000000}
REF=${REF:-/mnt/sargasso1/CxSV/refs/hs1.fa}
mkdir -p "$OUT"
samtools view -F 0x904 -m 1000 --reference "$REF" "$ALN" "$REGION" | awk 'NR%5==0' | head -n 20000 \
    | awk -F'\t' '{print "@"$1"\n"$10"\n+\n"$11}' > "$OUT/reads20k.fq" || true
minimap2 -t 8 -c -x "$PRESET" "$REF" "$OUT/reads20k.fq" > "$OUT/reads20k.paf" 2> "$OUT/mm.log"
badread error_model --reference "$REF" --reads "$OUT/reads20k.fq" --alignment "$OUT/reads20k.paf" > "$OUT/host_error_model" 2> "$OUT/em.log"
badread qscore_model --reference "$REF" --reads "$OUT/reads20k.fq" --alignment "$OUT/reads20k.paf" > "$OUT/host_qscore_model" 2> "$OUT/qm.log"
python3 - "$OUT" <<'PY'
import statistics, sys
out = sys.argv[1]
L = [len(l.strip()) for i, l in enumerate(open(f"{out}/reads20k.fq")) if i % 4 == 1]
best = {}
for l in open(f"{out}/reads20k.paf"):
    f = l.split("\t")
    idt = 100.0 * int(f[9]) / int(f[10])
    if f[0] not in best or int(f[10]) > best[f[0]][1]:
        best[f[0]] = (idt, int(f[10]))
I = sorted(v[0] for v in best.values())
p995 = I[int(0.995 * (len(I) - 1))]
rows = dict(n_reads=len(L), mean_len=round(statistics.mean(L)), sd_len=round(statistics.pstdev(L)),
            identity_mean=round(statistics.mean(I), 3), identity_sd=round(statistics.pstdev(I), 3),
            **{"identity_p99.5": round(p995, 3)}, identity_max_sim=round(min(p995, 99.5), 3))
with open(f"{out}/host_read_stats.tsv", "w") as fh:
    for k, v in rows.items():
        fh.write(f"{k}\t{v}\n")
print(rows)
PY
