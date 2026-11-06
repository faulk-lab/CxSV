#!/usr/bin/env bash
# =============================================================================
# run_phase1.sh — roadmap item 15 Phase 1 chain after the depth pass (all in the
# pipeline image; freya-sized: ~30 min for 100 genomes + LOFO):
#   1 make_genomes.py   genomes.tsv (100 subset genomes + PC + HG01868), sex checked on depth
#   2 segment_all.py    PoN per stratum (unflagged), LOFO segmentation, segments_all.tsv
#   3 link_all.py       candidates + augmented discovery trees (OUT/aug/samples)
#   4 subset manifest   the panel-v1.1 sample_manifest rows of the subset
#   5 cn_lofo.sh        test D (healthy Tier 1 per genome, LOFO, deployed scorer)
#
#   run_phase1.sh OUT [STAGES=12345]
# =============================================================================
set -euo pipefail
OUT=$1; STAGES=${2:-12345}
HERE="$(cd "$(dirname "$0")" && pwd)"; CODE=$(cd "$HERE/../.." && pwd)
B=/mnt/sargasso1/CxSV/calibration/cn_linkage
SIF=${SIF:-/mnt/sargasso1/CxSV/bundles/cxsv-bundle-20261105-8fc5bef725bd/container/cxsv-pipeline.sif}
X=${X:-singularity exec -B /mnt/sargasso1 -B /home/minknow/scratch $SIF}
R2=/mnt/sargasso1/CxSV/build/panel-v1.1-r2/results/panel
FAM=/mnt/sargasso1/CxSV/calibration/item10/families.tsv
PARBED=/mnt/sargasso1/CxSV/resources/scoring/par.bed
PC_SRC=${PC_SRC:-/home/minknow/scratch/item11r2_full/run/results/patient/samples/pos_control/pos_control}
HG_SRC=${HG_SRC:-/home/minknow/scratch/item_tag/bundle_validation/normal/results/patient/samples/standin_HG01868/HG01868}
mkdir -p "$OUT/logs"
t() { echo "$(date +%T) $*" | tee -a "$OUT/logs/run_phase1.log"; }
if [[ $STAGES == *1* ]]; then t "1 make_genomes"
  $X python3 "$HERE/make_genomes.py" --subset "$HERE/subset100.tsv" --depth-root "$B/depth" --panel-samples "$R2/samples" \
     --extra "extra/pos_control_L,LRSC_R10,$PC_SRC" --extra "extra/HG01868,LRSC_R9,$HG_SRC" \
     --out "$OUT/genomes.tsv" 2> "$OUT/logs/make_genomes.log"; fi
if [[ $STAGES == *2* ]]; then t "2 segment_all"
  $X python3 "$HERE/segment_all.py" --genomes "$OUT/genomes.tsv" --depth-root "$B/depth" --families "$FAM" \
     --par "$PARBED" --out "$OUT/seg" --threads 8 2> "$OUT/logs/segment_all.log"; fi
if [[ $STAGES == *3* ]]; then t "3 link_all"
  $X python3 "$HERE/link_all.py" --genomes "$OUT/genomes.tsv" --seg "$OUT/seg" --code "$CODE" \
     --fai /mnt/sargasso1/CxSV/refs/hs1.fa.fai --families "$FAM" --out "$OUT/aug" --threads 8 2> "$OUT/logs/link_all.log"; fi
if [[ $STAGES == *4* ]]; then t "4 subset manifest"
  python3 - "$R2/panel/sample_manifest.tsv" "$HERE/subset100.tsv" "$OUT/aug/subset_manifest.tsv" <<'PY'
import csv, sys
man, sub, out = sys.argv[1:]
keep = {(r["cohort"], r["sample_id"]) for r in csv.DictReader(open(sub), delimiter="\t")}
with open(man) as fh, open(out, "w") as fo:
    head = fh.readline(); fo.write(head); n = 0
    for l in fh:
        f = l.rstrip("\n").split("\t")
        if (f[0], f[1]) in keep:
            fo.write(l); n += 1
print(f"subset manifest rows {n} of {len(keep)}", file=sys.stderr)
PY
fi
if [[ $STAGES == *5* ]]; then t "5 cn_lofo"
  CODE=$CODE bash "$HERE/cn_lofo.sh" "$OUT/aug" "$OUT/lofo" 8 > "$OUT/logs/cn_lofo.log" 2>&1; fi
t "done"
