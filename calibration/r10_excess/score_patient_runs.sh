#!/usr/bin/env bash
# Score the candidate-rule patient reruns (patient_rerun_downstream.sh): positive control
# against its truth (score_pos_control.py, E01-E12 + E13 add-on) and host / stand-in
# Tier counts (loci not overlapping a planted event).
#   score_patient_runs.sh BUNDLE CODE OUTDIR RUN...   (RUN = pc_<rule> or hg01868_<rule>)
set -euo pipefail
BUNDLE=$1; CODE=$2; O=$3; shift 3
SIF=$BUNDLE/container/cxsv-pipeline.sif; PC=/mnt/sargasso1/CxSV/datasets/pos_control
for r in "$@"; do
  W=$O/$r
  case $r in
    pc_*) S=$W/run/results/patient/samples/pos_control/pos_control
          singularity exec --containall -B "$W:$W" -B "$CODE:/code:ro" -B "$PC:$PC:ro" "$SIF" \
            python /code/tests/positive_control/score_pos_control.py $PC/release/truth $S \
            /code/tests/positive_control/events.yaml $W/score --extra-truth $PC/addons/E13/truth \
            --known-failures /code/tests/positive_control/known_failures.tsv > $W/score.txt 2>&1 || true ;;
    *)    S=$(ls -d $W/run/results/patient/samples/*/*) ;;
  esac
  python3 - "$S/scoring/cxsv_scores.tsv" "$r" <<'PY'
import csv, sys, collections
truth = [("chrX",152320014,152349014),("chr17",32243849,32271849),("chr17",15140559,15166559),("chr22",51207890,51230915),
         ("chr2",50453329,50488904),("chr7",71340258,71406158),("chr2",199854293,199856293),("chr7",115868955,115869905),
         ("chr4",1931409,1931784),("chr9",149997337,149997837),("chr13",99195184,99195684),("chr15",90515727,90925728),
         ("chr5",25003550,25041800),("chr10",29098000,29115650),("chr12",74754175,74814325),("chr13",54982500,55017500)]
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
pc = sys.argv[2].startswith("pc_")
host = [x for x in rows if not (pc and any(x["Chr"] == c and int(x["Start"]) <= e and int(x["End"]) >= b for c, b, e in truth))]
t = collections.Counter(x["Tier"] for x in host)
print(f"{sys.argv[2]}\tloci {len(rows)}\thost/stand-in loci {len(host)}\tTier1 {t['Tier1_high_priority']}\tTier2 {t['Tier2']}\tTier3 {t['Tier3']}")
PY
  [ -f $W/score.txt ] && grep -E '^E[0-9]+\s' $W/score.txt | cut -c1-230 || true
done
