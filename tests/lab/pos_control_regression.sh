#!/usr/bin/env bash
# =============================================================================
# LAB-ONLY regression test with the real positive control (never in CI).
#
# The positive control is sample 3559 (human-subject data; lab NAS only) with
# a catalogue of 12 heterozygous planted CxSVs (tests/positive_control/events.yaml,
# docs/POSITIVE_CONTROL.md). Nothing derived from it beyond counts, coordinates
# and the planted-event truth set may leave the NAS.
#
# What it does (lab NAS / node-local scratch):
#   1. slices every planted event's replace window +/- 1 Mb (merged) from the
#      positive-control CRAM and from N background CRAMs (registry sheet, >=17x);
#   2. builds a panel from the background slices and runs the patient scorer on
#      the positive-control slice;
#   3. scores each event against the truth set (score_pos_control.py):
#      junction recovery, CxSV classification, private/public status;
#   4. with RERUN=1, runs the patient workflow a second time into a fresh
#      results root and requires byte-identical private tables (determinism).
#
#   bash tests/lab/pos_control_regression.sh [WORKDIR] [N_BACKGROUND]
# Env: CORES (16), RERUN (0), POS_CRAM / TRUTH to override the release paths.
# Needs the cxsv-pipeline env. Exit 0 = every event meets its expectation or is a
# listed known failure (tests/positive_control/known_failures.tsv) that still fails.
# =============================================================================
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
W="$(mkdir -p "${1:-$HOME/scratch/cxsv_posctrl_regression}" && cd "${1:-$HOME/scratch/cxsv_posctrl_regression}" && pwd)"
N="${2:-10}"
CORES="${CORES:-16}"
RERUN="${RERUN:-0}"
PC=/mnt/sargasso1/CxSV/datasets/pos_control
POS_CRAM="${POS_CRAM:-$PC/release/pos_control_3559.hs1.minimap2.cram}"
TRUTH="${TRUTH:-$PC/release/truth}"
REF=/mnt/sargasso1/CxSV/refs/hs1.fa
CACHE=/mnt/sargasso1/CxSV/refs/ref_cache
RES=/mnt/sargasso1/CxSV/resources
SHEET=/mnt/sargasso1/CxSV/datasets/registry/sample_sheets/1kGP_ONT_1019_CHM13.tsv
EVENTS="$REPO/tests/positive_control/events.yaml"
PAD=1000000
mkdir -p "$W/crams"
cd "$REPO"

# 1. slice regions: every replace window +/- 1 Mb, merged
awk -v pad=$PAD 'BEGIN{OFS="\t"} {s=$2-pad; if (s<0) s=0; print $1, s, $3+pad}' "$TRUTH/replace_all.bed" \
  | sort -k1,1 -k2,2n | bedtools merge -i - > "$W/regions.bed"
echo "slice: $(wc -l < "$W/regions.bed") regions, $(awk '{s+=$3-$2} END{printf "%.1f Mb", s/1e6}' "$W/regions.bed")"

slice() {  # slice CRAM_IN CRAM_OUT
    [ -f "$2.crai" ] && return 0
    samtools view -C -M -L "$W/regions.bed" -T "$REF" -o "$2" "$1" && samtools index "$2"
}
slice "$POS_CRAM" "$W/crams/pos_control.cram"

printf 'sample_id\talignment\tindividual_id\tplatform\n' > "$W/background.tsv"
# awk stops after N rows itself (a `| head` here would SIGPIPE awk under pipefail)
awk -F'\t' -v n="$N" 'NR==1{for(i=1;i<=NF;i++) c[$i]=i; next}
    $c["passes_17x"]=="yes"{print $c["sample_id"]"\t"$c["alignment"]; if (++k >= n) exit}' "$SHEET" \
  | while read -r sid cram; do
    slice "$cram" "$W/crams/bg_$sid.cram"
    printf '%s\t%s\t%s\tONT\n' "$sid" "$W/crams/bg_$sid.cram" "$sid" >> "$W/background.tsv"
done
printf 'sample_id\talignment\tindividual_id\tplatform\npos_control\t%s\tpos_control_3559\tONT\n' \
  "$W/crams/pos_control.cram" > "$W/patient.tsv"

site_yaml() {  # site_yaml PATIENT_ROOT > file
cat <<YAML
reference: "$REF"
ref_cache: "$CACHE"
annotation:
  genes_bed:   "$RES/genes.final.bed"
  exons_bed:   "$RES/exons.filtered.bed"
  repeats_bed: "$RES/repeats.final.bed"
  segdup_bed:  "$RES/segdup.final.bed"
threads: {samtools: 4, sv_callers: $CORES}
panel_results_root:   "$W/results/panel"
patient_results_root: "$1"
panel_release:        "$W/results/panel/release"
panel_cohorts:
  lab_background_pos_control_slices: {samples: "$W/background.tsv"}
patient_runs:
  pos_control: {samples: "$W/patient.tsv"}
YAML
}
site_yaml "$W/results/patient" > "$W/site.yaml"

# 2. both workflows
SM=(snakemake --cores "$CORES" --keep-going)
"${SM[@]}" -s workflow/panel/Snakefile   --configfile "$W/site.yaml" > "$W/panel.log" 2>&1   || { tail -30 "$W/panel.log"; exit 1; }
"${SM[@]}" -s workflow/patient/Snakefile --configfile "$W/site.yaml" > "$W/patient.log" 2>&1 || { tail -30 "$W/patient.log"; exit 1; }

# 3. score against the truth set
S="$W/results/patient/samples/pos_control/pos_control"
set +e
python tests/positive_control/score_pos_control.py "$TRUTH" "$S" "$EVENTS" "$W/score" \
    --known-failures tests/positive_control/known_failures.tsv | tee "$W/score.txt"
set -e
# FAIL = new failure; XPASS = a listed known failure now passes (update the list)
FAILS=$(awk -F'\t' 'NR>1 && ($NF=="FAIL" || $NF=="XPASS")' "$W/score.events.tsv" | wc -l)

# 4. optional determinism check
if [ "$RERUN" = "1" ]; then
    site_yaml "$W/results/patient_rerun" > "$W/site_rerun.yaml"
    "${SM[@]}" -s workflow/patient/Snakefile --configfile "$W/site_rerun.yaml" > "$W/patient_rerun.log" 2>&1 \
      || { tail -30 "$W/patient_rerun.log"; exit 1; }
    S2="$W/results/patient_rerun/samples/pos_control/pos_control"
    for f in merged_sv.vcf private/cxsv_panel_status.tsv private/private_cxsv.vcf; do
        # VCF '##' meta lines carry dates and command lines; compare records only
        if cmp -s <(grep -v '^##' "$S/$f") <(grep -v '^##' "$S2/$f"); then echo "rerun identical: $f"
        else echo "RERUN DIFFERS: $f"; FAILS=$((FAILS+1)); fi
    done
fi
echo "new failures or unexpected passes: $FAILS"
[ "$FAILS" -eq 0 ]
