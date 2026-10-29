#!/usr/bin/env bash
# =============================================================================
# Lab-fleet regression test with the real positive control (too large for CI).
#
# The positive control is sample 3559 (exempt; its ~105 GB CRAM is kept on the
# lab NAS for size) with a catalogue of 13 heterozygous planted events
# (tests/positive_control/events.yaml, docs/POSITIVE_CONTROL.md): E01-E12 are
# in the release CRAM; E13, a COMMON panel CxSV, is an add-on spliced here into
# the control AND into 3 background genomes from fixed, checksummed read sets
# (datasets/pos_control/addons/E13, made by make_spikein_reads.py).
#
# What it does (lab NAS / node-local scratch):
#   1. slices every planted event's replace window +/- 1 Mb (merged; release +
#      add-on windows) from the positive-control CRAM and from N background
#      CRAMs (registry sheet, first N with passes_17x);
#   2. splices E13 into the control slice and into its background carriers
#      (events.yaml E13.spike_in.backgrounds), after verifying the add-on MD5SUMS;
#   3. builds a panel from the background slices (junction sites + Sniffles2
#      population genotypes) and runs the patient scorer on the control slice;
#   4. scores each event against the truth set (score_pos_control.py):
#      junction recovery, CxSV classification, breakpoint-level panel status
#      (private / shared_junctions / public, E13 public with 3 carriers), and the
#      pathogenicity score (E01-E09 top tier; E10, E13 below them; E11/E12 not
#      scored; config/scoring.yaml, external catalogues from the NAS);
#   5. with RERUN=1, runs the patient workflow a second time into a fresh
#      results root and requires byte-identical records (determinism).
#
#   bash tests/lab/pos_control_regression.sh [WORKDIR] [N_BACKGROUND]
# Env: CORES (16), RERUN (0), POS_CRAM / TRUTH / ADDON to override the release paths.
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
ADDON="${ADDON:-$PC/addons/E13}"
REF=/mnt/sargasso1/CxSV/refs/hs1.fa
CACHE=/mnt/sargasso1/CxSV/refs/ref_cache
RES=/mnt/sargasso1/CxSV/resources
SHEET=/mnt/sargasso1/CxSV/datasets/registry/sample_sheets/1kGP_ONT_1019_CHM13.tsv
EVENTS="$REPO/tests/positive_control/events.yaml"
PAD=1000000
INDIVIDUALS=/mnt/sargasso1/CxSV/datasets/registry/individuals.tsv
mkdir -p "$W/crams/raw"
cd "$REPO"

# 0. the E13 add-on must be intact
(cd "$ADDON" && md5sum -c --quiet MD5SUMS) || { echo "E13 add-on fails its MD5SUMS: $ADDON"; exit 1; }
SPIKE_BED="$ADDON/reads/E13.replace.bed"
SPIKE_HOSTS=$(python -c "import yaml,sys; e=[x for x in yaml.safe_load(open('$EVENTS'))['events'] if x['id']=='E13'][0]; print(' '.join(e['spike_in']['backgrounds']))")

# 1. slice regions: every replace window (release + E13 add-on) +/- 1 Mb, merged
cat "$TRUTH/replace_all.bed" "$SPIKE_BED" \
  | awk -v pad=$PAD 'BEGIN{OFS="\t"} {s=$2-pad; if (s<0) s=0; print $1, s, $3+pad}' \
  | sort -k1,1 -k2,2n | bedtools merge -i - > "$W/regions.bed"
echo "slice: $(wc -l < "$W/regions.bed") regions, $(awk '{s+=$3-$2} END{printf "%.1f Mb", s/1e6}' "$W/regions.bed")"

slice() {  # slice CRAM_IN CRAM_OUT
    [ -f "$2.crai" ] && return 0
    samtools view -C -M -L "$W/regions.bed" -T "$REF" -o "$2" "$1" && samtools index "$2"
}
splice_spikein() {  # splice_spikein SLICE SIM_BAM CRAM_OUT: drop host reads in the E13 window, add simulated
    [ -f "$3.crai" ] && return 0
    samtools view --reference "$REF" -L "$SPIKE_BED" "$1" | cut -f1 | sort -u > "$3.dropped.txt"
    samtools view -u --reference "$REF" -N "^$3.dropped.txt" "$1" \
      | samtools merge -@ 4 --reference "$REF" -O cram -o "$3" - "$2"
    samtools index "$3"
    echo "spliced E13 into $(basename "$3"): dropped $(wc -l < "$3.dropped.txt") host reads," \
         "added $(samtools view -c "$2") simulated records"
}
slice "$POS_CRAM" "$W/crams/raw/pos_control.cram"
splice_spikein "$W/crams/raw/pos_control.cram" "$ADDON/reads/pos_control.E13.sim.bam" "$W/crams/pos_control.cram"

printf 'sample_id\talignment\tindividual_id\tplatform\n' > "$W/background.tsv"
# awk stops after N rows itself (a `| head` here would SIGPIPE awk under pipefail)
awk -F'\t' -v n="$N" 'NR==1{for(i=1;i<=NF;i++) c[$i]=i; next}
    $c["passes_17x"]=="yes"{print $c["sample_id"]"\t"$c["alignment"]; if (++k >= n) exit}' "$SHEET" \
  | while read -r sid cram; do
    if [[ " $SPIKE_HOSTS " == *" $sid "* ]]; then
        slice "$cram" "$W/crams/raw/bg_$sid.cram"
        splice_spikein "$W/crams/raw/bg_$sid.cram" "$ADDON/reads/$sid.E13.sim.bam" "$W/crams/bg_$sid.cram"
    else
        slice "$cram" "$W/crams/bg_$sid.cram"
    fi
    printf '%s\t%s\t%s\tONT\n' "$sid" "$W/crams/bg_$sid.cram" "$sid" >> "$W/background.tsv"
done
for h in $SPIKE_HOSTS; do
    grep -q -P "^$h\t" "$W/background.tsv" || { echo "E13 carrier $h is not among the $N backgrounds"; exit 1; }
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
panel_individuals:    "$INDIVIDUALS"
release:
  name: "${RELEASE_NAME:-TEST-posctrl-slices}"
  date: ${RELEASE_DATE:-null}
  test_release: true
  description: "TEST RELEASE: $N background genomes sliced to the positive-control regions; not for patient interpretation"
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
    --extra-truth "$ADDON/truth" \
    --known-failures tests/positive_control/known_failures.tsv | tee "$W/score.txt"
cp "$W/results/panel/release/panel_frequency_comparison_summary.tsv" "$W/frequency_comparison_summary.tsv"
cp "$S/scoring/cxsv_scores.tsv" "$W/score.cxsv_scores.tsv"; cp "$S/scoring/cxsv_score_report.md" "$W/score.report.md"
cp "$S/scoring/cxsv_score_features.tsv" "$W/score.features.tsv"
set -e
# FAIL = new failure; XPASS = a listed known failure now passes (update the list)
FAILS=$(awk -F'\t' 'NR>1 && ($NF=="FAIL" || $NF=="XPASS")' "$W/score.events.tsv" | wc -l)

# 4. optional determinism check
if [ "$RERUN" = "1" ]; then
    site_yaml "$W/results/patient_rerun" > "$W/site_rerun.yaml"
    "${SM[@]}" -s workflow/patient/Snakefile --configfile "$W/site_rerun.yaml" > "$W/patient_rerun.log" 2>&1 \
      || { tail -30 "$W/patient_rerun.log"; exit 1; }
    S2="$W/results/patient_rerun/samples/pos_control/pos_control"
    for f in merged_sv.vcf cxsv/cxsv_summary.tsv junctions/junctions.tsv private/cxsv_panel_status.tsv \
             private/junction_panel_status.tsv private/private_cxsv.vcf private/cxsv_catalogue_hits.tsv \
             private/junction_catalogue_hits.tsv scoring/locus_depth.tsv scoring/cxsv_scores.tsv \
             scoring/cxsv_score_features.tsv scoring/cxsv_score_report.md; do
        # VCF '##' meta lines carry dates and command lines; compare records only
        if cmp -s <(grep -v '^##' "$S/$f") <(grep -v '^##' "$S2/$f"); then echo "rerun identical: $f"
        else echo "RERUN DIFFERS: $f"; FAILS=$((FAILS+1)); fi
    done
fi
echo "new failures or unexpected passes: $FAILS"
[ "$FAILS" -eq 0 ]
