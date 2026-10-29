#!/usr/bin/env bash
# =============================================================================
# End-to-end test of both workflows on the fully synthetic dataset.
#
#   bash tests/synthetic/run_synthetic.sh [WORKDIR]     (default: .synthetic_test)
#
# 1. make_synthetic_data.py: seeded 1.2 Mb reference, 3 panel genomes + 1
#    patient with planted COMMON / RARE / PRIVATE CxSVs, ONT-like reads,
#    minimap2 alignments (panel3 as CRAM to exercise the CRAM path).
# 2. panel builder on panel1-3, patient scorer on the patient.
# 3. check_synthetic_results.py: carriers 3 / 1 / 0 and public / public /
#    private by breakpoint matching, although a HOTSPOT DECOY (a different
#    CxSV at the PRIVATE locus in panel1 and panel3) overlaps PRIVATE; the
#    SIMPLE deletion's junction has 3 panel carriers; compatible fingerprint;
#    private VCF confined to PRIVATE.
# 4. negative tests: a patient overlay that changes a fingerprinted parameter
#    must stop at check_panel_compatibility; a panel release with one byte
#    changed must be refused at verify_panel_release (manifest check).
# Needs the cxsv-pipeline env (it includes minimap2 for this test). Run from
# the repository root. Exit status 0 = all checks passed.
# =============================================================================
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
W="$(mkdir -p "${1:-$REPO/.synthetic_test}" && cd "${1:-$REPO/.synthetic_test}" && pwd)"
CORES="${CORES:-4}"
cd "$REPO"

if [ ! -f "$W/data/expected.json" ]; then
    python tests/synthetic/make_synthetic_data.py "$W/data" --threads "$CORES"
fi
D="$W/data"

printf 'sample_id\talignment\tindividual_id\tplatform\n' > "$W/panel_sheet.tsv"
for s in panel1 panel2; do printf '%s\t%s\t%s\tONT\n' $s "$D/bams/$s.bam" $s >> "$W/panel_sheet.tsv"; done
printf 'panel3\t%s\tpanel3\tONT\n' "$D/bams/panel3.cram" >> "$W/panel_sheet.tsv"
printf 'sample_id\talignment\tindividual_id\tplatform\npatient\t%s\tpatient\tONT\n' "$D/bams/patient.bam" > "$W/patient_sheet.tsv"

cat > "$W/site.yaml" <<YAML
reference: "$D/ref/chr20.fa"
ref_cache: "$W/ref_cache"
annotation:
  genes_bed:   "$D/resources/genes.bed"
  exons_bed:   "$D/resources/exons.bed"
  repeats_bed: "$D/resources/repeats.bed"
  segdup_bed:  "$D/resources/segdup.bed"
threads: {samtools: 2, sv_callers: $CORES}
panel_results_root:   "$W/results/panel"
patient_results_root: "$W/results/patient"
panel_release:        "$W/results/panel/release"
release: {name: "synthetic", date: "20260101", test_release: true, description: "CI synthetic panel"}
panel_cohorts:
  synthetic_panel: {samples: "$W/panel_sheet.tsv"}
patient_runs:
  synthetic_patient: {samples: "$W/patient_sheet.tsv"}
YAML
cat > "$W/mismatch.yaml" <<YAML
patient_results_root: "$W/results/patient_mismatch"
cxsv_params: {min_breakpoints: 4}
YAML

SM=(snakemake --cores "$CORES" --keep-going --printshellcmds)
echo "== panel";   "${SM[@]}" -s workflow/panel/Snakefile   --configfile "$W/site.yaml" > "$W/panel.log" 2>&1 \
    || { tail -40 "$W/panel.log"; echo "panel workflow FAILED (log: $W/panel.log)"; exit 1; }
echo "== patient"; "${SM[@]}" -s workflow/patient/Snakefile --configfile "$W/site.yaml" > "$W/patient.log" 2>&1 \
    || { tail -40 "$W/patient.log"; echo "patient workflow FAILED (log: $W/patient.log)"; exit 1; }

echo "== checks"
RC=0
python tests/synthetic/check_synthetic_results.py "$D/expected.json" "$W/results/patient" synthetic_patient patient || RC=1

echo "== negative test: changed fingerprinted parameter must be rejected"
if "${SM[@]}" -s workflow/patient/Snakefile --configfile "$W/site.yaml" "$W/mismatch.yaml" > "$W/mismatch.log" 2>&1; then
    echo "FAIL patient run with cxsv_params.min_breakpoints changed was NOT rejected"; RC=1
elif grep -q "params.cxsv_params.min_breakpoints" "$W/results/patient_mismatch/logs/panel_compatibility.log"; then
    echo "PASS mismatch rejected and reported: $(grep min_breakpoints "$W/results/patient_mismatch/logs/panel_compatibility.log" | xargs)"
else
    echo "FAIL run failed, but not at the compatibility check (log: $W/mismatch.log)"; RC=1
fi
echo "== negative test: a tampered panel release must be refused"
rm -rf "$W/release_tampered" "$W/results/patient_tampered"
cp -r "$W/results/panel/release" "$W/release_tampered"
# same-size change: one digit of the first carrier site's pos1
python - "$W/release_tampered/panel_junction_carriers.tsv" <<'PY'
import sys
p = sys.argv[1]; lines = open(p).read().split("\n"); f = lines[1].split("\t")
f[1] = f[1][:-1] + str((int(f[1][-1]) + 1) % 10); lines[1] = "\t".join(f)
open(p, "w").write("\n".join(lines))
PY
cat > "$W/tampered.yaml" <<YAML
patient_results_root: "$W/results/patient_tampered"
panel_release: "$W/release_tampered"
YAML
if "${SM[@]}" -s workflow/patient/Snakefile --configfile "$W/site.yaml" "$W/tampered.yaml" > "$W/tampered.log" 2>&1; then
    echo "FAIL patient run against a tampered release was NOT refused"; RC=1
elif grep -q "sha256 differs: panel_junction_carriers.tsv" "$W/tampered.log" && [ ! -d "$W/results/patient_tampered/samples" ]; then
    echo "PASS tampered release refused before any job ran: $(grep -m1 'sha256 differs' "$W/tampered.log" | xargs)"
else
    echo "FAIL run failed, but not at the release check (log: $W/tampered.log)"; RC=1
fi
exit $RC
