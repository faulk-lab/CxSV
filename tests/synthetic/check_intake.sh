#!/usr/bin/env bash
# =============================================================================
# Patient intake (workflow/intake/Snakefile, roadmap item 16k) on synthetic inputs.
#
#   bash tests/synthetic/check_intake.sh [WORKDIR]      (default: .intake_test)
#
# Four samples, one per kind of input the intake must accept:
#   S_modbam  unaligned dorado-style modBAM (MM/ML/MN/qs + a per-read RG:Z), with
#             unmappable and chimeric reads
#   S_fastq   gzipped FASTQ whose comments mix SAM tags and non-SAM tokens
#   S_multi   two BAMs (two flow cells) aligned by "another aligner" to the same
#             sequence under other contig names (A/B), with supplementary records
#   S_cram    a CRAM against that other reference, decoded with input_reference
# check_intake_results.py then requires, per sample: verdict PASS; every input read
# exactly once among the primary/unmapped records; MM/ML identical per read; one
# read group with SM = sample_id and no dorado RG copied; hs1-style contigs; the
# patient sheet loads with workflow/lib/cxsv_inputs.py.
# Negative tests: a k=15 index is refused; a CRAM from another reference without
# input_reference fails (no CRAM left behind); a PacBio row and a duplicate sample
# make the sheet fail at DAG build.
# Needs the cxsv-pipeline env (minimap2 2.24, samtools, pysam, snakemake).
# =============================================================================
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
W="$(mkdir -p "${1:-$REPO/.intake_test}" && cd "${1:-$REPO/.intake_test}" && pwd)"
CORES="${CORES:-4}"
cd "$REPO"
D="$W/data"
rm -rf "$W/results" "$W/neg_"*
python tests/synthetic/make_intake_data.py "$D"

# alignments "from another aligner / reference naming": other contig names, no --eqx, default k
for s in other1 other2; do
  samtools fastq -T MM,ML,MN,qs "$D/reads/$s.mod.bam" 2>/dev/null | minimap2 -t 2 -y -ax map-ont "$D/ref/other.fa" - 2>/dev/null \
    | samtools sort -o "$D/reads/$s.other.bam" - && samtools index "$D/reads/$s.other.bam"
done
samtools fastq -T MM,ML,MN,qs "$D/reads/othercram.mod.bam" 2>/dev/null | minimap2 -t 2 -y -ax map-ont "$D/ref/other.fa" - 2>/dev/null \
  | samtools sort -O cram --reference "$D/ref/other.fa" -o "$D/reads/othercram.other.cram" -

printf 'sample_id\tinput\tindividual_id\tplatform\tinput_reference\tdescription\n' > "$W/intake_sheet.tsv"
printf 'S_modbam\t%s\tS_modbam\tONT\t\tdorado 0.8.3 synthetic sup\n' "$D/reads/modbam.mod.bam" >> "$W/intake_sheet.tsv"
printf 'S_fastq\t%s\t\tONT\t\t\n' "$D/reads/fastq.fastq.gz" >> "$W/intake_sheet.tsv"
printf 'S_multi\t%s,%s\tindiv_multi\tONT\t\t\n' "$D/reads/other1.other.bam" "$D/reads/other2.other.bam" >> "$W/intake_sheet.tsv"
printf 'S_cram\t%s\tS_cram\tONT\t%s\t\n' "$D/reads/othercram.other.cram" "$D/ref/other.fa" >> "$W/intake_sheet.tsv"

cat > "$W/site.yaml" <<YAML
reference: "$D/ref/ref.fa"
intake:
  sheet: "$W/intake_sheet.tsv"
  results_root: "$W/results/intake"
  minimap2_index: null
  index_mem_mb: 1000
  minibatch: "20M"
  align: {threads: 2, sort_threads: 1, sort_mem: "100M", mem_mb: 2000}
YAML
snakemake -s workflow/intake/Snakefile --cores "$CORES" --configfile "$W/site.yaml" > "$W/intake.log" 2>&1 \
  || { tail -50 "$W/intake.log"; for l in "$W"/results/intake/*/logs/intake_align.log; do echo "== $l"; tail -20 "$l"; done; exit 1; }
python tests/synthetic/check_intake_results.py "$W"

fail=0
# negative 1: index with the wrong k is refused before any alignment
minimap2 -x map-ont -k 15 -d "$W/k15.mmi" "$D/ref/ref.fa" 2>/dev/null
if REF_FA="$D/ref/ref.fa" REF_MMI="$W/k15.mmi" bash workflow/scripts/intake_align.sh NEG "$W/neg_k15" "$D/reads/modbam.mod.bam" > "$W/neg_k15.log" 2>&1; then
  echo "FAIL negative: k=15 index accepted"; fail=1
elif ! grep -q 'recipe needs k=17' "$W/neg_k15.log"; then echo "FAIL negative: k=15 refused for the wrong reason"; cat "$W/neg_k15.log"; fail=1
else echo "PASS negative: k=15 index refused"; fi
# negative 2: CRAM from another reference without input_reference: empty MD5 cache (no
# download) and the FASTA named in the CRAM's @SQ UR: moved away, as at a site without it
mkdir -p "$W/empty_cache"; mv "$D/ref/other.fa" "$D/ref/other.fa.away"
if REF_PATH="$W/empty_cache/%2s/%2s/%s" REF_CACHE="$W/empty_cache/%2s/%2s/%s" REF_FA="$D/ref/ref.fa" \
   REF_MMI="$W/results/intake/ref/hs1.map-ont.k17.mmi" bash workflow/scripts/intake_align.sh NEG "$W/neg_cram" \
   "$D/reads/othercram.other.cram" > "$W/neg_cram.log" 2>&1; then
  echo "FAIL negative: undecodable CRAM accepted"; fail=1
elif [ -e "$W/neg_cram/NEG.hs1.cram" ]; then echo "FAIL negative: CRAM left behind after failure"; fail=1
elif ! grep -q 'needs input_reference' "$W/neg_cram.log"; then echo "FAIL negative: undecodable CRAM refused without the input_reference hint"; cat "$W/neg_cram.log"; fail=1
else echo "PASS negative: CRAM without its reference refused, nothing left behind"; fi
mv "$D/ref/other.fa.away" "$D/ref/other.fa"
# negative 3: PacBio row and duplicate sample fail the sheet at DAG build
{ head -2 "$W/intake_sheet.tsv"; printf 'S_pb\t%s\tS_pb\tPACBIO\t\t\n' "$D/reads/modbam.mod.bam";
  sed -n 2p "$W/intake_sheet.tsv"; } > "$W/bad_sheet.tsv"
sed "s|$W/intake_sheet.tsv|$W/bad_sheet.tsv|" "$W/site.yaml" > "$W/bad_site.yaml"
if snakemake -s workflow/intake/Snakefile -n --cores 1 --configfile "$W/bad_site.yaml" > "$W/bad_sheet.log" 2>&1; then
  echo "FAIL negative: bad sheet accepted"; fail=1
elif grep -q "platform 'PACBIO' not supported" "$W/bad_sheet.log" && grep -q "duplicate sample_id 'S_modbam'" "$W/bad_sheet.log"; then
  echo "PASS negative: PacBio row and duplicate sample refused"
else echo "FAIL negative: bad sheet refused without both reasons"; cat "$W/bad_sheet.log"; fail=1; fi
[ $fail -eq 0 ] && echo "INTAKE TEST PASSED" || { echo "INTAKE TEST FAILED"; exit 1; }
