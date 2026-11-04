#!/usr/bin/env bash
# =============================================================================
# Positive control, WHOLE GENOME, scored against the panel release shipped in an
# offline bundle (validation of the panel-v1.0 bundle; docs/POSITIVE_CONTROL.md
# section "Full panel").
#
#   bash tests/lab/pos_control_full_genome.sh BUNDLE WORKDIR [CORES]
#
# Unlike tests/lab/pos_control_regression.sh (10 sliced backgrounds, conda env),
# nothing is sliced and no panel is built: the patient workflow runs from the
# bundle (run_cxsv.sh patient, Singularity/Apptainer, network off) on the full
# positive-control CRAM against the bundle's own panel release.
#
#   0. provenance   MD5 of the release CRAM (datasets/pos_control/release/MD5SUMS),
#                   E13 add-on MD5SUMS, every @SQ M5 present in the bundle's cache
#   1. CRAM version a CRAM 3.1 input is rewritten as CRAM 3.0 exactly as the panel
#                   build did (scripts/panel_build/lib.sh cxsv_cram30): the pinned
#                   mosdepth 0.3.10 cannot decode 3.1 (empty callability, exit 0)
#   2. E13 splice   E13 exists only as an add-on (datasets/pos_control/addons/E13):
#                   host reads overlapping E13.replace.bed are dropped (by name,
#                   all their records) and the simulated E13 reads are merged in,
#                   as splice_spikein() in the regression, over the whole genome.
#                   Output CRAM 3.0. Record counts are checked:
#                   out = in - dropped records + simulated records
#   3. patient run  run_cxsv.sh patient (bundle only, --containall, --network none)
#   4. scoring      score_pos_control.py (truth + E13 add-on truth), in the SIF
#   5. release      the bundle's release still verifies
# All samtools/python steps run inside the bundle's SIF with no network.
# SKIP_MD5=1 skips the 100 GB MD5 of step 0. Exit 0 = the run finished (event
# verdicts are reported, not enforced: E13 is expected to FAIL against a real
# panel, whose genomes do not carry the spike-in).
# =============================================================================
set -euo pipefail
BUNDLE="$(cd "$1" && pwd)"; W="$(mkdir -p "$2" && cd "$2" && pwd)"; CORES="${3:-12}"
PC="${PC:-/mnt/sargasso1/CxSV/datasets/pos_control}"
POS_CRAM="${POS_CRAM:-$PC/release/pos_control_3559.hs1.minimap2.cram}"
ADDON="${ADDON:-$PC/addons/E13}"
# shellcheck disable=SC1091
source "$BUNDLE/bundle.env"
RUN="$BUNDLE/run_cxsv.sh"; SIF="$BUNDLE/container/cxsv-pipeline.sif"
RT=$(command -v singularity || command -v apptainer)
REF=/bundle/reference/hs1.fa; T="${SAMTOOLS_THREADS:-8}"
IN="$W/input"; mkdir -p "$IN"
exec > >(tee -a "$W/pos_control_full_genome.log") 2>&1
hdr() { echo; echo "== $* ($(date '+%F %T'))"; }
# in_sif CMD: bash -c CMD inside the bundle's SIF, network off, bundle read-only
in_sif() {
    "$RT" exec --containall --net --network none --bind "$BUNDLE:/bundle:ro" --bind "$W:$W" \
        --bind "$PC:$PC:ro" --env TMPDIR="$W" "$SIF" bash -c \
        "export REF_PATH=/bundle/reference/ref_cache/%2s/%2s/%s REF_CACHE=/bundle/reference/ref_cache/%2s/%2s/%s; $1"
}
cram_version() { local b; b=$(head -c 6 "$1" | od -An -tu1); set -- $b
                 [ "$1" = 67 ] && [ "$2" = 82 ] && [ "$3" = 65 ] && [ "$4" = 77 ] && echo "$5.$6"; }
total_records() { in_sif "samtools view -c -@ $T '$1'"; }   # all records, unmapped included

hdr "0. provenance"
echo "bundle $CXSV_BUNDLE_ID commit $CXSV_COMMIT release $CXSV_RELEASE_ID SIF $CXSV_SIF_SHA256"
echo "input $POS_CRAM ($(stat -c %s "$POS_CRAM") bytes)"
if [ "${SKIP_MD5:-0}" != 1 ]; then
    want=$(awk '$2=="pos_control_3559.hs1.minimap2.cram"{print $1}' "$PC/release/MD5SUMS")
    got=$(md5sum "$POS_CRAM" | cut -c1-32)
    [ "$got" = "$want" ] && echo "release CRAM md5 $got matches MD5SUMS" || { echo "release CRAM md5 $got != $want"; exit 1; }
fi
(cd "$ADDON" && md5sum -c --quiet MD5SUMS) && echo "E13 add-on matches its MD5SUMS"
n_m5=0
for m in $(in_sif "samtools view -H '$POS_CRAM'" | grep '^@SQ' | grep -o 'M5:[0-9a-f]*' | cut -c4-); do
    [ -f "$BUNDLE/reference/ref_cache/${m:0:2}/${m:2:2}/${m:4}" ] || { echo "M5 $m missing from the bundle cache"; exit 1; }
    n_m5=$((n_m5+1))
done
echo "all $n_m5 @SQ M5 of the input are in the bundle's MD5 cache"

hdr "1. CRAM version"
SRC="$POS_CRAM"; v=$(cram_version "$POS_CRAM"); echo "input is CRAM $v"
if [ "$v" = "3.1" ]; then
    if [ ! -f "$IN/pos_control.v30.cram.cram30.ok" ]; then
        in_sif "samtools view -@ $T -C --output-fmt-option version=3.0 -T $REF --write-index \
                -o '$IN/pos_control.v30.cram##idx##$IN/pos_control.v30.cram.crai' '$POS_CRAM'"
        a=$(total_records "$POS_CRAM"); b=$(total_records "$IN/pos_control.v30.cram")
        [ "$a" = "$b" ] || { echo "transcode record count mismatch $a vs $b"; exit 1; }
        printf 'source\t%s\nsource_cram_version\t3.1\ncram_version\t3.0\nrecords\t%s\n' "$POS_CRAM" "$a" > "$IN/pos_control.v30.cram.cram30.ok"
    fi
    SRC="$IN/pos_control.v30.cram"; echo "transcoded CRAM 3.1 -> 3.0 (lossless, records $(cut -f2 "$IN/pos_control.v30.cram.cram30.ok" | tail -1))"
else
    echo "no transcode needed"
fi

hdr "2. E13 splice (whole genome)"
OUT="$IN/pos_control.E13.cram"
if [ ! -f "$OUT.splice.ok" ]; then
    t=$(date +%s)
    in_sif "samtools view -@ $T -L '$ADDON/reads/E13.replace.bed' '$SRC' | cut -f1 | sort -u > '$IN/E13.dropped_names.txt'"
    in_sif "samtools view -@ $T -u -N '^$IN/E13.dropped_names.txt' -U '$IN/E13.dropped_records.bam' '$SRC' \
            | samtools merge -@ $T --reference $REF -O cram --output-fmt-option version=3.0 -o '$OUT.part' - \
                  '$ADDON/reads/pos_control.E13.sim.bam' && samtools index -@ $T '$OUT.part' '$OUT.part.crai'"
    n_in=$(total_records "$SRC"); n_drop=$(in_sif "samtools view -c '$IN/E13.dropped_records.bam'")
    n_sim=$(in_sif "samtools view -c '$ADDON/reads/pos_control.E13.sim.bam'"); n_out=$(total_records "$OUT.part")
    echo "records: input $n_in, dropped $n_drop ($(wc -l < "$IN/E13.dropped_names.txt") read names), simulated $n_sim, output $n_out"
    [ "$n_out" = $((n_in - n_drop + n_sim)) ] || { echo "splice record count mismatch"; exit 1; }
    mv "$OUT.part.crai" "$OUT.crai"; mv "$OUT.part" "$OUT"
    printf 'source\t%s\nsource_md5\t%s\nreplace_bed\t%s\ndropped_names\t%s\ndropped_records\t%s\nsimulated_records\t%s\ninput_records\t%s\noutput_records\t%s\ncram_version\t%s\nseconds\t%s\n' \
        "$SRC" "${got:-not checked}" "$ADDON/reads/E13.replace.bed" "$(wc -l < "$IN/E13.dropped_names.txt")" \
        "$n_drop" "$n_sim" "$n_in" "$n_out" "$(cram_version "$OUT")" "$(( $(date +%s)-t ))" > "$OUT.splice.ok"
fi
cat "$OUT.splice.ok"

hdr "3. patient run from the bundle (Singularity, network off)"
printf 'sample_id\talignment\tindividual_id\tplatform\npos_control\t%s\tpos_control_3559\tONT\n' "$OUT" > "$IN/patient.tsv"
cat > "$IN/pos_control_full.yaml" <<YAML
patient_runs:
  pos_control: {samples: "$IN/patient.tsv"}
threads: {samtools: 4, sv_callers: $CORES}
YAML
t=$(date +%s)
"$RUN" patient --runtime "$(basename "$RT")" --workdir "$W/run" --config "$IN/pos_control_full.yaml" \
    --bind "$IN" --cores "$CORES" > "$W/patient_run.log" 2>&1 || { tail -40 "$W/patient_run.log"; exit 1; }
echo "patient run finished: $(( $(date +%s)-t )) s"
S="$W/run/results/patient/samples/pos_control/pos_control"
sed 's/^/  /' "$W/run/results/patient/panel_release_check.txt"
grep -E '^(this run|panel|status)' "$W/run/results/patient/panel_compatibility.txt" | sed 's/^/  /'

hdr "4. score against the truth set (E01-E12 release truth + E13 add-on truth)"
in_sif "python /bundle/repo/tests/positive_control/score_pos_control.py '$PC/release/truth' '$S' \
        /bundle/repo/tests/positive_control/events.yaml '$W/score' --extra-truth '$ADDON/truth' \
        --known-failures /bundle/repo/tests/positive_control/known_failures.tsv" > "$W/score.txt" 2>&1 || true
tail -25 "$W/score.txt"

hdr "5. the bundle's panel release is unchanged"
"$RUN" verify-release | grep -E 'manifest sha256|files verified|status'
echo "done"
