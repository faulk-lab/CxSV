#!/usr/bin/env bash
# =============================================================================
# Reuse of caller outputs under a newer bundle (tests/lab/lib/reuse_callers.sh,
# used by tests/lab/pos_control_rerun_downstream.sh), on a toy workflow with the
# patient workflow's rule names and dependency shape:
#   md5_manifest.tsv + CRAM -> check_cram_reference -> run_{sniffles,cutesv,svim}
#   -> filter_* -> merge_sv_callers -> sample_vcf_to_bed -> all
#
#   bash tests/synthetic/check_rerun_touch.sh [WORKDIR]
#
# The finished run is used as the rerun script sees it: without .snakemake metadata.
#   1. control   a newer, identical manifest makes the plan re-run every caller
#   2. control   touching refcheck.ok (the first proposed fix) still re-runs them
#   3. fix       same_reference_manifest passes; touch_upstream; the plan with
#                --forcerun sample_vcf_to_bed holds no upstream rule, and the
#                caller outputs keep their content
#   4. negative  a manifest with other content is refused by same_reference_manifest
#   5. negative  a touch plan holding a downstream rule is refused (nothing touched)
# Needs snakemake (cxsv-pipeline env). Exit status 0 = all checks passed.
# =============================================================================
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
W="$(mkdir -p "${1:-$REPO/.rerun_touch_test}" && cd "${1:-$REPO/.rerun_touch_test}" && pwd)"
# shellcheck source=tests/lab/lib/reuse_callers.sh
source "$REPO/tests/lab/lib/reuse_callers.sh"
rm -rf "$W/run" "$W/bundle_old" "$W/bundle_new"; mkdir -p "$W/run" "$W/bundle_old" "$W/bundle_new"
fail() { echo "FAIL  $*"; exit 1; }
upstream_re="^(localrule|rule) ($(IFS='|'; echo "${CXSV_UPSTREAM_RULES[*]}")):"

printf 'name\tlength\tmd5\nchr1\t100\tabc\n' > "$W/bundle_old/md5_manifest.tsv"
echo reads > "$W/run/input.cram"
cat > "$W/run/Snakefile" <<'SMK'
MANIFEST = config["manifest"]
CALLERS = ["sniffles", "cutesv", "svim"]
rule all:
    input: "bed/merged_sv.bed"
rule check_cram_reference:
    input: MANIFEST, "input.cram"
    output: "refcheck.ok"
    shell: "echo ok > {output}"
for c in CALLERS:
    rule:
        name: f"run_{c}"
        input: "input.cram", "refcheck.ok"
        output: f"callers/{c}.vcf"
        shell: "echo call-$RANDOM > {output}"
    rule:
        name: f"filter_{c}"
        input: f"callers/{c}.vcf"
        output: f"callers/{c}_filtered.vcf"
        shell: "cp {input} {output}"
rule merge_sv_callers:
    input: expand("callers/{c}_filtered.vcf", c=CALLERS)
    output: "merged_sv.vcf"
    shell: "cat {input} > {output}"
rule sample_vcf_to_bed:
    input: "merged_sv.vcf"
    output: "bed/merged_sv.bed"
    shell: "cp {input} {output}"
SMK
SM_BASE=(snakemake -s "$W/run/Snakefile" --directory "$W/run" --cores 1 --rerun-triggers mtime)
plan() { "${SM_BASE[@]}" --config manifest="$1" -n "${@:2}" 2>&1; }

"${SM_BASE[@]}" --config manifest="$W/bundle_old/md5_manifest.tsv" -q > "$W/first_run.log" 2>&1 || { cat "$W/first_run.log"; fail "first run"; }
rm -rf "$W/run/.snakemake"   # the rerun script copies the run without .snakemake (no metadata)
sha_before=$(cat "$W/run/callers/"*.vcf "$W/run/merged_sv.vcf" | sha256sum | cut -c1-64)
sleep 1.1; cp "$W/bundle_old/md5_manifest.tsv" "$W/bundle_new/md5_manifest.tsv"   # newer bundle: same content, newer mtime
NEW="$W/bundle_new/md5_manifest.tsv"

plan "$NEW" > "$W/plan_control.txt"
grep -Eq '^rule run_sniffles:' "$W/plan_control.txt" || fail "control: the newer manifest should re-run the callers"
echo "PASS  1. control: a newer identical manifest re-runs the callers"

touch "$W/run/refcheck.ok"
plan "$NEW" > "$W/plan_touch_refcheck.txt"
grep -Eq '^rule run_sniffles:' "$W/plan_touch_refcheck.txt" || fail "control: touching refcheck.ok was expected not to help"
echo "PASS  2. control: touching refcheck.ok alone still re-runs the callers"

same_reference_manifest "$W/bundle_old/md5_manifest.tsv" "$NEW" || fail "identical manifests not recognised"
touch_upstream merged_sv.vcf "${SM_BASE[@]}" --config manifest="$NEW" > "$W/touch_upstream.txt" 2>&1 || { cat "$W/touch_upstream.txt"; fail "touch_upstream"; }
plan "$NEW" --forcerun sample_vcf_to_bed > "$W/plan_fixed.txt"
if grep -E "$upstream_re" "$W/plan_fixed.txt"; then fail "after touch_upstream the plan still holds an upstream rule"; fi
grep -Eq '^rule sample_vcf_to_bed:' "$W/plan_fixed.txt" || fail "downstream rule not planned"
sha_after=$(cat "$W/run/callers/"*.vcf "$W/run/merged_sv.vcf" | sha256sum | cut -c1-64)
[ "$sha_before" = "$sha_after" ] || fail "caller outputs changed"
echo "PASS  3. fix: after touch_upstream only downstream rules are planned; caller outputs unchanged"

printf 'name\tlength\tmd5\nchr1\t100\tdef\n' > "$W/bundle_other.tsv"
if same_reference_manifest "$W/bundle_old/md5_manifest.tsv" "$W/bundle_other.tsv"; then fail "a different manifest was accepted"; fi
if same_reference_manifest "$W/bundle_old/md5_manifest.tsv" "$W/missing.tsv"; then fail "a missing manifest was accepted"; fi
echo "PASS  4. negative: a different or missing manifest is refused"
touch "$W/run/merged_sv.vcf"; sleep 1.1; touch "$W/run/input.cram"
if touch_upstream bed/merged_sv.bed "${SM_BASE[@]}" --config manifest="$NEW" > "$W/touch_refused.txt" 2>&1; then
    fail "a touch plan with sample_vcf_to_bed was accepted"; fi
grep -q 'outside CXSV_UPSTREAM_RULES' "$W/touch_refused.txt" || fail "refusal message missing"
echo "PASS  5. negative: a touch plan with a rule outside CXSV_UPSTREAM_RULES touches nothing"
echo "all rerun-touch checks passed"
