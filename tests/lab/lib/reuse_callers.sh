# shellcheck shell=bash
# =============================================================================
# Reusing a finished run's caller outputs under a different bundle.
# Sourced by tests/lab/pos_control_rerun_downstream.sh; tested by
# tests/synthetic/check_rerun_touch.sh (CI) and tests/unit/test_reuse_callers.py.
#
# Problem: a newer bundle ships the same reference MD5 cache but with newer mtimes.
# Under --rerun-triggers mtime, check_cram_reference then re-runs, and every caller
# follows. Touching refcheck.ok, or forcing check_cram_reference, does not help:
# refcheck.ok is an input of every caller rule, so the callers still re-run.
#
# Fix: once the two md5_manifest.tsv files are shown to be byte-identical,
# `snakemake --touch <merged_sv.vcf>` marks the reference check, the callers,
# their filters and the merge as up to date, in DAG order, without running them.
# The caller VCFs keep their content; only mtimes change. (`--touch --until
# <rules>` is not used: Snakemake 9.23 stops it with "Out of jobs ready to be
# started" after the callers.)
# =============================================================================

# Rules whose outputs may be reused. Everything downstream of merge_sv_callers is
# re-made by the rerun (--forcerun sample_vcf_to_bed calling_fingerprint).
CXSV_UPSTREAM_RULES=(check_cram_reference run_sniffles run_cutesv run_svim
                     filter_sniffles filter_cutesv filter_svim merge_sv_callers)

# same_reference_manifest A B: status 0 when both files exist and are byte-identical.
same_reference_manifest() {
    [ -f "$1" ] && [ -f "$2" ] || return 1
    [ "$(sha256sum < "$1" | cut -c1-64)" = "$(sha256sum < "$2" | cut -c1-64)" ]
}

# touch_upstream TARGET SNAKEMAKE_COMMAND...: touch the jobs the plan would run to
# make TARGET (the run's merged_sv.vcf). The touch plan is listed first and must
# hold only CXSV_UPSTREAM_RULES; otherwise nothing is touched (status 1). The
# command must not contain --forcerun.
touch_upstream() {
    local target="$1"; shift
    local plan rules r ok
    plan=$("$@" --touch -n "$target" 2>&1) || { echo "$plan"; return 1; }
    echo "$plan"
    rules=$(echo "$plan" | sed -nE 's/^(local)?rule ([A-Za-z0-9_]+):$/\2/p' | sort -u)
    for r in $rules; do
        ok=0; for u in "${CXSV_UPSTREAM_RULES[@]}"; do [ "$r" = "$u" ] && ok=1; done
        [ $ok = 1 ] || { echo "touch plan holds rule $r outside CXSV_UPSTREAM_RULES; nothing touched"; return 1; }
    done
    "$@" --touch "$target"
}
