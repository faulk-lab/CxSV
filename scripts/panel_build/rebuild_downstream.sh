#!/usr/bin/env bash
# =============================================================================
# rebuild_downstream.sh — rebuild a panel release from an existing build's
# upstream outputs, re-running only the rules downstream of a calling change
# (panel-v1.1 = panel-v1.0 + roadmap item 13; docs/PANEL_BUILD.md,
# "Downstream-only rebuild").
#
#   bash rebuild_downstream.sh OLD_BUILD NEW_BUILD CORES [--dry-run | --resume]
#
# --resume continues an interrupted or partly failed run in the same tree:
# no --forcerun (Snakemake's own metadata in NEW_BUILD/aggregate now records
# the jobs that ran), --rerun-incomplete, and no check that planned outputs are
# absent; the held-rule check and the post-run checks still apply.
#
# Expects (made by the steps in docs/PANEL_BUILD.md):
#   NEW_BUILD/code               git archive of the release commit + COMMIT
#   NEW_BUILD/aggregate.yaml     OLD_BUILD/aggregate.yaml with
#                                panel_results_root = NEW_BUILD/results/panel
#                                and the new release description
#   NEW_BUILD/results/panel      upstream outputs of OLD_BUILD, hard-linked by
#                                link_upstream.py (inventory in
#                                NEW_BUILD/provenance/reused_upstream.tsv)
# Environment: CXSV_SIF (pinned image; its bundle.env gives the image ids for
# release.json), CXSV_REF_CACHE (htslib MD5 cache).
#
# The Snakemake work directory (its .snakemake metadata) is NEW_BUILD/aggregate,
# or CXSV_REBUILD_WORKDIR: on node-local disk the per-job metadata writes are
# much faster than on the CIFS NAS (about 50 jobs/min there for 12,975 short jobs).
#
# Held rules (never run; their outputs are the reused files): the callers and
# their filters, merge_sv_callers, check_cram_reference, build_ref_cache,
# tool_versions, callable_mask, panel_callability_masks and
# panel_sniffles_population, by --allowed-rules = every other rule. The
# aggregate work directory has no per-sample Snakemake metadata, so
# --forcerun sample_vcf_to_bed calling_fingerprint is what makes every
# downstream rule run. The dry run must plan no held-rule job and no planned
# output may be a hard link (so no job writes to an inode shared with the old
# build); after the run
# the reused files are re-checked against the inventory and the release is
# verified.
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # link_upstream.py lives next to this script
OLD="$1" NEW="$2" CORES="$3" DRY=false RESUME=false
[ "${4:-}" = --dry-run ] && DRY=true
[ "${4:-}" = --resume ] && RESUME=true
CODE="$NEW/code"
: "${CXSV_SIF:?set CXSV_SIF}" "${CXSV_REF_CACHE:=/mnt/sargasso1/CxSV/refs/ref_cache}"
export CXSV_SIF CXSV_REF_CACHE
# shellcheck source=/dev/null
source "$CODE/scripts/panel_build/lib.sh"
[ "$(cat "$CODE/COMMIT")" = "$(cat "$CODE/workflow/COMMIT" 2>/dev/null)" ] \
    || cp "$CODE/COMMIT" "$CODE/workflow/COMMIT"     # release.json source.commit (cxsv_release.source_info)
HELD='run_sniffles|run_cutesv|run_svim|filter_sniffles|filter_cutesv|filter_svim|merge_sv_callers|check_cram_reference|build_ref_cache|tool_versions|callable_mask|panel_callability_masks|panel_sniffles_population'
ulimit -Sn "$(ulimit -Hn)"
WD="${CXSV_REBUILD_WORKDIR:-$NEW/aggregate}" TMP="$HOME/scratch/cxsv_panel_build/tmp/aggregate-$(basename "$NEW")"
mkdir -p "$WD" "$NEW/provenance"
case "$WD" in "$CXSV_NAS"/*) ;; *) export CXSV_EXTRA_BINDS="${CXSV_EXTRA_BINDS:-} $WD" ;; esac
sm() { cxsv_container "$WD" "$TMP" snakemake -s "$CODE/workflow/panel/Snakefile" --directory "$WD" \
        --configfile "$NEW/aggregate.yaml" --cores "$CORES" "$@"; }
sm --list-rules 2>/dev/null | grep -vE '^[[:space:]]' | awk '{print $1}' | sort -u > "$WD/rules.txt"
grep -v -x -E "$HELD" "$WD/rules.txt" > "$WD/allowed.txt"
ARGS=(--allowed-rules $(cat "$WD/allowed.txt") --keep-going --printshellcmds)
if $RESUME; then ARGS+=(--rerun-incomplete); DRYOUT="$WD/dryrun.resume.txt"
else ARGS+=(--forcerun sample_vcf_to_bed calling_fingerprint); DRYOUT="$WD/dryrun.txt"; fi
sm -n "${ARGS[@]}" > "$DRYOUT" 2>&1
awk '/^Job stats/{f=1} f{print} f && /^$/{exit}' "$DRYOUT"
held=$(awk -v re="^($HELD)\$" '/^Job stats/{f=1} f && $1 ~ re {s+=$2} f && /^$/{exit} END{print s+0}' "$DRYOUT")
[ "$held" = 0 ] || { echo "dry run plans $held held-rule jobs; stop"; exit 1; }
$RESUME || python3 - "$DRYOUT" "$HERE" <<'PY'
import os, sys
sys.path.insert(0, sys.argv[2])
from link_upstream import planned_paths
planned, _ = planned_paths(sys.argv[1])
there = [p for p in planned if os.path.exists(p)]
linked = sorted(p for p in there if os.path.isfile(p) and os.stat(p).st_nlink > 1)
print(f"planned outputs/logs/benchmarks: {len(planned)}; present (from an earlier attempt): {len(there)}; hard-linked: {len(linked)}")
if linked:
    print("\n".join(linked[:20])); sys.exit("a planned output is a hard link (shared with the old build); stop")
PY
$DRY && exit 0
LOG="$WD/aggregate.log"; $RESUME && LOG="$WD/aggregate.resume.log"
sm "${ARGS[@]}" > "$LOG" 2>&1 || { echo "snakemake failed: $LOG"; tail -n 30 "$LOG"; exit 1; }
R="$NEW/results/panel"
bad=$(grep -lE 'Too many open files|Errno 24|ERROR sniffles\.worker' "$R"/panel/logs/*.log "$R"/cohorts/*/logs/*.log 2>/dev/null || true)
[ -z "$bad" ] || { echo "worker errors in: $bad"; exit 1; }
python3 "$HERE/link_upstream.py" --verify "$NEW/provenance/reused_upstream.tsv"
cxsv_container "$WD" "$TMP" python "$CODE/workflow/scripts/verify_panel_release.py" "$R/release" \
    --out "$NEW/provenance/verify_panel_release.txt"
echo "rebuild done: $R/release"
