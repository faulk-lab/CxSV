#!/usr/bin/env bash
# =============================================================================
# Regression tests of two plot scripts (small fixes 2026-11).
#
#   bash tests/synthetic/check_plot_scripts.sh [WORKDIR]
#
# 1. plot_complexity_class.R must draw the same PDF under every locale. Renders one small master table (class names chosen so that C and natural-
# language collation order them differently: "CN_plus_orientation_change" vs
# "Chromothripsis_like") under each available locale among C, POSIX, C.UTF-8
# and en_US.UTF-8, and compares the PDFs without their render timestamps
# (tests/lab/compare_patient_runs.py rules). Bundle validation check 7 failed
# on 2026-11-03 because the heatmap rows followed LC_COLLATE (host vs image).
# Needs Rscript with ggplot2/dplyr/tidyr/patchwork (the
# cxsv-pipeline env). Exit 0 = passed.
# =============================================================================
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
W="${1:-$(mktemp -d)}"; mkdir -p "$W"; W="$(cd "$W" && pwd)"
T="$W/plot_checks"; rm -rf "$T"; mkdir -p "$T"

printf 'Locus_ID\tCmplx_Class\tSV_Types\tCxSV_Criteria\n' > "$T/master.tsv"
i=0
for row in "CN_plus_orientation_change|DEL,INV|C1_min_breakpoints,C4_CN_plus_orient" \
           "Chromothripsis_like|DEL,DUP,INV,BND|C1_min_breakpoints,C5_chromoanagenesis_cluster" \
           "Complex_inversion|INV,DUP|C1_min_breakpoints" "complex_lowercase|DEL|C1_min_breakpoints" \
           "DEL+INV|DEL,INV|C1_min_breakpoints" "BND_cluster|BND,TRA|C1_min_breakpoints"; do
    IFS='|' read -r cls types crit <<< "$row"
    for _ in 1 2; do i=$((i+1)); printf 'L%d\t%s\t%s\t%s\n' "$i" "$cls" "$types" "$crit" >> "$T/master.tsv"; done
done

avail="$(locale -a 2>/dev/null | tr 'A-Z' 'a-z' || true)"
used=()
for loc in C POSIX C.UTF-8 en_US.UTF-8; do
    l="$(echo "$loc" | tr 'A-Z' 'a-z' | sed 's/utf-8/utf8/')"
    case "$loc" in C|POSIX) ;; *) echo "$avail" | sed 's/utf-8/utf8/' | grep -qx "$l" || { echo "skip $loc (not installed)"; continue; } ;; esac
    out="$T/plot_${loc}.pdf"
    LC_ALL="$loc" LANG="$loc" Rscript "$REPO/workflow/scripts/plot_complexity_class.R" "$T/master.tsv" "$out" \
        > "$T/plot_${loc}.log" 2>&1 || { echo "FAIL render under $loc"; cat "$T/plot_${loc}.log"; exit 1; }
    used+=("$loc")
done
[ "${#used[@]}" -ge 2 ] || { echo "FAIL fewer than two locales available: ${used[*]}"; exit 1; }

python - "$T" "${used[@]}" <<'PY'
import re, sys
d, locs = sys.argv[1], sys.argv[2:]
vol = re.compile(rb"/CreationDate \(D:[^)]*\)|/ModDate \(D:[^)]*\)")
norm = {l: vol.sub(b"", open(f"{d}/plot_{l}.pdf", "rb").read()) for l in locs}
ref = locs[0]
bad = [l for l in locs[1:] if norm[l] != norm[ref]]
if bad:
    sys.exit(f"FAIL plot_complexity_class.R output depends on the locale: {ref} differs from {bad}")
print(f"PASS plot_complexity_class.R: identical PDF (minus timestamps) under {', '.join(locs)}")
PY

