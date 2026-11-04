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
# 2. plot_cxsv_venn.R must keep calls of different samples apart: two samples
#    whose cxsv_only.vcf files carry the same record ids (each caller numbers
#    its calls per sample) must give 2 calls per caller, not 1, and the log
#    must name the samples <cohort>/<sample>.
# Needs Rscript with ggplot2/dplyr/tidyr/patchwork/ggVennDiagram (the
# cxsv-pipeline env). Exit 0 = both passed.
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

# ── 2. Venn call ids are per sample ──────────────────────────────────────────
V="$T/venn"; mkdir -p "$V"
printf 'Chr\tStart\tEnd\tIsCxSV\nchr1\t1000\t5000\tYES\n' > "$V/cxsv_summary.tsv"
vcfs=()
for s in sA sB; do
    d="$V/samples/COH/$s/cxsv"; mkdir -p "$d"
    { printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
      printf 'chr1\t1500\tsame_id_1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=2500;SUPP_VEC=111\n'; } > "$d/cxsv_only.vcf"
    vcfs+=("$d/cxsv_only.vcf")
done
list="$(IFS=,; echo "${vcfs[*]}")"
Rscript "$REPO/workflow/scripts/plot_cxsv_venn.R" "$list" "$V/cxsv_summary.tsv" "$V/venn.pdf" > "$V/venn.log" 2>&1 \
    || { echo "FAIL plot_cxsv_venn.R"; cat "$V/venn.log"; exit 1; }
if grep -q "Caller sets: CuteSV=2  SVIM=2  Sniffles2=2" "$V/venn.log"; then
    echo "PASS plot_cxsv_venn.R: same record id in two samples counted as two calls"
else
    echo "FAIL plot_cxsv_venn.R: call ids collide across samples"; cat "$V/venn.log"; exit 1
fi
# two inputs that resolve to one sample name must be refused, not merged
mkdir -p "$V/dup"; cp "${vcfs[0]}" "$V/dup/x_merged_sv.vcf"; mkdir -p "$V/dup2"; cp "${vcfs[0]}" "$V/dup2/x_merged_sv.vcf"
if Rscript "$REPO/workflow/scripts/plot_cxsv_venn.R" "$V/dup/x_merged_sv.vcf,$V/dup2/x_merged_sv.vcf" \
        "$V/cxsv_summary.tsv" "$V/venn_dup.pdf" > "$V/venn_dup.log" 2>&1; then
    echo "FAIL plot_cxsv_venn.R accepted two inputs with one sample name"; exit 1
fi
grep -q "same sample name" "$V/venn_dup.log" && echo "PASS plot_cxsv_venn.R: duplicate sample names refused" \
    || { echo "FAIL unexpected error"; cat "$V/venn_dup.log"; exit 1; }
