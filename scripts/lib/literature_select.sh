# shellcheck shell=bash
# =============================================================================
# literature_select.sh — choose the curated literature set that make_bundle.sh
# ships (roadmap item 11; docs/RELEASE_AND_BUNDLE.md). Sourced, not executed.
#
#   cxsv_select_literature RESOURCES LITERATURE NO_LITERATURE
#
#   RESOURCES      lab resources root (<RESOURCES>/literature/lit-v*/ are candidates)
#   LITERATURE     explicit set directory (make_bundle.sh --literature), or ""
#   NO_LITERATURE  1 = build without a literature set (make_bundle.sh --no-literature)
#
# Prints the chosen directory (nothing with NO_LITERATURE=1) and returns 0;
# otherwise prints the reason to stderr and returns 1.
#
# Only DIRECTORIES named lit-v* that hold a MANIFEST.tsv are candidates; the
# newest by version sort wins. Until panel-v1.2 the selection was
# `ls -d <RESOURCES>/literature/lit-v* | sort -V | tail -1`, which also matched
# the file lit-v1.0.build.log, sorted it last and so skipped the literature set
# without a message (bundles cxsv-bundle-20261103-33c075a22071,
# cxsv-bundle-20261106-169503679a60 and cxsv-bundle-20261108-c40e26528d5f carry
# no literature set).
# =============================================================================

cxsv_select_literature() {
    local resources="$1" literature="$2" no_literature="${3:-0}"
    if [ "$no_literature" = 1 ]; then
        if [ -n "$literature" ]; then
            echo "--literature and --no-literature are mutually exclusive" >&2
            return 1
        fi
        return 0
    fi
    if [ -n "$literature" ]; then
        if [ ! -d "$literature" ]; then
            echo "literature set $literature: not a directory" >&2; return 1
        fi
        if [ ! -f "$literature/MANIFEST.tsv" ]; then
            echo "literature set $literature: no MANIFEST.tsv" >&2; return 1
        fi
        (cd "$literature" && pwd)
        return 0
    fi
    local root="$resources/literature" chosen=""
    if [ -d "$root" ]; then
        chosen=$(find "$root" -mindepth 1 -maxdepth 1 -type d -name 'lit-v*' -print 2>/dev/null \
                 | while read -r d; do [ -f "$d/MANIFEST.tsv" ] && echo "$d"; done \
                 | sort -V | tail -n 1)
    fi
    if [ -z "$chosen" ]; then
        echo "no literature set found: no directory $root/lit-v*/ with a MANIFEST.tsv." >&2
        echo "Pass --literature DIR, or --no-literature to build a bundle without one" >&2
        echo "(such a bundle fails check 9 of tests/lab/validate_bundle.sh)." >&2
        return 1
    fi
    echo "$chosen"
}

# cxsv_verify_literature DIR: every file listed in DIR/MANIFEST.tsv (columns
# file, bytes, sha256; header line) has its sha256. Returns 0 when all match.
cxsv_verify_literature() {
    local dir="$1"
    [ -f "$dir/MANIFEST.tsv" ] || { echo "$dir: no MANIFEST.tsv" >&2; return 1; }
    local n; n=$(awk -F'\t' 'NR>1 && NF>=3' "$dir/MANIFEST.tsv" | wc -l)
    [ "$n" -gt 0 ] || { echo "$dir/MANIFEST.tsv lists no files" >&2; return 1; }
    awk -F'\t' 'NR>1 && NF>=3 {print $3"  "$1}' "$dir/MANIFEST.tsv" | (cd "$dir" && sha256sum --quiet -c -)
}
