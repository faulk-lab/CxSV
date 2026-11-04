#!/usr/bin/env bash
# =============================================================================
# The container environment that container/run_cxsv.sh gives a run (small
# fixes 2026-11; docs/RELEASE_AND_BUNDLE.md).
#
#   bash tests/lab/check_container_env.sh BUNDLE WORKDIR [RUNTIME...]
#
# Uses the repository's run_cxsv.sh with an existing bundle's image: WORKDIR/
# minibundle holds this checkout's run_cxsv.sh, BUNDLE's bundle.env and links
# to BUNDLE's SIF and Docker archive. RUNTIME: singularity, apptainer, docker
# (default: singularity and docker, those installed). Per runtime, inside
# `run_cxsv.sh shell`:
#   tmp   TMPDIR is /work/.tmp; Python tempfile and R tempdir() use it; 200 MB
#         written through Python, through R and to /tmp itself all succeed
#         (Singularity's --containall /tmp is a 64 MB tmpfs); the Python and
#         /tmp files appear in WORKDIR/run_<runtime>/.tmp on the host
#   ref   REF_PATH and REF_CACHE point at /bundle/reference/ref_cache; samtools
#         decodes every base of a CRAM (tests/fixtures/cram31; its header UR
#         path is not visible) with no --reference, i.e. through REF_PATH
#         alone, and fails once REF_PATH/REF_CACHE are unset (control)
# Writes WORKDIR/container_env.log. Exit 0 = all passed.
# =============================================================================
set -uo pipefail
BUNDLE="$(cd "$1" && pwd)"; W="$(mkdir -p "$2" && cd "$2" && pwd)"; shift 2
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
RUNTIMES=("$@")
if [ ${#RUNTIMES[@]} -eq 0 ]; then
    for r in singularity docker; do command -v "$r" >/dev/null && RUNTIMES+=("$r"); done
fi
exec > >(tee "$W/container_env.log") 2>&1
FAILS=0
pass() { echo "PASS  $*"; }
fail() { echo "FAIL  $*"; FAILS=$((FAILS+1)); }

MB="$W/minibundle"; rm -rf "$MB"; mkdir -p "$MB/container" "$MB/reference"
cp "$REPO/container/run_cxsv.sh" "$MB/"; chmod +x "$MB/run_cxsv.sh"; cp "$BUNDLE/bundle.env" "$MB/"
for f in cxsv-pipeline.sif cxsv-pipeline.docker.tar.gz; do ln -s "$BUNDLE/container/$f" "$MB/container/$f"; done
# the fixture's reference as the "bundle" MD5 cache; the CRAM itself goes into each run directory
python3 "$REPO/workflow/scripts/build_ref_cache.py" "$REPO/tests/fixtures/cram31/ref.fa" "$MB/reference/ref_cache" > /dev/null
echo "repo $(git -C "$REPO" rev-parse --short HEAD 2>/dev/null || echo '?'); image from $BUNDLE"

cat > "$W/probe.sh" <<'PROBE'
echo "TMPDIR=$TMPDIR"
python -c "import tempfile; print('python_tmp', tempfile.gettempdir())"
Rscript -e 'cat("R_tmp", dirname(tempdir()), "\n")'
python - <<'PY'
import os, tempfile
fd, p = tempfile.mkstemp(prefix="py_probe_")
with os.fdopen(fd, "wb") as fh:
    for _ in range(200):
        fh.write(b"\0" * (1 << 20))
print("python_write_ok", os.path.basename(p), os.path.getsize(p))
PY
Rscript -e 'f <- tempfile("R_probe_"); writeBin(raw(200 * 2^20), f); cat("R_write_ok", basename(f), file.size(f), "\n")'
dd if=/dev/zero of=/tmp/slash_tmp_probe bs=1M count=200 status=none && echo "slash_tmp_write_ok $(stat -c %s /tmp/slash_tmp_probe)"
echo "tmp_fs $(df -P /tmp | tail -1 | awk '{print $1, $2}')"
echo "REF_PATH=$REF_PATH"; echo "REF_CACHE=$REF_CACHE"
# decode the sequences (samtools view -c would not need the reference); the reads match the reference,
# so CRAM stores no bases and every base must come from the cache
dec() { if "$@" samtools view -o /tmp/dec.sam /work/fixture.cram 2>/dev/null; then
            awk '{n += length($10)} END {print NR, n}' /tmp/dec.sam; else echo failed; fi; rm -f /tmp/dec.sam; }
echo "cram_decoded $(dec env)"
echo "cram_decoded_without_ref_path $(dec env -u REF_PATH -u REF_CACHE)"
PROBE

for rt in "${RUNTIMES[@]}"; do
    echo; echo "== $rt ($(date '+%F %T'))"
    R="$W/run_$rt"; rm -rf "$R"; mkdir -p "$R"
    cp "$REPO/tests/fixtures/cram31/reads.v31.cram" "$R/fixture.cram"
    cp "$REPO/tests/fixtures/cram31/reads.v31.cram.crai" "$R/fixture.cram.crai"
    out="$("$MB/run_cxsv.sh" shell --runtime "$rt" --workdir "$R" < "$W/probe.sh" 2>&1)"; rc=$?
    echo "$out" | sed 's/^/  | /'
    [ $rc -eq 0 ] || { fail "$rt: probe exited $rc"; continue; }
    echo "$out" | grep -qx "TMPDIR=/work/.tmp" && pass "$rt: TMPDIR=/work/.tmp" || fail "$rt: TMPDIR"
    echo "$out" | grep -qx "python_tmp /work/.tmp" && pass "$rt: Python tempfile uses it" || fail "$rt: Python tempdir"
    echo "$out" | grep -qE "^R_tmp /work/.tmp *$" && pass "$rt: R tempdir() uses it" || fail "$rt: R tempdir"
    for k in python R slash_tmp; do
        echo "$out" | grep -qE "^${k}_write_ok .*209715200 *$" && pass "$rt: 200 MB written via $k" || fail "$rt: 200 MB via $k"
    done
    # (R removes its session directory, and the R file with it, when it exits)
    n=$(find "$R/.tmp" -maxdepth 1 -type f \( -name 'py_probe_*' -o -name slash_tmp_probe \) -size +199M | wc -l)
    [ "$n" -eq 2 ] && pass "$rt: the Python and /tmp files are in the run directory on the host ($R/.tmp)" \
        || fail "$rt: $n of 2 probe files found in $R/.tmp"
    c="/bundle/reference/ref_cache/%2s/%2s/%s"
    echo "$out" | grep -qx "REF_PATH=$c" && echo "$out" | grep -qx "REF_CACHE=$c" \
        && pass "$rt: REF_PATH and REF_CACHE = $c" || fail "$rt: REF_PATH/REF_CACHE"
    echo "$out" | grep -qx "cram_decoded 364 728000" && pass "$rt: CRAM decoded through REF_PATH alone (364 reads, 728,000 bases)" \
        || fail "$rt: CRAM not decoded through REF_PATH"
    echo "$out" | grep -qx "cram_decoded_without_ref_path failed" && pass "$rt: control: decoding fails without REF_PATH" \
        || fail "$rt: control: CRAM decoded without REF_PATH (reference found elsewhere?)"
    rm -rf "$R/.tmp"/*
done

echo; echo "failures: $FAILS"
[ "$FAILS" -eq 0 ]
