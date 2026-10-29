#!/usr/bin/env bash
# =============================================================================
# Validate an offline bundle (roadmap item 2; docs/RELEASE_AND_BUNDLE.md).
#
#   bash tests/lab/validate_bundle.sh BUNDLE WORKDIR STANDIN_DIR [NORMAL_RELEASE]
#
# STANDIN_DIR holds a public stand-in patient: standin.tsv (sample sheet,
# absolute alignment paths inside STANDIN_DIR) + its CRAM/BAM and index.
# NORMAL_RELEASE: the panel release as built on the lab fleet (default: the
# bundle's own copy). Needs the cxsv-pipeline conda env for step 6 only.
#
#   1. bundle integrity        run_cxsv.sh verify-bundle (sha256 of every file)
#   2. release integrity       run_cxsv.sh verify-release
#   3. network off             probe inside the container exactly as runs see it
#                              (Singularity --net --network none, Docker --network
#                              none): interfaces, TCP to 1.1.1.1:443, DNS, curl;
#                              host paths (repo clone, conda, NAS) invisible;
#                              /bundle read-only. Control: the same probe with
#                              --no-netns must reach the network.
#   4. patient run offline, Singularity/Apptainer (bundle only)
#   5. patient run offline, Docker (image loaded from the bundle's archive)
#   6. patient run outside the bundle: repository clone + conda env + lab paths
#   7. records identical: 4 vs 6 and 5 vs 6 (compare_patient_runs.py; files that
#      differ only by the run's absolute paths are reported as such)
#   8. tampered release refused: one byte changed / one unlisted file added
# Writes WORKDIR/validation.log and WORKDIR/compare_*.tsv. Exit 0 = all passed.
# =============================================================================
set -uo pipefail
BUNDLE="$(cd "$1" && pwd)"; W="$(mkdir -p "$2" && cd "$2" && pwd)"; STANDIN="$(cd "$3" && pwd)"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
source "$BUNDLE/bundle.env"
NORMAL_RELEASE="${4:-$BUNDLE/panel_release/$CXSV_RELEASE_ID}"
CORES="${CORES:-16}"
RUN="$BUNDLE/run_cxsv.sh"
FAILS=0
exec > >(tee "$W/validation.log") 2>&1
pass() { echo "PASS  $*"; }
fail() { echo "FAIL  $*"; FAILS=$((FAILS+1)); }
hdr()  { echo; echo "== $* ($(date '+%F %T'))"; }

hdr "0. context"
echo "host $(hostname -s); kernel $(uname -r); $(singularity --version 2>/dev/null || apptainer --version); $(docker --version)"
echo "bundle $BUNDLE"; echo "bundle id $CXSV_BUNDLE_ID; commit $CXSV_COMMIT; release $CXSV_RELEASE_ID (test: $CXSV_TEST_RELEASE)"
echo "image $CXSV_IMAGE $CXSV_IMAGE_ID; SIF sha256 $CXSV_SIF_SHA256"
echo "bundle size $(du -sb "$BUNDLE" | cut -f1) bytes ($(du -sh "$BUNDLE" | cut -f1)), $(wc -l < "$BUNDLE/MANIFEST.sha256") files in MANIFEST.sha256"
echo "stand-in: $(tail -n +2 "$STANDIN/standin.tsv" | cut -f1,2 | tr '\t' ' ')"

hdr "1. bundle integrity"
t=$(date +%s); if "$RUN" verify-bundle; then pass "every bundle file matches MANIFEST.sha256 ($(( $(date +%s)-t )) s)"; else fail "bundle checksums"; fi

hdr "2. panel release integrity"
if "$RUN" verify-release; then pass "panel release verifies against its manifest"; else fail "panel release"; fi

hdr "3. network off inside the container"
cat > "$W/probe.sh" <<'PROBE'
echo "interfaces: $(tail -n +3 /proc/net/dev | cut -d: -f1 | tr -d ' ' | tr '\n' ' ')"
python - <<'PY'
import socket
for host, port in (("1.1.1.1", 443), ("140.82.112.3", 443)):
    try:
        socket.create_connection((host, port), timeout=5).close(); print(f"tcp {host}:{port}: CONNECTED")
    except OSError as e:
        print(f"tcp {host}:{port}: failed ({e.__class__.__name__}: {e})")
try:
    print("dns github.com:", socket.gethostbyname("github.com"))
except OSError as e:
    print(f"dns github.com: failed ({e})")
PY
if command -v curl >/dev/null; then curl -sS -m 10 -o /dev/null https://github.com && echo "curl https://github.com: OK" || echo "curl https://github.com: failed (exit $?)"; fi
for p in /mnt/sargasso1 /home/minknow/miniforge3 /mnt/sargasso1/CxSV/repo; do [ -e "$p" ] && echo "host path visible: $p" || echo "host path absent: $p"; done
touch /bundle/.write_test 2>/dev/null && echo "/bundle WRITABLE" || echo "/bundle read-only"
echo "snakemake: $(command -v snakemake); CXSV_IMAGE=$CXSV_IMAGE; HOME=$HOME"
PROBE
for rt in singularity docker; do
    out=$("$RUN" shell --runtime "$rt" --workdir "$W/probe_$rt" < "$W/probe.sh" 2>&1); echo "--- $rt, as used for runs:"; echo "$out"
    if echo "$out" | grep -q "interfaces: lo $" && ! echo "$out" | grep -qE "CONNECTED|curl https://github.com: OK|dns github.com: [0-9]" \
       && ! echo "$out" | grep -q "host path visible" && echo "$out" | grep -q "/bundle read-only"; then
        pass "$rt: loopback only, no TCP/DNS/HTTPS, no host repo/conda/NAS, bundle read-only"
    else fail "$rt network/isolation probe"; fi
    ctl=$("$RUN" shell --runtime "$rt" --no-netns --workdir "$W/probe_${rt}_control" < "$W/probe.sh" 2>&1)
    echo "--- $rt control with --no-netns:"; echo "$ctl" | grep -E "interfaces|tcp|dns|curl"
    if echo "$ctl" | grep -q "CONNECTED"; then pass "$rt control: the probe does detect a network when one exists"
    else fail "$rt control probe saw no network (the probe would not detect one)"; fi
done

cat > "$W/standin_run.yaml" <<YAML
patient_runs:
  standin_HG00130: {samples: "$STANDIN/standin.tsv"}
threads: {samtools: 4, sv_callers: $CORES}
YAML

hdr "4. patient run offline: Singularity, bundle only"
rm -rf "$W/offline_singularity"; t=$(date +%s)
if "$RUN" patient --runtime singularity --workdir "$W/offline_singularity" --config "$W/standin_run.yaml" \
        --bind "$STANDIN" --cores "$CORES" > "$W/offline_singularity.log" 2>&1; then
    pass "offline Singularity run finished ($(( $(date +%s)-t )) s)"; else fail "offline Singularity run (log $W/offline_singularity.log)"; tail -20 "$W/offline_singularity.log"; fi

hdr "5. patient run offline: Docker, image loaded from the bundle archive"
docker image rm "$CXSV_IMAGE" > /dev/null 2>&1 && echo "removed local $CXSV_IMAGE: run_cxsv.sh must load it from the bundle"
rm -rf "$W/offline_docker"; t=$(date +%s)
if "$RUN" patient --runtime docker --workdir "$W/offline_docker" --config "$W/standin_run.yaml" \
        --bind "$STANDIN" --cores "$CORES" > "$W/offline_docker.log" 2>&1; then
    pass "offline Docker run finished ($(( $(date +%s)-t )) s)"; grep -m1 "loading" "$W/offline_docker.log" || true
else fail "offline Docker run (log $W/offline_docker.log)"; tail -20 "$W/offline_docker.log"; fi

hdr "6. patient run outside the bundle: repository + conda env + lab paths"
source ~/miniforge3/etc/profile.d/conda.sh && conda activate cxsv-pipeline
echo "repo $REPO at $(git -C "$REPO" rev-parse HEAD) (dirty files: $(git -C "$REPO" status --porcelain --untracked-files=no | wc -l)); $(command -v snakemake)"
rm -rf "$W/normal"; mkdir -p "$W/normal"
cat > "$W/normal/site.yaml" <<YAML
reference: "/mnt/sargasso1/CxSV/refs/hs1.fa"
ref_cache: "/mnt/sargasso1/CxSV/refs/ref_cache"
annotation:
  genes_bed:   "/mnt/sargasso1/CxSV/resources/genes.final.bed"
  exons_bed:   "/mnt/sargasso1/CxSV/resources/exons.filtered.bed"
  repeats_bed: "/mnt/sargasso1/CxSV/resources/repeats.final.bed"
  segdup_bed:  "/mnt/sargasso1/CxSV/resources/segdup.final.bed"
panel_release: "$NORMAL_RELEASE"
patient_results_root: "$W/normal/results/patient"
YAML
t=$(date +%s)
if (cd "$REPO" && snakemake -s workflow/patient/Snakefile --cores "$CORES" --keep-going --printshellcmds \
        --configfile "$W/normal/site.yaml" "$W/standin_run.yaml") > "$W/normal.log" 2>&1; then
    pass "normal run finished ($(( $(date +%s)-t )) s)"; else fail "normal run (log $W/normal.log)"; tail -20 "$W/normal.log"; fi

hdr "7. outputs identical (records)"
for x in offline_singularity offline_docker; do
    # same files, different absolute paths: results root and reference inside the container vs on the host
    if python "$REPO/tests/lab/compare_patient_runs.py" "$W/$x/results/patient" "$W/normal/results/patient" \
            --alias-a /work/results/patient=@RESULTS@ --alias-a /bundle/reference=@REF@ --alias-a /work=@WORKDIR@ \
            --alias-a /bundle/resources=@RES@ --alias-a /bundle/catalogues=@CAT@ \
            --alias-b "$W/normal/results/patient=@RESULTS@" --alias-b /mnt/sargasso1/CxSV/refs=@REF@ \
            --alias-b /mnt/sargasso1/CxSV/resources=@RES@ \
            --alias-b /mnt/sargasso1/CxSV/datasets/external_catalogues=@CAT@ \
            --alias-b "$REPO=@WORKDIR@" \
            --out "$W/compare_${x}_vs_normal.tsv" | tail -1; then pass "$x vs normal: no unexpected difference"
    else fail "$x vs normal differ (see $W/compare_${x}_vs_normal.tsv)"; grep -vE "identical|expected" "$W/compare_${x}_vs_normal.tsv"; fi
done
echo "release check record (offline Singularity):"; sed 's/^/  /' "$W/offline_singularity/results/patient/panel_release_check.txt"
echo "fingerprint (offline Singularity): $(grep -E '^(this run|panel|status)' "$W/offline_singularity/results/patient/panel_compatibility.txt" | xargs)"

hdr "8. tampered panel release must be refused"
for kind in byte extra; do
    T="$W/tamper_$kind"; rm -rf "$T"; mkdir -p "$T"
    cp -r "$BUNDLE/panel_release/$CXSV_RELEASE_ID" "$T/tampered_release"
    if [ $kind = byte ]; then
        python - "$T/tampered_release/panel_junction_carriers.tsv" <<'PY'
import sys
p = sys.argv[1]; lines = open(p).read().split("\n"); f = lines[1].split("\t")
f[1] = f[1][:-1] + str((int(f[1][-1]) + 1) % 10); lines[1] = "\t".join(f); open(p, "w").write("\n".join(lines))
print(f"changed one digit of pos1 in row 1: {f[0]}:{f[1]}")
PY
        want="sha256 differs: panel_junction_carriers.tsv"
    else
        printf 'chr1\t1\tchr1\t2\tDEL\t1\n' > "$T/tampered_release/extra_sites.tsv"; want="not in manifest: extra_sites.tsv"
    fi
    printf 'panel_release: "/work/tampered_release"\n' > "$T/tamper.yaml"
    if "$RUN" patient --runtime singularity --workdir "$T" --config "$W/standin_run.yaml" --config "$T/tamper.yaml" \
            --bind "$STANDIN" --cores "$CORES" > "$T/run.log" 2>&1; then
        fail "tampered release ($kind) was NOT refused"
    elif grep -q "$want" "$T/run.log" && [ ! -d "$T/results/patient/samples" ]; then
        pass "tampered release ($kind) refused before any job ran: $(grep -m1 "$want" "$T/run.log" | xargs)"
    else fail "tampered run ($kind) failed, but not with '$want' (log $T/run.log)"; tail -10 "$T/run.log"; fi
done

hdr "summary"
echo "failures: $FAILS"
[ "$FAILS" -eq 0 ]
