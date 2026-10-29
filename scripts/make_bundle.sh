#!/usr/bin/env bash
# =============================================================================
# make_bundle.sh — build the offline bundle for the restricted environment
# (roadmap item 2; docs/RELEASE_AND_BUNDLE.md), reproducibly from a commit.
#
#   scripts/make_bundle.sh image  --commit REV [--out ROOT]
#       build the container image from the commit's container/Dockerfile and
#       envs/cxsv-pipeline.explicit.txt, save it (Docker archive) and convert
#       it to a SIF, under ROOT/_images/<lock sha12>/ (image.env records the
#       image id and SIF sha256). Reused when it already exists: the image
#       holds the environment only, so it is keyed by the lock, not the commit.
#
#   scripts/make_bundle.sh bundle --commit REV --release DIR [--out ROOT]
#           [--bundle-id ID] [--refs DIR] [--resources DIR] [--catalogues DIR]
#       assemble ROOT/<bundle id>/: repo at the commit (git archive + git
#       bundle), the image (SIF + Docker archive), hs1 + .fai + MD5 cache,
#       annotation resources, the panel release DIR (must verify), the external
#       catalogue junction-site tables with their manifests (each verified),
#       run_cxsv.sh, config/bundle_site.yaml, README.md, bundle.json and a
#       checksum manifest over everything (MANIFEST.sha256, sha256sum -c format;
#       bundle_manifest.tsv with sizes and components).
#
# Defaults: ROOT /mnt/sargasso1/CxSV/bundles, refs /mnt/sargasso1/CxSV/refs,
# resources /mnt/sargasso1/CxSV/resources, catalogues
# /mnt/sargasso1/CxSV/datasets/external_catalogues. Bundle id default:
# cxsv-bundle-<YYYYMMDD>-<commit12>. Needs git, docker, singularity|apptainer,
# rsync, sha256sum. Run it from a clone of the repository. Never pushes
# anything and creates no tags (a temporary local ref for `git bundle` is
# deleted again).
# =============================================================================
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MODE="${1:-}"; shift || true
COMMIT_REV=""; RELEASE=""; OUT=/mnt/sargasso1/CxSV/bundles; BUNDLE_ID=""
REFS=/mnt/sargasso1/CxSV/refs; RESOURCES=/mnt/sargasso1/CxSV/resources
CATALOGUES=/mnt/sargasso1/CxSV/datasets/external_catalogues
while [ $# -gt 0 ]; do
    case "$1" in
        --commit) COMMIT_REV="$2"; shift 2 ;;   --release) RELEASE="$2"; shift 2 ;;
        --out) OUT="$2"; shift 2 ;;             --bundle-id) BUNDLE_ID="$2"; shift 2 ;;
        --refs) REFS="$2"; shift 2 ;;           --resources) RESOURCES="$2"; shift 2 ;;
        --catalogues) CATALOGUES="$2"; shift 2 ;;
        *) echo "unknown option $1"; sed -n '4,30p' "$0"; exit 1 ;;
    esac
done
[[ "$MODE" == image || "$MODE" == bundle ]] || { sed -n '4,30p' "$0"; exit 1; }
[ -n "$COMMIT_REV" ] || { echo "--commit is required"; exit 1; }
COMMIT=$(git -C "$REPO" rev-parse --verify "${COMMIT_REV}^{commit}")
C12=${COMMIT:0:12}
SINGULARITY=$(command -v singularity || command -v apptainer)
T_START=$(date +%s)
log() { echo "[$(date +%H:%M:%S)] $*"; }

# exported tree of the commit (code and Dockerfile come from the commit, never the worktree)
EXPORT=$(mktemp -d "${TMPDIR:-/tmp}/cxsv_export.XXXXXX"); trap 'rm -rf "$EXPORT"' EXIT
git -C "$REPO" archive --format=tar "$COMMIT" | tar -x -C "$EXPORT"
LOCK_SHA=$(sha256sum "$EXPORT/envs/cxsv-pipeline.explicit.txt" | cut -c1-64)
IMG_DIR="$OUT/_images/lock-${LOCK_SHA:0:12}"
IMAGE_TAG="cxsv-pipeline:lock-${LOCK_SHA:0:12}"

build_image() {
    if [ -f "$IMG_DIR/image.env" ]; then log "image: reusing $IMG_DIR"; return 0; fi
    mkdir -p "$IMG_DIR.partial"
    local t0; t0=$(date +%s)
    log "image: docker build $IMAGE_TAG from commit $C12"
    docker build -f "$EXPORT/container/Dockerfile" --build-arg ENV_LOCK_SHA256="$LOCK_SHA" \
        --build-arg IMAGE_TAG="$IMAGE_TAG" -t "$IMAGE_TAG" "$EXPORT" > "$IMG_DIR.partial/docker_build.log" 2>&1
    local id; id=$(docker image inspect -f '{{.Id}}' "$IMAGE_TAG")
    log "image: $id; saving and converting to SIF"
    docker save -o "$IMG_DIR.partial/cxsv-pipeline.docker.tar" "$IMAGE_TAG"
    "$SINGULARITY" build "$IMG_DIR.partial/cxsv-pipeline.sif" \
        "docker-archive://$IMG_DIR.partial/cxsv-pipeline.docker.tar" > "$IMG_DIR.partial/sif_build.log" 2>&1
    gzip -n "$IMG_DIR.partial/cxsv-pipeline.docker.tar"
    cp "$EXPORT/container/Dockerfile" "$EXPORT/envs/cxsv-pipeline.explicit.txt" "$IMG_DIR.partial/"
    cat > "$IMG_DIR.partial/image.env" <<EOF
CXSV_IMAGE=$IMAGE_TAG
CXSV_IMAGE_ID=$id
CXSV_SIF_SHA256=$(sha256sum "$IMG_DIR.partial/cxsv-pipeline.sif" | cut -c1-64)
CXSV_DOCKER_ARCHIVE_SHA256=$(sha256sum "$IMG_DIR.partial/cxsv-pipeline.docker.tar.gz" | cut -c1-64)
CXSV_ENV_LOCK_SHA256=$LOCK_SHA
CXSV_IMAGE_BASE=$(grep -m1 '^FROM ' "$EXPORT/container/Dockerfile" | cut -d' ' -f2)
CXSV_IMAGE_BUILT_FROM_COMMIT=$COMMIT
CXSV_IMAGE_BUILD_SECONDS=$(( $(date +%s) - t0 ))
EOF
    mv "$IMG_DIR.partial" "$IMG_DIR"
    log "image: done ($(( $(date +%s) - t0 )) s)"
}

build_image
[ "$MODE" = image ] && { cat "$IMG_DIR/image.env"; exit 0; }

# ─── bundle ──────────────────────────────────────────────────────────────────
[ -n "$RELEASE" ] && [ -f "$RELEASE/release.json" ] || { echo "--release DIR with release.json is required"; exit 1; }
# shellcheck disable=SC1091
source "$IMG_DIR/image.env"
in_image() { "$SINGULARITY" exec --containall --bind "$EXPORT:/repo:ro" --bind "$RELEASE:/release:ro" \
                 "$IMG_DIR/cxsv-pipeline.sif" "$@"; }
log "release: verifying $RELEASE"
in_image python /repo/workflow/scripts/verify_panel_release.py /release
RELEASE_ID=$(in_image python -c "import json; print(json.load(open('/release/release.json'))['release_id'])")
TEST_RELEASE=$(in_image python -c "import json; print(str(json.load(open('/release/release.json'))['test_release']).lower())")
BUNDLE_ID="${BUNDLE_ID:-cxsv-bundle-$(date +%Y%m%d)-$C12}"
DEST="$OUT/$BUNDLE_ID"; STAGE="$DEST.partial"
[ ! -e "$DEST" ] || { echo "$DEST exists; choose another --bundle-id"; exit 1; }
rm -rf "$STAGE"; mkdir -p "$STAGE"/{source,container,reference,resources,panel_release,catalogues,config}
declare -A SECS
step() { SECS[$1]=$(( $(date +%s) - $2 )); }

t=$(date +%s); log "repo: git archive + git bundle of $COMMIT"
cp -a "$EXPORT" "$STAGE/repo"
git -C "$REPO" archive --format=tar.gz -o "$STAGE/source/cxsv-$C12.tar.gz" "$COMMIT"
TMPREF="refs/cxsv-bundle-tmp/$C12"
git -C "$REPO" update-ref "$TMPREF" "$COMMIT"
git -C "$REPO" bundle create "$STAGE/source/cxsv-$C12.bundle" "$TMPREF" 2> /dev/null
git -C "$REPO" update-ref -d "$TMPREF"
git -C "$REPO" bundle verify "$STAGE/source/cxsv-$C12.bundle" > /dev/null 2>&1 && log "repo: git bundle verifies"
step repo "$t"

t=$(date +%s); log "container: SIF + Docker archive"
cp "$IMG_DIR"/{cxsv-pipeline.sif,cxsv-pipeline.docker.tar.gz,image.env,Dockerfile,cxsv-pipeline.explicit.txt} "$STAGE/container/"
step container "$t"

t=$(date +%s); log "reference: hs1 + .fai + MD5 cache (mtimes kept: the cache must stay newer than hs1.fa)"
rsync -a "$REFS/hs1.fa" "$REFS/hs1.fa.fai" "$STAGE/reference/"
rsync -a "$REFS/ref_cache/" "$STAGE/reference/ref_cache/"
step reference "$t"

t=$(date +%s); log "resources"
rsync -a "$RESOURCES"/*.bed "$STAGE/resources/"
step resources "$t"

t=$(date +%s); log "panel release $RELEASE_ID"
rsync -a "$RELEASE/" "$STAGE/panel_release/$RELEASE_ID/"
step release "$t"

t=$(date +%s); log "catalogues: junction-site tables + manifests"
cp "$CATALOGUES/external_catalogues_summary.tsv" "$CATALOGUES/README.md" "$STAGE/catalogues/" 2>/dev/null || true
CAT_LIST=()
for man in "$CATALOGUES"/*/*/*.manifest.tsv; do
    d=$(dirname "$man"); rel=${d#"$CATALOGUES"/}
    out=$(awk -F'\t' '$1=="output_file"{print $2}' "$man"); want=$(awk -F'\t' '$1=="output_sha256"{print $2}' "$man")
    got=$(sha256sum "$d/$out" | cut -c1-64)
    [ "$got" = "$want" ] || { echo "catalogue $rel: $out sha256 $got != manifest $want"; exit 1; }
    mkdir -p "$STAGE/catalogues/$rel"
    for f in "$out" "$out.tbi" "$(basename "$man")" SOURCES.tsv; do cp -p "$d/$f" "$STAGE/catalogues/$rel/"; done
    for f in "$d"/README* "$d"/*liftover_stats.tsv; do [ -f "$f" ] && cp -p "$f" "$STAGE/catalogues/$rel/"; done
    CAT_LIST+=("$rel/$out")
    log "  $rel/$out verified against its manifest"
done
step catalogues "$t"

# ─── run entry point, site config, metadata, README ─────────────────────────
cp "$EXPORT/container/run_cxsv.sh" "$STAGE/run_cxsv.sh"; chmod +x "$STAGE/run_cxsv.sh"
cat > "$STAGE/config/bundle_site.yaml" <<EOF
# Bundle site config: every path is inside the container (bundle at /bundle,
# run directory at /work). Applied first by run_cxsv.sh; put patient_runs /
# threads in your own --config file.
reference: "/bundle/reference/hs1.fa"
ref_cache: "/bundle/reference/ref_cache"
annotation:
  genes_bed:   "/bundle/resources/genes.final.bed"
  exons_bed:   "/bundle/resources/exons.filtered.bed"
  repeats_bed: "/bundle/resources/repeats.final.bed"
  segdup_bed:  "/bundle/resources/segdup.final.bed"
panel_release: "/bundle/panel_release/$RELEASE_ID"
patient_results_root: "/work/results/patient"
panel_results_root:   "/work/results/panel"
EOF
cat > "$STAGE/bundle.env" <<EOF
CXSV_BUNDLE_ID=$BUNDLE_ID
CXSV_COMMIT=$COMMIT
CXSV_RELEASE_ID=$RELEASE_ID
CXSV_TEST_RELEASE=$TEST_RELEASE
$(grep -E '^CXSV_(IMAGE|IMAGE_ID|SIF_SHA256|DOCKER_ARCHIVE_SHA256|ENV_LOCK_SHA256)=' "$IMG_DIR/image.env")
EOF
sed -e "s|@BUNDLE_ID@|$BUNDLE_ID|g" -e "s|@COMMIT@|$COMMIT|g" -e "s|@RELEASE_ID@|$RELEASE_ID|g" \
    -e "s|@TEST_RELEASE@|$TEST_RELEASE|g" -e "s|@IMAGE@|$CXSV_IMAGE|g" -e "s|@IMAGE_ID@|$CXSV_IMAGE_ID|g" \
    -e "s|@SIF_SHA256@|$CXSV_SIF_SHA256|g" "$EXPORT/container/BUNDLE_README.md" > "$STAGE/README.md"

TOTAL_SECS=$(( $(date +%s) - T_START ))
# bundle.json first (sizes from a file walk), so MANIFEST.sha256 covers it too
python3 - "$STAGE" "$BUNDLE_ID" "$COMMIT" "$RELEASE_ID" "$TEST_RELEASE" "$TOTAL_SECS" \
    "repo=${SECS[repo]} container=${SECS[container]} reference=${SECS[reference]} resources=${SECS[resources]} release=${SECS[release]} catalogues=${SECS[catalogues]}" \
    "${CAT_LIST[@]}" <<'PY'
import json, os, sys, datetime
stage, bid, commit, rid, test, total, secs, *cats = sys.argv[1:]
env = dict(l.strip().split("=", 1) for l in open(os.path.join(stage, "container", "image.env")) if "=" in l)
comp_bytes, n = {}, 0
for root, _, files in os.walk(stage):
    for f in files:
        path = os.path.relpath(os.path.join(root, f), stage)
        comp = path.split("/")[0] if "/" in path else "top"
        comp_bytes[comp] = comp_bytes.get(comp, 0) + os.path.getsize(os.path.join(root, f))
        n += 1
rel = json.load(open(os.path.join(stage, "panel_release", rid, "release.json")))
meta = {
    "schema": "cxsv-offline-bundle/1", "bundle_id": bid,
    "created_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    "source": {"commit": commit, "repo": "https://github.com/faulk-lab/CxSV", "tree": "repo/ (git archive)",
               "archive": f"source/cxsv-{commit[:12]}.tar.gz", "git_bundle": f"source/cxsv-{commit[:12]}.bundle"},
    "container": {"image": env["CXSV_IMAGE"], "image_id": env["CXSV_IMAGE_ID"], "sif": "container/cxsv-pipeline.sif",
                  "sif_sha256": env["CXSV_SIF_SHA256"], "docker_archive": "container/cxsv-pipeline.docker.tar.gz",
                  "docker_archive_sha256": env["CXSV_DOCKER_ARCHIVE_SHA256"],
                  "env_lock_sha256": env["CXSV_ENV_LOCK_SHA256"], "base_image": env["CXSV_IMAGE_BASE"],
                  "image_built_from_commit": env.get("CXSV_IMAGE_BUILT_FROM_COMMIT")},
    "panel_release": {"release_id": rid, "test_release": test == "true", "path": f"panel_release/{rid}",
                      "fingerprint": rel["calling_fingerprint"]["sha256"], "source_commit": rel["source"]["commit"],
                      "built_in_image": rel["runtime"].get("image_id"), "built_in_sif": rel["runtime"].get("sif_sha256")},
    "reference": {"fasta": "reference/hs1.fa", "fai": "reference/hs1.fa.fai", "md5_cache": "reference/ref_cache"},
    "resources": "resources/", "catalogues": [f"catalogues/{c}" for c in cats],
    "integrity": "MANIFEST.sha256 (sha256sum -c format) over every file but itself and bundle_manifest.tsv "
                 "(the same list with sizes and components)",
    "size_bytes_before_manifests": sum(comp_bytes.values()), "n_files": n + 1,
    "size_bytes_by_component": dict(sorted(comp_bytes.items())),
    "build_seconds_before_checksums": {"total": int(total), **{k: int(v) for k, v in (x.split("=") for x in secs.split())}},
}
json.dump(meta, open(os.path.join(stage, "bundle.json"), "w"), indent=2)
print(json.dumps({k: meta[k] for k in ("bundle_id", "size_bytes_before_manifests", "n_files")}))
PY

t=$(date +%s); log "checksums"
( cd "$STAGE" && find . -type f ! -name MANIFEST.sha256 ! -name bundle_manifest.tsv -printf '%P\n' \
    | LC_ALL=C sort | xargs -d '\n' sha256sum > MANIFEST.sha256.tmp && mv MANIFEST.sha256.tmp MANIFEST.sha256 )
log "checksums: $(( $(date +%s) - t )) s"
# sizes are stat-ed relative to the bundle root
( cd "$STAGE" && awk 'BEGIN{OFS="\t"; print "path","bytes","sha256","component"}
     {sha=$1; p=substr($0, 67); c=(index(p,"/") ? substr(p,1,index(p,"/")-1) : "top");
      cmd="stat -c %s \"" p "\""; cmd | getline sz; close(cmd); print p, sz, sha, c}' MANIFEST.sha256 > bundle_manifest.tsv )
mv "$STAGE" "$DEST"
log "bundle: $DEST ($(du -sb "$DEST" | cut -f1) bytes, $(du -sh "$DEST" | cut -f1)); total $(( $(date +%s) - T_START )) s"
log "verify with: $DEST/run_cxsv.sh verify-bundle"
