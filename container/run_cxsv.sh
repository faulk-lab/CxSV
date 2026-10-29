#!/usr/bin/env bash
# =============================================================================
# run_cxsv.sh — run CxSV from the offline bundle, inside its container, with
# networking disabled. Copied to the bundle root by scripts/make_bundle.sh;
# everything it uses lives in the bundle directory it sits in.
#
#   ./run_cxsv.sh patient --workdir DIR --config patient.yaml [--bind PATH]... [--cores N]
#   ./run_cxsv.sh panel   --workdir DIR --config panel.yaml   [--bind PATH]... [--cores N]
#   ./run_cxsv.sh verify-bundle            sha256 of every bundle file (MANIFEST.sha256)
#   ./run_cxsv.sh verify-release [DIR]     the shipped panel release against its manifest
#   ./run_cxsv.sh netcheck                 show that the container has no network
#   ./run_cxsv.sh shell --workdir DIR      interactive shell in the container
#
#   --workdir DIR     writable run directory (results/, .snakemake/); seen as /work
#   --config FILE     your run config (patient_runs / panel_cohorts, threads, ...);
#                     repeatable, applied after config/bundle_site.yaml
#   --bind PATH       a host directory with inputs (alignments, sample sheets),
#                     mounted READ-ONLY at the same path; repeatable
#   --cores N         snakemake --cores (default 16)
#   --runtime NAME    singularity | apptainer | docker (default: first found,
#                     in that order)
#   --dry-run         resolve the DAG only (snakemake -n)
#   --no-netns        do not create a private network namespace (only if the
#                     site's runtime forbids it; the bundle needs no network)
#   -- ARGS           passed to snakemake
#
# Inside the container: bundle -> /bundle (read-only), workdir -> /work, the
# --bind paths, nothing else (Singularity/Apptainer --containall: no $HOME, no
# host environment; Docker: no host mounts besides these). Network: Singularity /
# Apptainer --net --network none, Docker --network none (loopback only).
# =============================================================================
set -euo pipefail
BUNDLE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "$BUNDLE/bundle.env"          # CXSV_BUNDLE_ID CXSV_COMMIT CXSV_RELEASE_ID CXSV_IMAGE ...

usage() { sed -n '4,30p' "$0" | sed 's/^# \{0,1\}//'; }
[ $# -ge 1 ] || { usage; exit 1; }
CMD="$1"; shift
WORKDIR=""; CONFIGS=(); BINDS=(); CORES=16; RUNTIME="${CXSV_RUNTIME:-}"; DRY=false; NETNS=true; EXTRA=()
while [ $# -gt 0 ]; do
    case "$1" in
        --workdir) WORKDIR="$2"; shift 2 ;;
        --config)  CONFIGS+=("$2"); shift 2 ;;
        --bind)    BINDS+=("$(cd "$2" && pwd)"); shift 2 ;;
        --cores)   CORES="$2"; shift 2 ;;
        --runtime) RUNTIME="$2"; shift 2 ;;
        --dry-run|-n) DRY=true; shift ;;
        --no-netns) NETNS=false; shift ;;
        --) shift; EXTRA=("$@"); break ;;
        -h|--help) usage; exit 0 ;;
        *) if [ "$CMD" = verify-release ] && [ -z "${REL_ARG:-}" ]; then REL_ARG="$1"; shift
           else echo "unknown option: $1"; usage; exit 1; fi ;;
    esac
done

if [ "$CMD" = verify-bundle ]; then
    cd "$BUNDLE" && sha256sum --quiet -c MANIFEST.sha256 && echo "bundle $CXSV_BUNDLE_ID: all $(wc -l < MANIFEST.sha256) files verified"
    exit $?
fi

if [ -z "$RUNTIME" ]; then
    for r in singularity apptainer docker; do command -v "$r" >/dev/null && { RUNTIME=$r; break; }; done
fi
[ -n "$RUNTIME" ] || { echo "no container runtime found (singularity, apptainer or docker)"; exit 1; }

case "$CMD" in
    patient|panel|shell)
        [ -n "$WORKDIR" ] || { echo "--workdir is required"; exit 1; }
        mkdir -p "$WORKDIR/config"; WORKDIR="$(cd "$WORKDIR" && pwd)" ;;
    verify-release|netcheck) WORKDIR="$(mktemp -d)"; trap 'rm -rf "$WORKDIR"' EXIT ;;
    *) echo "unknown command: $CMD"; usage; exit 1 ;;
esac

# run configs are copied into the workdir so the container sees them as /work/config/<name>
CFG_IN=(/bundle/config/bundle_site.yaml)
for f in "${CONFIGS[@]}"; do
    cp "$f" "$WORKDIR/config/$(basename "$f")"; CFG_IN+=("/work/config/$(basename "$f")")
done

ENVS=(CXSV_BUNDLE_ID="$CXSV_BUNDLE_ID" CXSV_COMMIT="$CXSV_COMMIT" CXSV_IMAGE_ID="$CXSV_IMAGE_ID"
      CXSV_SIF_SHA256="$CXSV_SIF_SHA256" CXSV_RUNTIME="$RUNTIME")

case "$CMD" in
    patient|panel)
        INNER=(snakemake -s "/bundle/repo/workflow/$CMD/Snakefile" --directory /work --cores "$CORES"
               --keep-going --printshellcmds --configfile "${CFG_IN[@]}")
        $DRY && INNER+=(--dry-run)
        INNER+=("${EXTRA[@]}") ;;
    verify-release)
        INNER=(python /bundle/repo/workflow/scripts/verify_panel_release.py "${REL_ARG:-/bundle/panel_release/$CXSV_RELEASE_ID}") ;;
    netcheck)
        INNER=(bash -c 'echo "interfaces:"; tail -n +3 /proc/net/dev | cut -d: -f1 | sed "s/^ */  /";
               if python -c "import socket; socket.create_connection((\"1.1.1.1\", 443), timeout=5)" 2>/dev/null;
               then echo "NETWORK: ON (reached 1.1.1.1:443)"; exit 1;
               else echo "NETWORK: OFF (no route to 1.1.1.1:443)"; fi') ;;
    shell) INNER=(bash) ;;
esac

if [ "$RUNTIME" = docker ]; then
    if ! docker image inspect "$CXSV_IMAGE" >/dev/null 2>&1; then
        echo "loading $CXSV_IMAGE from the bundle ..." >&2
        docker load -i "$BUNDLE/container/cxsv-pipeline.docker.tar.gz" >&2
    fi
    got=$(docker image inspect -f '{{.Id}}' "$CXSV_IMAGE")
    [ "$got" = "$CXSV_IMAGE_ID" ] || { echo "image $CXSV_IMAGE is $got, bundle expects $CXSV_IMAGE_ID"; exit 1; }
    D=(docker run --rm -i --user "$(id -u):$(id -g)" -e HOME=/tmp -w /work
       -v "$BUNDLE:/bundle:ro" -v "$WORKDIR:/work")
    $NETNS && D+=(--network none)
    [ "$CMD" = shell ] && [ -t 0 ] && D+=(-t)
    for b in "${BINDS[@]}"; do D+=(-v "$b:$b:ro"); done
    for e in "${ENVS[@]}"; do D+=(-e "$e"); done
    exec "${D[@]}" "$CXSV_IMAGE" "${INNER[@]}"
else
    S=("$RUNTIME" exec --containall --pwd /work --bind "$BUNDLE:/bundle:ro" --bind "$WORKDIR:/work")
    $NETNS && S+=(--net --network none)
    for b in "${BINDS[@]}"; do S+=(--bind "$b:$b:ro"); done
    for e in "${ENVS[@]}"; do S+=(--env "$e"); done
    exec "${S[@]}" "$BUNDLE/container/cxsv-pipeline.sif" "${INNER[@]}"
fi
