#!/usr/bin/env bash
# =============================================================================
# control.sh — pause, resume, stop or inspect panel build nodes. Run anywhere
# the NAS is mounted (the runners poll their control files every 20 s).
#
#   bash scripts/panel_build/control.sh BUILD pause  NODE|all   drain: start nothing new
#   bash scripts/panel_build/control.sh BUILD resume NODE|all   remove PAUSE
#   bash scripts/panel_build/control.sh BUILD stop   NODE|all   terminate running samples, exit
#   bash scripts/panel_build/control.sh BUILD slots  NODE N     change concurrent samples live
#   bash scripts/panel_build/control.sh BUILD bwlimit NODE MB    throttle staging (0 = off), next copy on
#   bash scripts/panel_build/control.sh BUILD retry  NODE COHORT/SAMPLE...  retry gave_up samples
#   bash scripts/panel_build/control.sh BUILD status            = build_status.py BUILD
#
# After `stop`, `launch.sh BUILD NODE` continues the node (interrupted samples
# are redone from their last complete step).
# =============================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD="$1"; ACT="$2"; shift 2
targets() { if [ "$1" = all ]; then ls "$BUILD/nodes"; else echo "$1"; fi; }
case "$ACT" in
    pause)  for n in $(targets "$1"); do date -Is > "$BUILD/nodes/$n/PAUSE"; echo "$n: PAUSE set"; done ;;
    resume) for n in $(targets "$1"); do rm -f "$BUILD/nodes/$n/PAUSE"; echo "$n: resumed"; done ;;
    stop)   for n in $(targets "$1"); do date -Is > "$BUILD/nodes/$n/STOP"; echo "$n: STOP set"; done ;;
    slots|bwlimit)
        key=slots; [ "$ACT" = bwlimit ] && key=stage_bwlimit_mb_s
        for n in $(targets "$1"); do python3 - "$BUILD/nodes/$n/runner.json" "$key" "$2" <<'PY'
import json, os, sys
p, k, v = sys.argv[1], sys.argv[2], int(sys.argv[3])
c = json.load(open(p)); c[k] = v
json.dump(c, open(p + ".tmp", "w"), indent=2); os.replace(p + ".tmp", p)
print(f"{p}: {k} = {v}")
PY
        done ;;
    retry)  n="$1"; shift; printf '%s\n' "$@" > "$BUILD/nodes/$n/RETRY.tmp"; mv "$BUILD/nodes/$n/RETRY.tmp" "$BUILD/nodes/$n/RETRY"
            echo "$n: RETRY $*" ;;
    status) python3 "$HERE/build_status.py" "$BUILD" "$@" ;;
    *) sed -n '4,16p' "$0"; exit 1 ;;
esac
