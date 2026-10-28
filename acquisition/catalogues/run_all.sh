#!/usr/bin/env bash
# Download (verified) and convert every external SV catalogue end to end.
#
#   conda activate cxsv-catalogues          # envs/cxsv-catalogues.yaml
#   bash acquisition/catalogues/run_all.sh [ROOT] [FAI] [THREADS]
#
# Defaults are the lab NAS paths. Re-running is safe: files whose published
# checksum already matches are kept, and outputs are replaced atomically.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${1:-/mnt/sargasso1/CxSV/datasets/external_catalogues}"
FAI="${2:-/mnt/sargasso1/CxSV/refs/hs1.fa.fai}"
THREADS="${3:-6}"
CODE_VERSION="$(git -C "$HERE" rev-parse --short HEAD 2>/dev/null || echo unknown)"
if ! git -C "$HERE" diff --quiet -- . 2>/dev/null; then CODE_VERSION="${CODE_VERSION}-dirty"; fi

for t in bcftools bgzip tabix liftOver curl; do
  command -v "$t" >/dev/null || { echo "missing tool: $t (activate cxsv-catalogues)" >&2; exit 1; }
done

echo "== free space on $(df -h "$ROOT" | tail -1)"
python3 "$HERE/download_catalogues.py" --root "$ROOT" --dry-run
python3 "$HERE/download_catalogues.py" --root "$ROOT"

for cat in CoLoRSdb HGSVC3 HPRC gnomAD-SV; do
  python3 "$HERE/convert_catalogue.py" --catalogue "$cat" --root "$ROOT" --fai "$FAI" \
      --threads "$THREADS" --code-version "$CODE_VERSION"
done
python3 "$HERE/make_summary.py" --root "$ROOT"
