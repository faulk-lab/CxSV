#!/usr/bin/env bash
# Library QC deployment, step 2a: 8-check offline validation of the new bundle on stand-in HG01868.
set -euo pipefail
B=/mnt/sargasso1/CxSV/bundles/cxsv-bundle-20261105-95126e14f036
REL=/mnt/sargasso1/CxSV/build/panel-v1.1-r2/results/panel/release
OUT=/home/minknow/scratch/libqc2/deploy/bundle_validation
echo "tagged bundle MANIFEST.sha256: $(sha256sum /mnt/sargasso1/CxSV/bundles/cxsv-bundle-20261105-8fc5bef725bd/MANIFEST.sha256 | cut -c1-64)"
echo "release manifest: $(sha256sum $REL/release_manifest.tsv | cut -c1-64)"
cd /home/minknow/scratch/libqc2/repo; echo "checkout $(git rev-parse HEAD)"
source ~/miniforge3/etc/profile.d/conda.sh
t0=$(date +%s)
CORES=12 STANDIN_RUN=standin_HG01868 bash tests/lab/validate_bundle.sh $B $OUT /home/minknow/scratch/cxsv_validation/standin $REL && rc=0 || rc=$?
echo "validation exit $rc; wall $(( $(date +%s)-t0 )) s"; exit $rc
