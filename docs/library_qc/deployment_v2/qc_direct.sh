#!/usr/bin/env bash
# Library QC follow-ups, step 2c: the bundle's library_qc.py (inside the bundle SIF, network off) on
# genomes without a full patient run here. Arguments: BUNDLE LABEL ALN MERGED_VCF PLATFORM ...
# The params JSON is what rule library_qc_params writes: the bundle's library_qc config block plus
# group_override '' (no sample-sheet qc_group) and the platform.
set -euo pipefail
B=$1; shift; O=/home/minknow/scratch/libqc3/deploy/qc_direct; mkdir -p $O
SIF=$B/container/cxsv-pipeline.sif; RT=$(command -v singularity || command -v apptainer)
while [ $# -ge 4 ]; do
  L=$1; A=$2; V=$3; P=$4; shift 4; W=$O/$L; mkdir -p $W
  "$RT" exec --containall --net --network none --bind "$B:/bundle:ro" --bind "$W:$W" \
    --bind "$(dirname "$A"):$(dirname "$A"):ro" --bind "$(dirname "$V"):$(dirname "$V"):ro" "$SIF" bash -c "
    export REF_PATH=/bundle/reference/ref_cache/%2s/%2s/%s REF_CACHE=/bundle/reference/ref_cache/%2s/%2s/%s
    python -c \"import json, yaml; c = yaml.safe_load(open('/bundle/repo/config/library_qc.yaml'))['library_qc']; print(json.dumps({**c, 'group_override': '', 'platform': '$P'}, sort_keys=True))\" > $W/library_qc_params.json
    python /bundle/repo/workflow/scripts/library_qc.py --aln '$A' --ref /bundle/reference/hs1.fa \
      --windows /bundle/repo/config/library_qc_windows.hs1.bed --merged-vcf '$V' --sample $L --min-mapq 20 \
      --thresholds $W/library_qc_params.json --group-override= --platform=$P --out $W/library_qc.tsv 2> $W/library_qc.log"
  python3 -c "
import csv; r = next(csv.DictReader(open('$W/library_qc.tsv'), delimiter='\t'))
print('$L', *(f'{k}={r.get(k)}' for k in ('qc_group', 'qc_group_source', 'qc_flag', 'z_svim_only_bnd_per_gb', 'z_svim_only_any_per_gb', 'reads_median_mean_q', 'qc_group_evidence')), sep=' | ')"
done
