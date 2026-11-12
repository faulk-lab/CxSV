#!/bin/bash
# panel-v1.3 design check, step D2 (run on freya)
set -euo pipefail
DC=/home/minknow/scratch/v13/design_check
SIF=/mnt/sargasso1/CxSV/bundles/cxsv-bundle-20261109-84d380e08cc4/container/cxsv-pipeline.sif
OUT=$DC/d2; mkdir -p $OUT
[ -d $DC/main_tree ] || { mkdir -p $DC/main_tree; git -C $DC/repo archive origin/main | tar -x -C $DC/main_tree; }
git -C $DC/repo rev-parse origin/main > $DC/main_tree/COMMIT_USED
X="singularity exec --containall -B /mnt/sargasso1 -B /home/minknow/scratch $SIF"
$X python -c "import yaml,json;json.dump(yaml.safe_load(open('$DC/repo/config/calling.yaml')),open('$OUT/calling.json','w'))"
$X python $DC/tools/d2_prep.py $OUT
$X python $DC/tools/d2_check.py --old $DC/main_tree --new $DC/repo --design $OUT/design.tsv --bounds $OUT/bounds.tsv \
   --calling-json $OUT/calling.json --out $OUT
mkdir -p /mnt/sargasso1/CxSV/calibration/panel_v13_design/d2 && cp $OUT/*.tsv $OUT/*.json /mnt/sargasso1/CxSV/calibration/panel_v13_design/d2/
