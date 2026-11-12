#!/bin/bash
# run the D5 pass-1 tasks of shard I of N (task file lines: key stratum depth), PAR at a time
I=$1; N=$2; PAR=${3:-12}
D5=/mnt/sargasso1/CxSV/calibration/panel_v13_design/d5
awk -v i=$I -v n=$N 'NR%n==i' $D5/seg_tasks.txt | xargs -P $PAR -L 1 bash $D5/tools/d5_seg.sh
echo "shard $I/$N done on $(hostname): $(ls $D5/seg/z6.0/*.json 2>/dev/null | wc -l) z6.0 json total"
