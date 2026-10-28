#!/usr/bin/env bash
# betsy_stress_test.sh - post-BIOS-update acceptance test for betsy (CxSV project)
#
# Why: betsy hard-reset 43 times (2026-10-21..26) under 24-thread minimap2 load, each
# reset leaving a fatal machine check on logical CPU 0/16 (physical core 0). After the
# BIOS update + optimized defaults, this test reproduces that exact workload and checks
# both stability (no reset) and correctness (bit-identical output vs freya).
#
# Phases (default ~6.2 h):
#   1. minimap2 loop, MM_HOURS h: production command (-t 24 -x map-ont -a --eqx -K 10g),
#      SAM (minus @PG) md5 compared to reference_sam.md5 computed on freya.
#   2. hwprobe cpu CPU_SECS: one pinned thread per logical CPU, per-core consistency.
#   3. hwprobe mem MEM_GB x2 (fresh mapping each): DRAM verify.
# Logs live on betsy local disk (~/stress/run_*) so they survive a reset; the run dir
# is copied to the NAS scratch dir at the end.
set +e
MM_HOURS="${MM_HOURS:-4.5}"; CPU_SECS="${CPU_SECS:-3600}"; MEM_GB="${MEM_GB:-96}"; MEM_SECS="${MEM_SECS:-1200}"
S=/mnt/sargasso1/CxSV/scratch/betsy_stress
IDX=~/refs/hs1/hs1.map-ont.k17.mmi
PROBE=~/stress/hwprobe
W=~/stress/run_$(date +%Y%m%d_%H%M); mkdir -p "$W"; LOG="$W/progress.log"
log(){ echo "$(date -Is) $*" | tee -a "$LOG"; }

source ~/miniforge3/etc/profile.d/conda.sh; conda activate cxsvalign
BOOT0=$(cat /proc/sys/kernel/random/boot_id)
log "START boot_id=$BOOT0 bios=$(cat /sys/devices/virtual/dmi/id/bios_version) microcode=$(grep -m1 microcode /proc/cpuinfo | awk '{print $3}') minimap2=$(minimap2 --version)"
[ "$(minimap2 --version)" = "2.24-r1122" ] || { log "ABORT minimap2 version mismatch"; exit 4; }
REF=$(cat "$S/reference_sam.md5"); log "reference md5 (freya) = $REF"
cp "$S/test_reads.fq" "$W/test_reads.fq"
[ "$(md5sum < "$W/test_reads.fq" | awk '{print $1}')" = "$(awk '{print $1}' "$S/test_reads.fq.md5")" ] || { log "ABORT test reads copy md5 mismatch"; exit 5; }

# temperature logger (k10temp Tctl + all CCD sensors), every 60 s
HW=$(for h in /sys/class/hwmon/hwmon*; do [ "$(cat $h/name)" = k10temp ] && echo $h; done | head -1)
( while true; do echo "$(date +%H:%M:%S) $(for t in $HW/temp*_input; do printf '%s=%d ' $(cat ${t%_input}_label 2>/dev/null || basename $t) $(( $(cat $t)/1000 )); done)"; sleep 60; done ) >> "$W/temps.log" &
TPID=$!

# ---- phase 1: production minimap2 workload ----
END=$(( SECONDS + $(awk -v h="$MM_HOURS" 'BEGIN{printf "%d", h*3600}') ))
i=0; ok=0; bad=0; fail=0
log "PHASE1 minimap2 loop for ${MM_HOURS} h"
while [ $SECONDS -lt $END ]; do
  i=$((i+1)); t=$SECONDS
  h=$(minimap2 -t 24 -x map-ont -a --eqx -K 10g "$IDX" "$W/test_reads.fq" 2>>"$W/mm2.err" | grep -v '^@PG' | md5sum | awk '{print $1}'; exit ${PIPESTATUS[0]})
  rc=$?
  if [ $rc -ne 0 ]; then fail=$((fail+1)); m=FAIL; elif [ "$h" = "$REF" ]; then ok=$((ok+1)); m=match; else bad=$((bad+1)); m=MISMATCH; fi
  log "mm2 iter=$i rc=$rc secs=$((SECONDS-t)) md5=$h $m"
done
log "PHASE1 done iters=$i match=$ok mismatch=$bad failed=$fail"

# ---- phase 2: per-core CPU consistency ----
log "PHASE2 hwprobe cpu ${CPU_SECS}s"
$PROBE cpu "$CPU_SECS" 2>&1 | tee "$W/hwprobe_cpu.txt" | grep -E 'SUMMARY|RESULT|TOTAL' | while read -r l; do log "cpu: $l"; done

# ---- phase 3: DRAM verify, two independent mappings ----
for k in 1 2; do
  log "PHASE3 hwprobe mem ${MEM_GB}GB pass $k"
  $PROBE mem "$MEM_GB" 16 "$MEM_SECS" 2>&1 | tee "$W/hwprobe_mem$k.txt" | grep -E 'RESULT|TOTAL' | while read -r l; do log "mem$k: $l"; done
done

kill $TPID 2>/dev/null
BOOT1=$(cat /proc/sys/kernel/random/boot_id)
MCE=$(journalctl -b 0 -k --no-pager 2>/dev/null | grep -cE 'Machine Check|mce: \[Hardware Error\]')
log "END boot_id_unchanged=$([ "$BOOT0" = "$BOOT1" ] && echo yes || echo NO) mce_lines_this_boot=$MCE max_tctl=$(awk '{for(i=2;i<=NF;i++){split($i,a,"="); if(a[2]>m)m=a[2]}}END{print m}' "$W/temps.log")C"
rm -f "$W/test_reads.fq"
mkdir -p "$S/results"; cp -r "$W" "$S/results/"
