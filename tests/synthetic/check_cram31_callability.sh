#!/usr/bin/env bash
# =============================================================================
# Callability masks from CRAM 3.1 (workflow/scripts/callable_mask.sh).
#
#   bash tests/synthetic/check_cram31_callability.sh [WORKDIR]
#
# Fixture: tests/fixtures/cram31/ (reads.v31.cram, CRAM 3.1, 50 kb reference;
# make_fixture.sh). Checks:
#   1. BAM made from the fixture -> a mask with > 0 callable bp (reference)
#   2. CRAM 3.1, cram31=transcode -> the same mask, byte for byte; the input
#      CRAM unchanged; no temporary CRAM left behind
#   3. CRAM 3.1, cram31=fail -> exit != 0, message names CRAM 3.1, no mask
#   4. min_depth above every depth -> exit != 0 ("0 callable bp"), no mask
#   5. (information) mosdepth run directly on the CRAM 3.1: callable bp it
#      reports (0 with the pinned mosdepth 0.3.10, which cannot decode 3.1)
# Needs samtools and mosdepth (the cxsv-pipeline env). Exit 0 = 1-4 passed.
# =============================================================================
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
W="${1:-$(mktemp -d)}"; mkdir -p "$W"; W="$(cd "$W" && pwd)"
T="$W/cram31"; rm -rf "$T"; mkdir -p "$T"
FX="$REPO/tests/fixtures/cram31"
cp "$FX/ref.fa" "$FX/ref.fa.fai" "$FX/reads.v31.cram" "$FX/reads.v31.cram.crai" "$T/"
CM="$REPO/workflow/scripts/callable_mask.sh"
cd "$T"
echo "mosdepth $(mosdepth --version 2>&1 | head -1); $(samtools --version | head -1)"
v=$(head -c 6 reads.v31.cram | od -An -tu1 | awk '{print $5"."$6}')
[ "$v" = 3.1 ] || { echo "FAIL fixture is CRAM $v, not 3.1"; exit 1; }
sum_before=$(sha256sum reads.v31.cram | cut -c1-64)
bp_of() { zcat "$1" | awk -F'\t' '{s += $3 - $2} END {printf "%d\n", s}'; }

# 1. BAM reference mask
samtools view -b -T ref.fa -o reads.bam reads.v31.cram && samtools index reads.bam
bash "$CM" reads.bam ref.fa bam/mosdepth bam.callable.bed.gz 10 20 2 > bam.log 2>&1 \
    || { echo "FAIL mask from BAM"; cat bam.log; exit 1; }
bp_bam=$(bp_of bam.callable.bed.gz)
[ "$bp_bam" -gt 0 ] && [ "$bp_bam" -lt 50000 ] || { echo "FAIL BAM mask has $bp_bam bp"; exit 1; }
echo "PASS 1 BAM: $bp_bam callable bp of 50000"

# 2. CRAM 3.1, transcode
bash "$CM" reads.v31.cram ref.fa c31/mosdepth c31.callable.bed.gz 10 20 2 transcode > c31.log 2>&1 \
    || { echo "FAIL mask from CRAM 3.1 (transcode)"; cat c31.log; exit 1; }
grep -q "transcoded CRAM 3.1 -> 3.0" c31.log || { echo "FAIL no transcode recorded"; cat c31.log; exit 1; }
cmp -s <(zcat bam.callable.bed.gz) <(zcat c31.callable.bed.gz) || { echo "FAIL CRAM 3.1 mask differs from BAM mask"; exit 1; }
[ "$(sha256sum reads.v31.cram | cut -c1-64)" = "$sum_before" ] || { echo "FAIL input CRAM was modified"; exit 1; }
[ -z "$(ls c31/ | grep -E 'v30\.cram' || true)" ] || { echo "FAIL temporary CRAM left: $(ls c31/)"; exit 1; }
echo "PASS 2 CRAM 3.1 transcode: identical mask ($(bp_of c31.callable.bed.gz) bp), input unchanged, no temporary left"

# 3. CRAM 3.1, fail
if bash "$CM" reads.v31.cram ref.fa f31/mosdepth f31.callable.bed.gz 10 20 2 fail > f31.log 2>&1; then
    echo "FAIL cram31=fail accepted a CRAM 3.1"; exit 1
fi
grep -q "is CRAM 3.1" f31.log && [ ! -e f31.callable.bed.gz ] || { echo "FAIL cram31=fail message/mask"; cat f31.log; exit 1; }
echo "PASS 3 CRAM 3.1 fail: refused ($(grep -m1 -o 'is CRAM 3.1[^;]*' f31.log))"

# 4. 0 callable bp
if bash "$CM" reads.bam ref.fa z/mosdepth z.callable.bed.gz 1000 20 2 > z.log 2>&1; then
    echo "FAIL a 0-bp mask was accepted"; exit 1
fi
grep -q "0 callable bp" z.log && [ ! -e z.callable.bed.gz ] || { echo "FAIL 0-bp message/mask"; cat z.log; exit 1; }
echo "PASS 4 empty mask refused"

# 5. information: mosdepth directly on CRAM 3.1
mkdir -p raw
MOSDEPTH_Q0=LOW MOSDEPTH_Q1=CALLABLE mosdepth -t 2 -n -x --mapq 20 --quantize 0:10: -f ref.fa raw/m reads.v31.cram > raw.log 2>&1 || true
raw_bp=$( (zcat raw/m.quantized.bed.gz 2>/dev/null || true) | awk -F'\t' '$4 == "CALLABLE" {s += $3 - $2} END {printf "%d\n", s}')
echo "INFO mosdepth directly on CRAM 3.1: $raw_bp callable bp (BAM: $bp_bam); decode errors logged: $(grep -ciE 'decode failure|\[E::' raw.log || true)"
