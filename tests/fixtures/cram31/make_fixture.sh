#!/usr/bin/env bash
# Build the committed fixture files from make_fixture.py's output (run in this directory):
#   ref.fa.fai, reads.v31.cram (CRAM 3.1) + reads.v31.cram.crai. reads.sam is not kept.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
python make_fixture.py
samtools faidx ref.fa
samtools view -C -T ref.fa --output-fmt-option version=3.1 --no-PG -o reads.v31.cram reads.sam
# htslib writes the absolute -T path into @SQ UR; keep the relative UR:ref.fa of reads.sam
samtools view -H --no-PG reads.v31.cram | sed 's|\tUR:[^\t]*|\tUR:ref.fa|' > header.sam
samtools reheader --no-PG -i header.sam reads.v31.cram
rm -f header.sam
samtools index reads.v31.cram
! samtools view -H reads.v31.cram | grep -q "UR:/" || { echo "absolute UR left in the header" >&2; exit 1; }
rm -f reads.sam
v=$(head -c 6 reads.v31.cram | od -An -tu1 | awk '{print $5"."$6}')
[ "$v" = 3.1 ] || { echo "expected CRAM 3.1, got $v" >&2; exit 1; }
echo "reads.v31.cram: CRAM $v, $(samtools view -c -T ref.fa reads.v31.cram) records, $(stat -c %s reads.v31.cram) bytes"
