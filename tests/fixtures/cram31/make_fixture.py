#!/usr/bin/env python3
"""Write the CRAM 3.1 callability fixture's reference and reads (stdlib only, seeded).

    python make_fixture.py && bash make_fixture.sh      (in this directory; needs samtools >= 1.12)

ref.fa: two contigs, fx1 30 kb and fx2 20 kb. reads.sam: 2 kb reads that match
the reference exactly, MAPQ 60, coordinate-sorted. fx1 0-15 kb and all of fx2
are tiled every 100 bp (depth ~20), fx1 15-30 kb every 400 bp (depth ~5), so a
min_depth-10 mask covers part of the genome, not all of it.
"""
import hashlib
import random

rng = random.Random(20261104)
contigs = {"fx1": 30000, "fx2": 20000}
seqs = {c: "".join(rng.choice("ACGT") for _ in range(n)) for c, n in contigs.items()}
with open("ref.fa", "w") as fh:
    for c, s in seqs.items():
        fh.write(f">{c}\n" + "".join(s[i:i + 60] + "\n" for i in range(0, len(s), 60)))

L = 2000
with open("reads.sam", "w") as fh:
    fh.write("@HD\tVN:1.6\tSO:coordinate\n")
    for c, n in contigs.items():
        # M5 as htslib computes it; UR relative, so no build path ends up in the CRAM header
        fh.write(f"@SQ\tSN:{c}\tLN:{n}\tM5:{hashlib.md5(seqs[c].encode()).hexdigest()}\tUR:ref.fa\n")
    k = 0
    for c, n in contigs.items():
        pos = 0
        while pos + L <= n:
            k += 1
            flag = 16 if k % 2 else 0
            fh.write(f"r{k}\t{flag}\t{c}\t{pos + 1}\t60\t{L}M\t*\t0\t0\t{seqs[c][pos:pos + L]}\t*\n")
            pos += 400 if (c == "fx1" and pos >= 15000) else 100
