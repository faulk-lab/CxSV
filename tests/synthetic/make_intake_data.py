#!/usr/bin/env python3
"""Synthetic inputs for the patient intake test (tests/synthetic/check_intake.sh).

    python tests/synthetic/make_intake_data.py OUTDIR [--seed 7]

Writes, from one seeded two-contig reference (chrA 150 kb, chrB 100 kb, hs1-style names):
  ref/ref.fa(.fai)                 the intake reference
  ref/other.fa(.fai)               the same sequence under other names (A, B), for inputs
                                   "from another aligner / reference naming"
  reads/<name>.mod.bam                unaligned BAMs, dorado style: MM/ML/MN/qs plus a per-read
                                   RG:Z tag (which intake must NOT copy)
     modbam   400 reads (12 unmappable random reads, 10 chimeric)
     other1   150 reads, other2 150 reads  (one sample on two flow cells)
     othercram 200 reads
  reads/fastq.fastq.gz             250 reads; comments mix SAM tags and non-SAM tokens,
                                   space-separated (runid=... qs:f:... MM:Z:... ML:B:C,...)
  expected.json                    per input: read names and their MM/ML strings
"""
import argparse
import array
import gzip
import json
import os
import random

import pysam

COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(COMP)[::-1]


def mm_ml(seq, rng):
    """MM/ML for every third C of the read (original orientation, SAM spec MM skip counts)."""
    n_c = seq.count("C")
    mods = list(range(2, n_c, 3))
    skips, prev = [], -1
    for m in mods:
        skips.append(m - prev - 1)
        prev = m
    mm = "C+m?" + "".join(f",{k}" for k in skips) + ";"
    ml = [rng.randrange(0, 256) for _ in mods]
    return mm, ml


def make_reads(ref, n, rng, prefix, n_random=0, n_chimeric=0):
    names = list(ref)
    reads = []
    for i in range(n):
        if i < n_random:
            seq = "".join(rng.choice("ACGT") for _ in range(rng.randrange(1500, 4000)))
        else:
            parts = 2 if i < n_random + n_chimeric else 1
            seq = ""
            for _ in range(parts):
                c = rng.choice(names)
                L = rng.randrange(2000, 7000)
                s = rng.randrange(0, len(ref[c]) - L)
                frag = ref[c][s:s + L]
                frag = "".join(b if rng.random() > 0.01 else rng.choice("ACGT") for b in frag)
                seq += revcomp(frag) if rng.random() < 0.5 else frag
        qual = [rng.randrange(8, 40) for _ in seq]
        mm, ml = mm_ml(seq, rng)
        reads.append({"name": f"{prefix}_{i:04d}", "seq": seq, "qual": qual, "mm": mm, "ml": ml,
                      "qs": round(sum(qual) / len(qual), 3)})
    return reads


def write_ubam(path, reads):
    header = {"HD": {"VN": "1.6", "SO": "unknown"},
              "RG": [{"ID": "run1_dorado", "PL": "ONT", "DS": "basecall_model=synthetic", "SM": "lab_name_not_sample"}],
              "PG": [{"ID": "basecaller", "PN": "dorado", "VN": "0.8.3"}]}
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for r in reads:
            a = pysam.AlignedSegment(out.header)
            a.query_name, a.flag = r["name"], 4
            a.query_sequence, a.query_qualities = r["seq"], pysam.qualitystring_to_array("".join(chr(q + 33) for q in r["qual"]))
            a.set_tag("MM", r["mm"], "Z")
            a.set_tag("ML", array.array("B", r["ml"]))          # B:C (uint8) array
            a.set_tag("MN", len(r["seq"]), "i")
            a.set_tag("qs", float(r["qs"]), "f")
            a.set_tag("RG", "run1_dorado", "Z")
            out.write(a)


def write_fastq(path, reads, rng):
    with gzip.open(path, "wt") as fh:
        for r in reads:
            ml = "ML:B:C" + "".join(f",{v}" for v in r["ml"])
            toks = [f"runid=abc{rng.randrange(1000)}", f"qs:f:{r['qs']}", f"MM:Z:{r['mm']}", ml, "ch=12", "start_time=2026-01-01T00:00:00Z"]
            fh.write(f"@{r['name']} {' '.join(toks)}\n{r['seq']}\n+\n{''.join(chr(q + 33) for q in r['qual'])}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--seed", type=int, default=7)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    os.makedirs(f"{a.out}/ref", exist_ok=True)
    os.makedirs(f"{a.out}/reads", exist_ok=True)
    ref = {"chrA": "".join(rng.choice("ACGT") for _ in range(150_000)),
           "chrB": "".join(rng.choice("ACGT") for _ in range(100_000))}
    for fname, names in (("ref.fa", {"chrA": "chrA", "chrB": "chrB"}), ("other.fa", {"chrA": "A", "chrB": "B"})):
        with open(f"{a.out}/ref/{fname}", "w") as fh:
            for c, s in ref.items():
                fh.write(f">{names[c]}\n" + "\n".join(s[i:i + 80] for i in range(0, len(s), 80)) + "\n")
        pysam.faidx(f"{a.out}/ref/{fname}")
    expected = {}
    sets = {"modbam": make_reads(ref, 400, rng, "mod", n_random=12, n_chimeric=10),
            "other1": make_reads(ref, 150, rng, "fc1", n_chimeric=10),
            "other2": make_reads(ref, 150, rng, "fc2", n_chimeric=10),
            "othercram": make_reads(ref, 200, rng, "cr", n_random=5, n_chimeric=10),
            "fastq": make_reads(ref, 250, rng, "fq", n_random=4)}
    for name, reads in sets.items():
        if name == "fastq":
            write_fastq(f"{a.out}/reads/fastq.fastq.gz", reads, rng)
        else:
            write_ubam(f"{a.out}/reads/{name}.mod.bam", reads)
        expected[name] = {r["name"]: {"mm": r["mm"], "ml": r["ml"]} for r in reads}
    json.dump(expected, open(f"{a.out}/expected.json", "w"))
    print({k: len(v) for k, v in expected.items()})


if __name__ == "__main__":
    main()
