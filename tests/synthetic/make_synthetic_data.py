#!/usr/bin/env python3
"""Generate the synthetic end-to-end test dataset (fully artificial; no human data).

A seeded random 1.2 Mb contig ("chr20") is the reference. Four diploid-looking
(homozygous) genomes are derived from it by applying structural events, then
long reads are simulated from each with ONT-like errors and aligned with
minimap2 (map-ont). The expected outcome of the panel + patient workflows is
written to expected.json and checked by check_synthetic_results.py.

Events (reference coordinates, 0-based half-open):
  COMMON  (panel1, panel2, panel3, patient)  CxSV at 200-226 kb:
          DEL 200,000-205,000 + INV 210,000-218,000 + tandem DUP 222,000-226,000
  RARE    (panel2, patient)                   CxSV at 500-526 kb, same shape
  PRIVATE (patient only)                      nested CxSV at 800-829 kb:
          DEL 806,000-809,000 inside INV 800,000-820,000, + DUP 826,000-829,000
  DECOY   (panel1, panel3)                    a DIFFERENT CxSV at the PRIVATE locus (hotspot):
          DEL 802,000-804,000 + INV 812,000-818,000 + tandem DUP 821,000-824,000
  SIMPLE  (everyone)                          DEL 1,000,000-1,003,000 (not a CxSV)

Expected (breakpoint-level matching, roadmap item 3): the patient's CxSV loci
overlapping COMMON, RARE and PRIVATE have panel carriers 3, 1 and 0, so COMMON
and RARE are public and PRIVATE is private, although DECOY gives PRIVATE 2
region-overlap carriers (it would have been called public by region overlap).
The SIMPLE deletion's junction has 3 panel carriers.

Usage: make_synthetic_data.py <outdir> [--coverage 25] [--seed 7] [--threads 4]
Needs minimap2 and samtools on PATH. Writes:
  ref/chr20.fa(.fai)  resources/{genes,exons,repeats,segdup}.bed
  bams/<sample>.bam(.bai)  (panel3 additionally as CRAM)  expected.json
"""
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

CONTIG, LENGTH = "chr20", 1_200_000
COMP = bytes.maketrans(b"ACGT", b"TGCA")

EVENTS = {
    "COMMON":  [("DEL", 200_000, 205_000), ("INV", 210_000, 218_000), ("DUP", 222_000, 226_000)],
    "RARE":    [("DEL", 500_000, 505_000), ("INV", 510_000, 518_000), ("DUP", 522_000, 526_000)],
    # INVDEL: the 800-820 kb segment is inverted after deleting 806-809 kb inside it
    "PRIVATE": [("INVDEL", 800_000, 820_000, 806_000, 809_000), ("DUP", 826_000, 829_000)],
    "DECOY":   [("DEL", 802_000, 804_000), ("INV", 812_000, 818_000), ("DUP", 821_000, 824_000)],
    "SIMPLE":  [("DEL", 1_000_000, 1_003_000)],
}
SAMPLES = {
    "panel1":  ["COMMON", "DECOY", "SIMPLE"],
    "panel2":  ["COMMON", "RARE", "SIMPLE"],
    "panel3":  ["COMMON", "DECOY", "SIMPLE"],
    "patient": ["COMMON", "RARE", "PRIVATE", "SIMPLE"],
}
SPAN = {name: (min(e[1] for e in ev), max(e[2] for e in ev)) for name, ev in EVENTS.items()}


def revcomp(s: bytes) -> bytes:
    return s.translate(COMP)[::-1]


def apply_events(ref: bytes, groups) -> bytes:
    """Apply event groups, and the events inside each group, right to left,
    so every event's reference coordinates are still valid when it is applied."""
    seq = bytearray(ref)
    for name in sorted(groups, key=lambda g: -SPAN[g][0]):
        for ev in sorted(EVENTS[name], key=lambda ev: -ev[1]):
            kind, s, e = ev[:3]
            if kind == "DEL":
                del seq[s:e]
            elif kind == "INV":
                seq[s:e] = revcomp(bytes(seq[s:e]))
            elif kind == "DUP":                      # tandem: copy inserted right after
                seq[e:e] = seq[s:e]
            elif kind == "INVDEL":
                ds, de = ev[3], ev[4]
                seq[s:e] = revcomp(bytes(seq[s:ds]) + bytes(seq[de:e]))
            else:
                raise ValueError(kind)
    return bytes(seq)


def add_errors(r: np.ndarray, rng, bases) -> np.ndarray:
    """~2% substitutions, 1.5% deletions, 1.5% insertions, independently per base."""
    u = rng.random(len(r))
    r = r.copy()
    sub = u < 0.02
    r[sub] = bases[rng.integers(0, 4, int(sub.sum()))]
    keep = ~((u >= 0.02) & (u < 0.035))
    ins_after = ((u >= 0.035) & (u < 0.05))[keep]
    kept = r[keep]
    n_ins = int(ins_after.sum())
    out = np.empty(len(kept) + n_ins, dtype=np.uint8)
    ins_pos = np.flatnonzero(ins_after) + np.arange(n_ins) + 1   # slot after each flagged base
    is_ins = np.zeros(len(out), bool)
    is_ins[ins_pos] = True
    out[~is_ins] = kept
    out[is_ins] = bases[rng.integers(0, 4, n_ins)]
    return out


def simulate_reads(hap: bytes, coverage: float, rng, out_fq, name):
    """ONT-like reads: lognormal lengths (median ~8 kb, 2-40 kb), random strand."""
    bases = np.frombuffer(b"ACGT", dtype=np.uint8)
    n_bases, i = 0, 0
    target = coverage * len(hap)
    with open(out_fq, "w") as fq:
        while n_bases < target:
            L = int(np.clip(rng.lognormal(np.log(8000), 0.5), 2000, 40000))
            start = int(rng.integers(0, max(1, len(hap) - L)))
            frag = hap[start:start + L]
            if rng.random() < 0.5:
                frag = revcomp(frag)
            seq = add_errors(np.frombuffer(frag, dtype=np.uint8), rng, bases).tobytes().decode()
            fq.write(f"@{name}_{i}\n{seq}\n+\n{'5' * len(seq)}\n")
            n_bases += len(seq)
            i += 1
    return i


def run(cmd, **kw):
    print("+", cmd if isinstance(cmd, str) else " ".join(cmd), file=sys.stderr)
    subprocess.run(cmd, check=True, shell=isinstance(cmd, str), **kw)


def write_resources(res: Path):
    res.mkdir(parents=True, exist_ok=True)
    (res / "genes.bed").write_text(
        f"{CONTIG}\t150000\t260000\tSYNA\n{CONTIG}\t790000\t840000\tSYNP\n{CONTIG}\t990000\t1010000\tSYNS\n")
    (res / "exons.bed").write_text(
        f"{CONTIG}\t211000\t212000\tSYNA\n{CONTIG}\t805000\t807500\tSYNP\n{CONTIG}\t812000\t813000\tSYNP\n")
    (res / "repeats.bed").write_text(
        f"{CONTIG}\t300000\t301000\tLINE\n{CONTIG}\t600000\t600300\tSINE\n")
    (res / "segdup.bed").write_text(f"{CONTIG}\t650000\t660000\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("outdir", type=Path)
    ap.add_argument("--coverage", type=float, default=25)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()
    out = a.outdir
    for d in ("ref", "bams", "reads", "resources"):
        (out / d).mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(a.seed)
    gc = 0.41
    p = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    ref = bytes(np.frombuffer(b"ACGT", dtype=np.uint8)[rng.choice(4, LENGTH, p=p)])
    fa = out / "ref" / "chr20.fa"
    with open(fa, "w") as fh:
        fh.write(f">{CONTIG}\n")
        for k in range(0, LENGTH, 60):
            fh.write(ref[k:k + 60].decode() + "\n")
    run(["samtools", "faidx", str(fa)])
    write_resources(out / "resources")

    for sample, groups in SAMPLES.items():
        hap = apply_events(ref, groups)
        fq = out / "reads" / f"{sample}.fq"
        n = simulate_reads(hap, a.coverage, rng, fq, sample)
        bam = out / "bams" / f"{sample}.bam"
        run(f"minimap2 -t {a.threads} -ax map-ont -R '@RG\\tID:{sample}\\tSM:{sample}' {fa} {fq} 2>/dev/null"
            f" | samtools sort -@ 2 -o {bam} - && samtools index {bam}")
        print(f"{sample}: {len(hap):,} bp haplotype, {n:,} reads", file=sys.stderr)
        fq.unlink()
    run(f"samtools view -C -T {fa} -o {out}/bams/panel3.cram {out}/bams/panel3.bam "
        f"&& samtools index {out}/bams/panel3.cram && rm {out}/bams/panel3.bam {out}/bams/panel3.bam.bai")

    expected = {
        "contig": CONTIG,
        "events": {k: {"start": v[0], "end": v[1]} for k, v in SPAN.items()},
        "patient_expectations": {
            "COMMON":  {"cxsv": True, "panel_carriers": 3, "status": "public"},
            "RARE":    {"cxsv": True, "panel_carriers": 1, "status": "public"},
            "PRIVATE": {"cxsv": True, "panel_carriers": 0, "status": "private",
                        "region_carriers": 2, "hotspot": "yes"},
        },
        "junction_expectations": {
            "SIMPLE": {"svtype": "DEL", "pos1": 1_000_000, "pos2": 1_003_000, "panel_carriers": 3},
        },
        "samples": SAMPLES,
        "seed": a.seed, "coverage": a.coverage,
    }
    (out / "expected.json").write_text(json.dumps(expected, indent=2) + "\n")
    print(f"synthetic dataset written to {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
