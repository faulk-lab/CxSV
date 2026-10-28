#!/usr/bin/env python3
"""Rebuild per-sample integrity records from ground truth.

The shared align_integrity.tsv was corrupted by concurrent cross-host appends
(NUL-filled rows, rows lost). This rebuilds datasets/logs/align/integrity/,
one file per accession, without trusting that file:

1. SALVAGE a legacy row only if every check passes: exactly 6 fields, no NUL
   bytes, the sample.accession names a CRAM that exists, both read counts are
   integers and equal, the verdict is read_count_match, and - the decisive
   check - the recorded size equals the CRAM's current size on disk. The size
   ties the row to that exact file, so a row spliced from two writes or left
   over from an earlier version of the file cannot pass.
2. RE-VERIFY every other CRAM from scratch: samtools quickcheck, rebuild a
   missing index, then primary-alignment count vs the read count sra-stat
   reports for the .sra.

CRAMs written within --min-age-h are skipped: their own align script may still
be indexing them, and two writers on one .crai is its own corruption risk.
Re-run until it reports zero skipped.

Records are written to a temp name and renamed, so this is safe to run while
alignments are completing.
"""
from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path("/mnt/sargasso1/CxSV/datasets")
OUT = ROOT / "1kGP_ONT_1019_CHM13"
SRA = ROOT / "1kGP_ONT_1019_sra"
LEGACY = ROOT / "logs/align/align_integrity.tsv"
IDIR = ROOT / "logs/align/integrity"
REF = "/home/minknow/refs/hs1/hs1.fa"
ENV = "source ~/miniforge3/etc/profile.d/conda.sh && conda activate sratools && "
NAME = re.compile(r"^([A-Z]{2}\d{5})\.(ERR\d+)\.hs1\.minimap2\.cram$")


def write_record(sample, acc, sra_n, cram_n, verdict, size) -> None:
    tmp = IDIR / f".{acc}.{os.getpid()}.tmp"
    tmp.write_text(f"{sample}\t{acc}\t{sra_n}\t{cram_n}\t{verdict}\t{size}\n")
    os.replace(tmp, IDIR / f"{acc}.tsv")


def sh(cmd: str, timeout: int = 5400) -> subprocess.CompletedProcess:
    return subprocess.run(["bash", "-lc", ENV + cmd], capture_output=True,
                          text=True, timeout=timeout)


def verify(sample: str, acc: str, cram: Path) -> tuple[str, str]:
    if sh(f"samtools quickcheck {cram}").returncode != 0:
        return acc, "QUICKCHECK_FAIL"
    crai = Path(f"{cram}.crai")
    if not crai.exists() or crai.stat().st_size == 0:
        if sh(f"samtools index -@ 4 {cram}").returncode != 0:
            return acc, "INDEX_FAIL"
    p1 = sh(f"sra-stat --quick {SRA}/{acc}.sra | awk -F'|' '{{print $3}}' | awk -F: '{{print $1}}'")
    p2 = sh(f"samtools view -c -F 0x900 --reference {REF} {cram}")
    a, b = p1.stdout.strip(), p2.stdout.strip()
    if not (a.isdigit() and b.isdigit()):
        return acc, "UNVERIFIED"
    verdict = "read_count_match" if a == b else "READ_COUNT_MISMATCH"
    write_record(sample, acc, a, b, verdict, cram.stat().st_size)
    return acc, verdict


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--min-age-h", type=float, default=3.0)
    a = ap.parse_args()
    IDIR.mkdir(parents=True, exist_ok=True)

    crams = {}
    for f in OUT.glob("*.cram"):
        m = NAME.match(f.name)
        if m:
            crams[m.group(2)] = (m.group(1), f)

    have = {p.stem for p in IDIR.glob("ERR*.tsv")}

    # --- 1. strict salvage of the legacy shared file ---
    salvaged = rejected = 0
    if LEGACY.exists():
        for raw in LEGACY.read_bytes().split(b"\n"):
            if not raw:
                continue
            if b"\x00" in raw:
                rejected += 1
                continue
            f = raw.decode(errors="replace").split("\t")
            ok = (len(f) == 6 and f[1] in crams and crams[f[1]][0] == f[0]
                  and f[2].isdigit() and f[2] == f[3] and f[4] == "read_count_match"
                  and f[5].isdigit() and int(f[5]) == crams[f[1]][1].stat().st_size)
            if not ok:
                rejected += 1
                continue
            if f[1] not in have:
                write_record(*f[:5], f[5])
                have.add(f[1])
                salvaged += 1

    # --- 2. re-verify everything still lacking a record ---
    now = time.time()
    todo, young = [], []
    for acc, (sample, f) in sorted(crams.items()):
        if acc in have:
            continue
        if now - f.stat().st_mtime < a.min_age_h * 3600:
            young.append(acc)
        else:
            todo.append((sample, acc, f))

    print(f"CRAMs={len(crams)} records_before={len(have) - salvaged} "
          f"salvaged={salvaged} legacy_rows_rejected={rejected} "
          f"to_verify={len(todo)} skipped_too_young={len(young)}", flush=True)

    results = {}
    with ThreadPoolExecutor(a.jobs) as ex:
        for acc, v in ex.map(lambda t: verify(*t), todo):
            results[acc] = v
            print(f"  {acc}: {v}", flush=True)

    recs = list(IDIR.glob("ERR*.tsv"))
    match = sum(1 for r in recs if "\tread_count_match\t" in r.read_text())
    print(f"--- records={len(recs)} read_count_match={match} "
          f"not_match={len(recs) - match} CRAMs_without_record="
          f"{len(set(crams) - {r.stem for r in recs})} skipped_young={young}")
    bad = {k: v for k, v in results.items() if v != "read_count_match"}
    if bad:
        print("NEEDS ATTENTION:", bad)
        sys.exit(1)


if __name__ == "__main__":
    main()
