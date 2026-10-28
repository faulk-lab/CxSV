#!/usr/bin/env python3
"""Delete the .sra archive for samples whose CRAM is verified - and nothing else.

Implements Chris's day-one storage decision: alignments are kept as CRAM and the
raw .sra is dropped once its alignment is verified. A .sra is removed only when
ALL of these hold at the moment of deletion:

  1. datasets/logs/align/integrity/<ACC>.tsv exists and says read_count_match,
     with the .sra read count equal to the CRAM primary-alignment count;
  2. the final CRAM named in that record exists and passes samtools quickcheck;
  3. the CRAM's size on disk equals the size recorded at verification - this
     ties the record to this exact file, so a CRAM replaced or truncated since
     it was verified cannot clear the gate;
  4. no alignment is in progress for the accession (.cram.tmp absent).

Default is a dry run. Pass --delete to act. Deleted files go to the NAS
recycle bin (/mnt/sargasso1/#recycle) and free no space until it is emptied.
Every deletion is appended to logs/sra_deletions.<host>.tsv - per-host,
because several hosts appending to one file on SMB corrupts it.
"""
from __future__ import annotations
import argparse, os, re, socket, subprocess, sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path("/mnt/sargasso1/CxSV/datasets")
SRA = ROOT / "1kGP_ONT_1019_sra"
OUT = ROOT / "1kGP_ONT_1019_CHM13"
IDIR = ROOT / "logs/align/integrity"
LOG = ROOT / f"logs/sra_deletions.{socket.gethostname()}.tsv"
ENV = "source ~/miniforge3/etc/profile.d/conda.sh && conda activate sratools && "


def quickcheck(p: Path) -> bool:
    return subprocess.run(["bash", "-lc", ENV + f"samtools quickcheck {p}"],
                          capture_output=True).returncode == 0


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--delete", action="store_true", help="actually delete (default: dry run)")
    a = ap.parse_args()

    tally = {}
    freed = 0
    new = not LOG.exists()
    logf = open(LOG, "a") if a.delete else None
    if logf and new:
        logf.write("ts\taccession\tsample\tsra_bytes\tcram\tcram_bytes\tsra_reads\tcram_reads\n")

    for sra in sorted(SRA.glob("ERR*.sra")):
        acc = sra.stem
        rec = IDIR / f"{acc}.tsv"

        def keep(reason):
            tally[reason] = tally.get(reason, 0) + 1

        if not rec.exists():
            keep("kept: no integrity record yet"); continue
        f = rec.read_text().strip().split("\t")
        if len(f) != 6 or f[4] != "read_count_match" or f[2] != f[3] or not f[2].isdigit():
            keep("kept: record is not an exact read-count match"); continue
        sample, rec_size = f[0], int(f[5])
        cram = OUT / f"{sample}.{acc}.hs1.minimap2.cram"
        if not cram.exists():
            keep("kept: CRAM missing"); continue
        if Path(f"{cram}.tmp").exists():
            keep("kept: alignment in progress"); continue
        if cram.stat().st_size != rec_size:
            keep("kept: CRAM size differs from verified size"); continue
        if not quickcheck(cram):
            keep("kept: CRAM fails quickcheck"); continue

        size = sra.stat().st_size
        if a.delete:
            sra.unlink()
            logf.write(f"{datetime.now(timezone.utc).isoformat()}\t{acc}\t{sample}\t{size}"
                       f"\t{cram.name}\t{rec_size}\t{f[2]}\t{f[3]}\n")
            logf.flush()
        freed += size
        keep("deleted" if a.delete else "would delete")

    if logf:
        logf.close()
    for k, v in sorted(tally.items()):
        print(f"  {k}: {v}")
    print(f"{'freed (into recycle bin)' if a.delete else 'would free'}: {freed/1e12:.2f} TB")


if __name__ == "__main__":
    main()
