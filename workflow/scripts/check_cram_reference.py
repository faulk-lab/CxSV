#!/usr/bin/env python3
"""Fail fast if a CRAM was not made against the pipeline's reference.

Compares every ``@SQ`` line's M5 tag in the CRAM header with the MD5 manifest
written by build_ref_cache.py. Reading a CRAM header needs no reference, so this
runs in seconds and turns an obscure mid-run decode failure (or, worse, a CRAM
silently decoded against a different assembly) into a clear error up front.

Usage: check_cram_reference.py <md5_manifest.tsv> <sample.cram>
Prints a one-line OK summary on success; exits 1 with the mismatches otherwise.
"""
import subprocess
import sys


def main(manifest, cram):
    ref = {}
    with open(manifest) as fh:
        next(fh)
        for line in fh:
            sn, ln, m5 = line.rstrip("\n").split("\t")
            ref[m5] = (sn, int(ln))
    hdr = subprocess.run(["samtools", "view", "-H", cram], check=True,
                         capture_output=True, text=True).stdout
    sq = [dict(f.split(":", 1) for f in l.split("\t")[1:] if ":" in f)
          for l in hdr.splitlines() if l.startswith("@SQ")]
    if not sq:
        sys.exit(f"{cram}: no @SQ lines in header")
    problems = []
    for t in sq:
        m5 = t.get("M5", "").lower()
        if not m5:
            problems.append(f"{t.get('SN')}: no M5 tag (cannot verify; re-encode with samtools view -T <ref>)")
        elif m5 not in ref:
            problems.append(f"{t.get('SN')}: M5 {m5} is not in the pipeline reference")
        elif ref[m5][1] != int(t.get("LN", -1)):
            problems.append(f"{t.get('SN')}: length {t.get('LN')} != reference {ref[m5][1]}")
    if problems:
        sys.exit(f"{cram}: {len(problems)} @SQ line(s) do not match the reference "
                 f"in {manifest}:\n  " + "\n  ".join(problems[:10]))
    print(f"OK {cram}: all {len(sq)} @SQ M5 tags match {manifest}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
