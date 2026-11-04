#!/usr/bin/env python3
"""Fail a sample whose SVIM run stopped reading its alignment early.

    check_svim_output.py --svim-vcf V --svim-log L --reference-vcf R --out OUT
                         [--min-reference-calls 20]

SVIM 2.0.0 reads CRAMs through pysam without --reference: htslib finds each
contig's sequence by its @SQ M5 tag in REF_PATH (the MD5 cache). When a cache
file was missing (chr7 of hs1, 2026-10-31, panel-v1.0 build), SVIM logged a
CRAM decode / "truncated file" error and ~269 runs either failed or ended with
no calls on the contigs after the failure point, still exiting 0. The caller
rules and outputs are fingerprinted and stay as they are; this check runs after
them and fails the sample instead (exit 1, no OUT), when

  1. log: the SVIM log has a decode / truncation / reference error, or lacks
     SVIM's closing "Done." line;
  2. contigs: a contig on which the reference caller (Sniffles2) has
     >= --min-reference-calls records has none in the SVIM VCF. Contigs are
     taken in SVIM's own order (its ##contig header lines = the alignment's
     @SQ order); the report names the last contig with SVIM calls, i.e. the
     point after which the output was cut. Contigs the reference caller hardly
     calls on (chrY of a female genome, chrM, contigs outside a sliced test
     alignment) are not judged.

OUT (TSV): one row per contig (svim_calls, reference_calls, status) after
`# key<TAB>value` summary lines. Not part of the calling fingerprint.
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import Counter

BAD_LOG = re.compile(r"truncated file|failure to decode|decode failure|failed to decode|"
                     r"failed to populate reference|failed to open reference|unable to fetch reference|"
                     r"EOF marker is absent|\[E::|Traceback \(most recent call last\)", re.I)
DONE = re.compile(r"\]\s+Done\.\s*$")


def contig_counts(path):
    """(contig order from ##contig lines, Counter of records per contig)."""
    order, n = [], Counter()
    with open(path, "rb") as fh:
        for line in fh:
            if line[:1] == b"#":
                if line.startswith(b"##contig=<"):
                    m = re.search(rb"ID=([^,>]+)", line)
                    if m:
                        order.append(m.group(1).decode())
                continue
            tab = line.find(b"\t")
            if tab > 0:
                n[line[:tab].decode()] += 1
    return order, n


def check_log(path):
    """Problems found in the SVIM log (empty = clean)."""
    try:
        text = open(path, errors="replace").read()
    except OSError as e:
        return [f"SVIM log unreadable: {e}"]
    probs = []
    bad = [ln.strip() for ln in text.splitlines() if BAD_LOG.search(ln)]
    if bad:
        probs.append(f"SVIM log has {len(bad)} decode/truncation/reference error line(s), first: {bad[0][:200]}")
    tail = [ln for ln in text.splitlines() if ln.strip()][-5:]
    if not any(DONE.search(ln) for ln in tail):
        probs.append("SVIM log does not end with SVIM's 'Done.' line")
    return probs


def check(svim_vcf, svim_log, ref_vcf, min_ref_calls=20):
    order, sv = contig_counts(svim_vcf)
    ref_order, ref = contig_counts(ref_vcf)
    for c in ref_order + sorted(set(sv) | set(ref)):
        if c not in order:
            order.append(c)
    rows, missing = [], []
    for c in order:
        judged = ref[c] >= min_ref_calls
        status = "ok" if sv[c] > 0 else ("MISSING" if judged else "not_judged")
        if status == "MISSING":
            missing.append(c)
        rows.append((c, sv[c], ref[c], status))
    with_calls = [c for c in order if sv[c] > 0]
    last = with_calls[-1] if with_calls else None
    probs = check_log(svim_log)
    if missing:
        after = [c for c in missing if last is None or order.index(c) > order.index(last)]
        where = (f"; output ends at {last}: no calls on {len(after)} later contig(s) ({','.join(after[:6])}"
                 f"{',...' if len(after) > 6 else ''})" if after else "")
        probs.append(f"SVIM has no calls on {len(missing)} contig(s) where the reference caller has "
                     f">= {min_ref_calls}: {','.join(missing[:10])}{',...' if len(missing) > 10 else ''}{where}")
    summary = {"svim_vcf": svim_vcf, "svim_log": svim_log, "reference_vcf": ref_vcf,
               "min_reference_calls": min_ref_calls, "svim_calls": sum(sv.values()),
               "reference_calls": sum(ref.values()), "last_contig_with_svim_calls": last or "",
               "contigs_missing": ",".join(missing), "status": "FAIL" if probs else "PASS"}
    return probs, summary, rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--svim-vcf", required=True)
    ap.add_argument("--svim-log", required=True)
    ap.add_argument("--reference-vcf", required=True, help="Sniffles2 VCF of the same sample")
    ap.add_argument("--min-reference-calls", type=int, default=20)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    probs, summary, rows = check(a.svim_vcf, a.svim_log, a.reference_vcf, a.min_reference_calls)
    for k, v in summary.items():
        print(f"{k}\t{v}", file=sys.stderr)
    if probs:
        for p in probs:
            print(f"FAIL: {p}", file=sys.stderr)
        print("The SVIM output of this sample is incomplete. For a CRAM, check that REF_PATH / the MD5 "
              "reference cache holds every contig (build_ref_cache.py), then rerun run_svim.", file=sys.stderr)
        return 1
    with open(a.out, "w") as fh:
        for k, v in summary.items():
            fh.write(f"# {k}\t{v}\n")
        fh.write("contig\tsvim_calls\treference_calls\tstatus\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
