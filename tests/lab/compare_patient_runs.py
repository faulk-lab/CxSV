#!/usr/bin/env python3
"""Compare two patient-workflow result trees record by record (bundle validation, roadmap item 2).

    compare_patient_runs.py <results_root_A> <results_root_B> [--out table.tsv]

Every file below either root is compared, except logs/ and .snakemake/:
  VCF   records only ('##' meta lines carry dates and command lines)
  PDF   after dropping /CreationDate, /ModDate and /ID lines (render timestamps)
  other byte for byte
Two run records are expected to differ and are listed as such: panel_release_check.txt
(release path, runtime image) and panel_compatibility.txt (release path); their
status lines must still agree. Exit 1 when any other file differs or exists on one side only.
"""
import argparse
import os
import re
import sys

EXPECTED = {"panel_release_check.txt": "status:", "panel_compatibility.txt": "status:"}
PDF_VOLATILE = re.compile(rb"^(/CreationDate|/ModDate|/ID|<< /CreationDate).*$", re.M)


def files(root):
    out = set()
    for d, dirs, fs in os.walk(root):
        dirs[:] = [x for x in dirs if x not in ("logs", ".snakemake")]
        for f in fs:
            out.add(os.path.relpath(os.path.join(d, f), root))
    return out


def content(path):
    data = open(path, "rb").read()
    if path.endswith(".vcf"):
        return b"\n".join(l for l in data.split(b"\n") if not l.startswith(b"##")), "records"
    if path.endswith(".pdf"):
        return PDF_VOLATILE.sub(b"", re.sub(rb"/CreationDate \(D:[^)]*\)|/ModDate \(D:[^)]*\)", b"", data)), \
            "pdf minus timestamps"
    return data, "bytes"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--out")
    a = ap.parse_args()
    fa, fb = files(a.a), files(a.b)
    rows, bad = [], 0
    for rel in sorted(fa | fb):
        if rel not in fa or rel not in fb:
            rows.append((rel, "only in " + ("A" if rel in fa else "B"), "-"))
            bad += 1
            continue
        ca, how = content(os.path.join(a.a, rel))
        cb, _ = content(os.path.join(a.b, rel))
        base = os.path.basename(rel)
        if ca == cb:
            rows.append((rel, "identical", how))
        elif base in EXPECTED:
            sa = [l for l in ca.decode().splitlines() if l.startswith(EXPECTED[base])]
            sb = [l for l in cb.decode().splitlines() if l.startswith(EXPECTED[base])]
            ok = sa == sb and sa
            rows.append((rel, "differs as expected (paths/runtime); same status" if ok else "STATUS DIFFERS", how))
            bad += 0 if ok else 1
        else:
            rows.append((rel, "DIFFERS", how))
            bad += 1
    text = "file\tresult\tcompared\n" + "".join("\t".join(r) + "\n" for r in rows)
    if a.out:
        open(a.out, "w").write(text)
    n_id = sum(r[1] == "identical" for r in rows)
    print(text, end="")
    print(f"# {len(rows)} files: {n_id} identical, {len(rows) - n_id - bad} expected run-record differences, "
          f"{bad} unexpected differences")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
