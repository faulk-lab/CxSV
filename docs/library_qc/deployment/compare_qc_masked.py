#!/usr/bin/env python3
"""Compare two patient result trees file by file, masking only the library-QC additions.

Usage: compare_qc_masked.py OLD_DIR NEW_DIR OUT_TSV [--old-root S --new-root S]...

A file is
  identical          same bytes (after replacing the run's absolute root paths),
  identical_masked   equal once the library-QC parts are removed:
                       *.tsv   columns Library_QC_Gate and the QC fields,
                       *.md    the "library QC" bullet and its sub-lines, and the
                               "Library QC: loci demoted by conservative mode" section,
                       *.json  the "library_qc" key,
  qc_new             present only in NEW and part of the library QC (path has library_qc or /qc/),
  DIFFERENT / ONLY_OLD / ONLY_NEW_UNEXPECTED otherwise (these fail the comparison).
Skipped: logs/, benchmarks/, .snakemake/, .tmp/, tmp/ (run bookkeeping, timings).
Exit 0 when nothing fails."""
import csv, gzip, hashlib, io, json, os, re, sys

SKIP_DIRS = {"logs", "benchmarks", ".snakemake", ".tmp", "tmp", "log"}
QC_COLS = {"Library_QC_Gate", "qc_group", "qc_flag", "qc_reason", "frac_sa", "frac_sa_interchrom",
           "svim_only_bnd_per_gb", "svim_only_any_per_gb", "z_svim_only_bnd_per_gb",
           "z_svim_only_any_per_gb", "coverage_est", "reads_sampled"}


def walk(root):
    out = {}
    for d, dirs, files in os.walk(root):
        dirs[:] = [x for x in dirs if x not in SKIP_DIRS]
        for f in files:
            p = os.path.join(d, f)
            out[os.path.relpath(p, root)] = p
    return out


def read(p, roots):
    b = open(p, "rb").read()
    if p.endswith(".gz"):
        try:
            b = gzip.decompress(b)
        except OSError:
            pass
    for r in roots:
        b = b.replace(r.encode(), b"<ROOT>")
    return b


def mask_tsv(b):
    rows = list(csv.reader(io.StringIO(b.decode()), delimiter="\t"))
    if not rows:
        return rows
    keep = [i for i, c in enumerate(rows[0]) if c not in QC_COLS]
    return [[r[i] for i in keep if i < len(r)] for r in rows]


def mask_md(b):
    out, skip_section, in_bullet = [], False, False
    for line in b.decode().splitlines():
        if line.startswith("## "):
            skip_section = line.startswith("## Library QC")
        if skip_section:
            continue
        if line.startswith("- library QC"):
            in_bullet = True
            continue
        if in_bullet and line.startswith("  "):
            continue
        in_bullet = False
        if line.startswith("- library QC gate:"):
            continue
        # PR #28 appended the conservative-mode clause to the QC bullet; nothing else to strip
        out.append(line)
    return "\n".join(out).strip()


def mask_json(b):
    try:
        d = json.loads(b)
    except ValueError:
        return None
    if isinstance(d, dict):
        d.pop("library_qc", None)
    return json.dumps(d, sort_keys=True)


def main():
    old, new, out_tsv = sys.argv[1:4]
    roots_old, roots_new = [os.path.abspath(old)], [os.path.abspath(new)]
    a = sys.argv[4:]
    while a:
        k, v = a[0], a[1]; a = a[2:]
        (roots_old if k == "--old-root" else roots_new).append(v)
    O, N = walk(old), walk(new)
    rows, fails = [], 0
    for rel in sorted(set(O) | set(N)):
        if rel not in N:
            st = "ONLY_OLD"
        elif rel not in O:
            st = "qc_new" if ("library_qc" in rel or "/qc/" in "/" + rel) else "ONLY_NEW_UNEXPECTED"
        else:
            bo, bn = read(O[rel], roots_old), read(N[rel], roots_new)
            if bo == bn:
                st = "identical"
            else:
                base = rel[:-3] if rel.endswith(".gz") else rel
                if base.endswith(".tsv"):
                    st = "identical_masked" if mask_tsv(bo) == mask_tsv(bn) else "DIFFERENT"
                elif base.endswith(".md"):
                    st = "identical_masked" if mask_md(bo) == mask_md(bn) else "DIFFERENT"
                elif base.endswith(".json"):
                    mo, mn = mask_json(bo), mask_json(bn)
                    st = "identical_masked" if mo is not None and mo == mn else "DIFFERENT"
                else:
                    st = "DIFFERENT"
        fails += st.isupper()
        rows.append((rel, st))
    with open(out_tsv, "w") as fh:
        fh.write("file\tstatus\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    from collections import Counter
    c = Counter(s for _, s in rows)
    print(f"{old} vs {new}: " + ", ".join(f"{k} {v}" for k, v in sorted(c.items())))
    for rel, st in rows:
        if st.isupper() or st in ("identical_masked", "qc_new"):
            print(f"  {st}\t{rel}")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
