#!/usr/bin/env python3
"""Snapshot of a build tree, to prove a downstream-only rebuild left the old build
unchanged (docs/PANEL_BUILD.md, "Downstream-only rebuild").

    tree_snapshot.py ROOT --out SNAPSHOT.tsv [--sha256 SUBDIR ...]
    tree_snapshot.py --compare BEFORE.tsv AFTER.tsv

Every file and directory under ROOT: relative path, type, size, mtime_ns, inode
(lstat, links not followed); for files under each --sha256 SUBDIR (relative to
ROOT, e.g. results/panel/release) also the sha256. --compare prints every
difference and exits 1 if there is one.
"""
import argparse
import hashlib
import os
import sys

COLS = ("path", "type", "size", "mtime_ns", "inode", "sha256")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def snapshot(root, sha_dirs=()):
    root = os.path.abspath(root)
    sha_dirs = [os.path.normpath(d) for d in sha_dirs]
    rows = []
    for d, dirs, files in os.walk(root):
        dirs.sort()
        for name in sorted(dirs) + sorted(files):
            p = os.path.join(d, name)
            rel = os.path.relpath(p, root)
            st = os.lstat(p)
            kind = "dir" if os.path.isdir(p) and not os.path.islink(p) else ("link" if os.path.islink(p) else "file")
            h = ""
            if kind == "file" and any(rel == s or rel.startswith(s + os.sep) for s in sha_dirs):
                h = sha256(p)
            rows.append((rel, kind, str(st.st_size if kind != "dir" else 0), str(st.st_mtime_ns), str(st.st_ino), h))
    return rows


def read(path):
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        assert tuple(head) == COLS, head
        return {l.split("\t", 1)[0]: l.rstrip("\n") for l in fh}


def compare(a, b):
    A, B = read(a), read(b)
    diffs = [f"only before: {k}" for k in sorted(set(A) - set(B))]
    diffs += [f"only after: {k}" for k in sorted(set(B) - set(A))]
    diffs += [f"changed: {A[k]}  ->  {B[k]}" for k in sorted(set(A) & set(B)) if A[k] != B[k]]
    return diffs


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("root", nargs="?")
    ap.add_argument("--out")
    ap.add_argument("--sha256", action="append", default=[])
    ap.add_argument("--compare", nargs=2, metavar=("BEFORE", "AFTER"))
    a = ap.parse_args(argv)
    if a.compare:
        d = compare(*a.compare)
        print("\n".join(d) if d else "identical")
        print(f"{len(d)} difference(s)", file=sys.stderr)
        return 1 if d else 0
    if not (a.root and a.out):
        ap.error("ROOT and --out are required (or --compare)")
    rows = snapshot(a.root, a.sha256)
    tmp = a.out + ".tmp"
    with open(tmp, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        fh.writelines("\t".join(r) + "\n" for r in rows)
    os.replace(tmp, a.out)
    print(f"{len(rows)} entries ({sum(1 for r in rows if r[5])} hashed) -> {a.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
