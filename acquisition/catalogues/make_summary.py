#!/usr/bin/env python3
"""Collect every <catalogue>.<version>.manifest.tsv under --root into one
summary table, <root>/external_catalogues_summary.tsv (columns below).

    python3 make_summary.py --root /mnt/sargasso1/CxSV/datasets/external_catalogues
"""
import argparse
import glob
import os

COLS = ["catalogue", "version", "native_build", "output_build", "licence",
        "n_DEL", "n_DUP", "n_INV", "n_INS", "n_BND", "n_total",
        "liftover_loss_pct", "junction_sites_path", "output_sha256", "code_version"]


def read_manifest(path):
    with open(path) as fh:
        next(fh)
        return dict(line.rstrip("\n").split("\t", 1) for line in fh if "\t" in line)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", required=True)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    out = a.out or os.path.join(a.root, "external_catalogues_summary.tsv")
    rows = []
    for mp in sorted(glob.glob(os.path.join(a.root, "*", "*", "*.manifest.tsv"))):
        m = read_manifest(mp)
        d = os.path.dirname(mp)
        rows.append([m["source"], m["version"], m["native_build"], m["output_build"],
                     m["licence"]] + [m.get(f"n_sites_{t}", "0") for t in
                                      ("DEL", "DUP", "INV", "INS", "BND")] +
                    [m["n_sites_total"], m.get("liftover_loss_pct", "NA (native hs1)"),
                     os.path.join(d, m["output_file"]), m["output_sha256"], m["code_version"]])
    with open(out + ".tmp", "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    os.replace(out + ".tmp", out)
    print(out)


if __name__ == "__main__":
    main()
