#!/usr/bin/env python3
"""genomes.tsv for segment_all.py / link_all.py (roadmap item 15, Phase 1).

  make_genomes.py --subset subset100.tsv --depth-root DEPTH --panel-samples DIR \\
      --extra KEY,STRATUM,SRC [--extra ...] [--only-done] --out genomes.tsv

Subset genomes (role panel): key <cohort>/<sample_id>, src DIR/<cohort>/<sample_id>.
Extras (role patient, e.g. the positive control and HG01868): key extra/<name>.
Sex: the subset sheet's sex, checked against the depth (chrX / autosome median
ratio, male ~0.5, female ~1.0); extras take the depth-inferred sex. Mismatches are
reported on stderr and in column sex_check.
--only-done: skip genomes whose depth pass has not finished (no depth.ok).
"""
import argparse
import csv
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cnlib  # noqa: E402


def depth_sex(path):
    fine = cnlib.read_bins(path, contigs=cnlib.AUTOSOMES + ["chrX"])
    co = cnlib.to_coarse(fine)
    g = cnlib.genome_median(co)
    x = co["chrX"][(co["chrX"] > 0) & np.isfinite(co["chrX"])]
    r = float(np.median(x)) / g
    return ("male" if r < 0.75 else "female"), round(r, 3)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subset", required=True)
    ap.add_argument("--depth-root", required=True)
    ap.add_argument("--panel-samples", required=True)
    ap.add_argument("--extra", action="append", default=[])
    ap.add_argument("--only-done", action="store_true")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rows = []
    for r in csv.DictReader(open(a.subset), delimiter="\t"):
        key = f"{r['cohort']}/{r['sample_id']}"
        rows.append(dict(key=key, cohort=r["cohort"], sample_id=r["sample_id"], individual_id=r["individual_id"],
                         stratum=r["stratum"], sex=r["sex"], flagged=r["flagged"], role="panel",
                         src=f"{a.panel_samples}/{key}"))
    for e in a.extra:
        key, stratum, src = e.split(",")
        name = key.split("/")[-1]
        rows.append(dict(key=key, cohort="extra", sample_id=name, individual_id=name, stratum=stratum,
                         sex="", flagged="False", role="patient", src=src))
    out = []
    for r in rows:
        d = f"{a.depth_root}/{r['key']}"
        if not os.path.exists(f"{d}/depth.ok"):
            if a.only_done:
                continue
            sys.exit(f"depth pass not finished: {d}")
        s, ratio = depth_sex(f"{d}/depth.regions.bed.gz")
        r["x_ratio"] = ratio
        if not r["sex"]:
            r["sex"], r["sex_check"] = s, "inferred"
        else:
            r["sex_check"] = "ok" if s == r["sex"] else f"MISMATCH(depth={s})"
            if s != r["sex"]:
                print(f"sex mismatch {r['key']}: sheet {r['sex']}, depth {s} (chrX ratio {ratio})", file=sys.stderr)
        out.append(r)
    cols = ["key", "cohort", "sample_id", "individual_id", "stratum", "sex", "flagged", "role", "src", "x_ratio", "sex_check"]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out)
    print(f"{len(out)} genomes -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
