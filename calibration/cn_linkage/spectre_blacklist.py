#!/usr/bin/env python3
"""Blacklist BED for Spectre from a Phase 1 PoN (roadmap item 15 Phase 2, step 4).

    spectre_blacklist.py PON.npz OUT.bed [--bin 10000]

The Phase 1 PoN (segment_all.py; per contig a bins x genomes matrix of normalised 10 kb
depth) is masked with the copy-number track's rules (cxsv_cn.DEFAULTS): PoN median
< 0.25, relative MAD > 0.30 or fewer than 4 genomes with data. Adjacent masked bins
are merged. hs1 has no N runs, so without this Spectre sees centromeric and satellite
arrays at depth ~0 as deletions.
"""
import argparse

import numpy as np


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("pon")
    ap.add_argument("out")
    ap.add_argument("--bin", type=int, default=10000)
    ap.add_argument("--min-median", type=float, default=0.25)
    ap.add_argument("--max-relmad", type=float, default=0.30)
    ap.add_argument("--min-genomes", type=int, default=4)
    a = ap.parse_args()
    z = np.load(a.pon, allow_pickle=False)
    sex = z["sex"] if "sex" in z.files else None
    tot = msk = 0
    with open(a.out, "w") as fh:
        for c in [k for k in z.files if k.startswith("chr")]:
            M = z[c].astype(float)
            if c == "chrY" and sex is not None:
                M = M[:, sex == "male"]
            with np.errstate(all="ignore"):
                k = np.sum(np.isfinite(M) & (M > 0), axis=1)
                md = np.nanmedian(M, axis=1)
                rel = 1.4826 * np.nanmedian(np.abs(M - md[:, None]), axis=1) / md
            m = ~np.isfinite(md) | (md < a.min_median) | ~(rel <= a.max_relmad) | (k < a.min_genomes)
            tot += len(m); msk += int(m.sum())
            i = 0
            while i < len(m):
                if m[i]:
                    j = i
                    while j + 1 < len(m) and m[j + 1]:
                        j += 1
                    fh.write(f"{c}\t{i * a.bin}\t{(j + 1) * a.bin}\n")
                    i = j + 1
                else:
                    i += 1
    print(f"{a.pon}: {msk}/{tot} bins masked ({msk / tot:.3%}) -> {a.out}")


if __name__ == "__main__":
    main()
