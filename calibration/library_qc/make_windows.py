#!/usr/bin/env python3
"""Fixed read-sampling windows for library QC (docs/LIBRARY_QC.md).

    make_windows.py --fai hs1.fa.fai --segdups segdup.final.bed --repeats repeats.final.bed \
        --out config/library_qc_windows.hs1.bed [--n 20] [--size 1000000]

One window per autosome chr1..chr20 (chr21/chr22 are short acrocentrics), size
--size, in unique sequence:
  - < 0.5% of bases in segmental duplications;
  - < 3% Satellite / Simple_repeat / Low_complexity bases;
  - >= 10 Mb from either chromosome end.
Among the eligible 1-Mb bins of a chromosome, the one nearest 30% of its length is
taken (deterministic; no random seed).
"""
import argparse
from collections import defaultdict

LOWC = {"Satellite", "Simple_repeat", "Low_complexity"}


def covered(rows, s, e):
    """Bases of [s, e) covered by the sorted, possibly overlapping intervals rows."""
    tot, cur = 0, s
    for a, b in rows:
        if b <= cur:
            continue
        if a >= e:
            break
        a2, b2 = max(a, cur), min(b, e)
        if b2 > a2:
            tot += b2 - a2
            cur = b2
    return tot


def load(path, keep=None):
    out = defaultdict(list)
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if keep is not None and (len(f) < 4 or f[3] not in keep):
            continue
        out[f[0]].append((int(f[1]), int(f[2])))
    for c in out:
        out[c].sort()
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--segdups", required=True)
    ap.add_argument("--repeats", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--size", type=int, default=1_000_000)
    ap.add_argument("--end-margin", type=int, default=10_000_000)
    ap.add_argument("--max-lowc", type=float, default=0.03)
    ap.add_argument("--max-segdup", type=float, default=0.005)
    a = ap.parse_args()
    lens = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(a.fai)}
    sd, lc = load(a.segdups), load(a.repeats, LOWC)
    rows = []
    for i in range(1, 21):
        c = f"chr{i}"
        L = lens[c]
        best = None
        for s in range(a.end_margin, L - a.end_margin - a.size + 1, a.size):
            e = s + a.size
            if covered(sd.get(c, []), s, e) / a.size >= a.max_segdup:
                continue
            f_lc = covered(lc.get(c, []), s, e) / a.size
            if f_lc >= a.max_lowc:
                continue
            d = abs((s + e) / 2 - 0.3 * L)
            if best is None or d < best[0]:
                best = (d, s, e, f_lc)
        if best is None:
            raise SystemExit(f"no eligible window on {c}")
        rows.append((c, best[1], best[2], f"{c}_qc lowc={best[3]:.4f}"))
    with open(a.out, "w") as fh:
        fh.write("# library QC read-sampling windows (calibration/library_qc/make_windows.py); hs1\n")
        for c, s, e, n in rows:
            fh.write(f"{c}\t{s}\t{e}\t{n}\n")
    print(f"{len(rows)} windows -> {a.out}")


if __name__ == "__main__":
    main()
