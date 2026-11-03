#!/usr/bin/env python3
"""Chance-match null for the graded privacy penalty (item 10).

A graded penalty also hits a truly private event whose junctions happen to
match panel sites by coincidence. To measure how often that happens, every
panel locus is moved, unchanged in structure, to an unrelated place: all of
its junction coordinates are shifted by one random offset d (|d| uniform in
[--min-shift, --max-shift] on the same chromosome(s), sign chosen so every end
stays on its contig) and the shifted junctions are matched against the panel
with the owner's family removed (LOFO; same matcher, tolerance and INS rule as
the scorer). The matched fraction of a shifted locus is what an unrelated
event of that shape would show by chance.

Output (--out, gzip TSV): cohort sample_id Locus_ID n_junctions shift
null_seen null_matched_fraction null_best_fraction null_status
"""
import argparse
import gc
import gzip
import multiprocessing as mp
import os
import random
import sys
import time
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gp_common import ExcludingIndex, excluded, load_families, read_tsv  # noqa: E402
import cxsv_junctions as cj  # noqa: E402
import cxsv_matching as cm   # noqa: E402

G = {}


def one(task):
    cohort, sample, individual = task
    a = G["a"]
    S = f"{a.samples_root}/{cohort}/{sample}"
    excl = excluded(individual, G["families"], "LOFO") & G["panel"]
    idx = ExcludingIndex(G["index"], excl)
    by_locus = defaultdict(list)
    for j in cj.read_table(f"{S}/junctions/junctions.tsv"):
        for L in j["loci"]:
            by_locus[L].append(j)
    out = []
    for L in sorted(by_locus):
        js = by_locus[L]
        rng = random.Random(f"{a.seed}:{cohort}:{sample}:{L}")
        d = None
        for _ in range(20):
            x = rng.randint(a.min_shift, a.max_shift)
            for cand in ((x, -x) if rng.random() < 0.5 else (-x, x)):
                if all(1 <= j["pos1"] + cand <= G["len"].get(j["chrom1"], 0) and
                       1 <= j["pos2"] + cand <= G["len"].get(j["chrom2"], 0) for j in js):
                    d = cand
                    break
            if d is not None:
                break
        if d is None:
            out.append((cohort, sample, L, len(js), "NA", "NA", "NA", "NA", "unshiftable"))
            continue
        sj = [dict(j, pos1=j["pos1"] + d, pos2=j["pos2"] + d) for j in js]
        res = cm.classify_locus(sj, idx, a.tolerance, a.ins_len_ratio, a.min_junction_fraction, a.min_carriers)
        out.append((cohort, sample, L, res["n_junctions"], d, res["n_seen"],
                    f"{res['n_seen'] / res['n_junctions']:.4f}", f"{res['best_fraction']:.3f}", res["status"]))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--release-dir", required=True)
    ap.add_argument("--sample-manifest", required=True)
    ap.add_argument("--samples-root", required=True)
    ap.add_argument("--families", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--min-shift", type=int, default=1_000_000)
    ap.add_argument("--max-shift", type=int, default=5_000_000)
    ap.add_argument("--seed", default="20261103")
    ap.add_argument("--tolerance", type=int, default=500)
    ap.add_argument("--ins-len-ratio", type=float, default=0.5)
    ap.add_argument("--min-junction-fraction", type=float, default=0.75)
    ap.add_argument("--min-carriers", type=int, default=1)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    t0 = time.time()
    G["a"] = a
    G["families"] = load_families(a.families)
    G["panel"] = {r["individual_id"] for r in read_tsv(f"{a.release_dir}/panel_individuals.tsv")}
    G["index"] = cm.SiteIndex(cm.read_carriers(f"{a.release_dir}/panel_junction_carriers.tsv"))
    G["len"] = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(a.fai)}
    tasks = [(r["cohort"], r["sample_id"], r["individual_id"]) for r in read_tsv(a.sample_manifest)]
    gc.collect()
    gc.freeze()
    n = 0
    tmp = a.out + ".tmp"
    with gzip.open(tmp, "wt") as fh, mp.get_context("fork").Pool(a.threads) as pool:
        fh.write("cohort\tsample_id\tLocus_ID\tn_junctions\tshift\tnull_seen\tnull_matched_fraction\t"
                 "null_best_fraction\tnull_status\n")
        for rows in pool.imap_unordered(one, tasks, chunksize=2):
            for r in rows:
                fh.write("\t".join(map(str, r)) + "\n")
            n += 1
            if n % 100 == 0:
                print(f"{n}/{len(tasks)} samples {time.time() - t0:.0f} s", file=sys.stderr, flush=True)
    os.replace(tmp, a.out)
    print(f"done {len(tasks)} samples {time.time() - t0:.0f} s", file=sys.stderr)


if __name__ == "__main__":
    main()
