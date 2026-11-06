#!/usr/bin/env python3
"""Panel-of-normals and segmentation of every genome (roadmap item 15, Phase 1).

  segment_all.py --genomes genomes.tsv --depth-root DEPTH --families families.tsv \
      --par par.bed --out OUT [--threads 8]

genomes.tsv: key, cohort, sample_id, individual_id, stratum, sex, flagged, role
  role "panel": a subset genome. It is segmented against its stratum's PoN with its
                own family left out (LOFO), and, unless flagged (item 14 library QC,
                calibration/library_qc/flagged_v1.1.tsv), it is a PoN member.
  role "patient": segmented against the full PoN of its stratum (PC, HG01868).
Depth: DEPTH/<key>/depth.regions.bed.gz (depth_pass.sh, 1 kb bins).

Outputs in OUT:
  segments/<key>.tsv     abnormal segments (chrom start end state ploidy cn mean_ratio n_bins z)
  segments_all.tsv       all genomes, with key / family / stratum / sex
  genome_stats.tsv       per genome: genome median depth, noise sigma, PoN size used,
                         masked fraction of autosomal 10 kb bins, segments, wall seconds
  pon/<stratum>.npz      the PoN matrices (10 kb, normalised), for tests A and B
"""
import argparse
import csv
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cnlib  # noqa: E402

SEG_COLS = ["chrom", "start", "end", "state", "ploidy", "cn", "mean_ratio", "n_bins", "z", "nested",
            "start_mask_before", "start_mask_after", "end_mask_before", "end_mask_after"]
G = {}


def read_tsv(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader([l for l in fh if not l.startswith("#")], delimiter="\t"))


def load(key):
    fine = cnlib.read_bins(f"{G['depth']}/{key}/depth.regions.bed.gz")
    coarse = cnlib.to_coarse(fine)
    return key, fine, coarse, cnlib.genome_median(coarse)


def run_one(g):
    t0 = time.time()
    key = g["key"]
    pon = G["pon"][g["stratum"]]
    fam = g["family"] if g["role"] == "panel" else None
    fine = G["fine"][key]
    segs, st = cnlib.segment_genome(fine, g["sex"], G["par"], pon, family=fam)
    with open(f"{G['out']}/segments/{key.replace('/', '__')}.tsv", "w") as fh:
        fh.write("\t".join(SEG_COLS) + "\n")
        for s in segs:
            fh.write("\t".join(str(round(s.get(c, 0), 4)) if isinstance(s.get(c), float) else str(s.get(c, False if c == "nested" else 0)) for c in SEG_COLS) + "\n")
    # masked fraction (autosomes)
    tot = msk = 0
    for c in cnlib.AUTOSOMES:
        med, rel, n = pon.reference(c, g["sex"], fam)
        m = (med < cnlib.DEFAULTS["min_pon_median"]) | ~np.isfinite(med) | (rel > cnlib.DEFAULTS["max_pon_relmad"]) | (n < cnlib.DEFAULTS["min_pon_genomes"])
        tot += len(m)
        msk += int(m.sum())
    n_pon = int(np.sum(pon.family != fam)) if fam is not None else len(pon.genomes)
    return key, segs, dict(key=key, genome_median=round(st["genome_median"], 3), sigma=round(st["sigma"], 4),
                           pon_genomes=n_pon, masked_autosomal=round(msk / tot, 4), segments=len(segs),
                           wall_s=round(time.time() - t0, 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--depth-root", required=True)
    ap.add_argument("--families", required=True)
    ap.add_argument("--par", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    os.makedirs(f"{a.out}/segments", exist_ok=True)
    os.makedirs(f"{a.out}/pon", exist_ok=True)
    fam = {r["individual_id"]: r["family_id"] for r in read_tsv(a.families)}
    genomes = read_tsv(a.genomes)
    for g in genomes:
        g["family"] = fam.get(g["individual_id"], g["individual_id"]) if g["role"] == "panel" else "patient:" + g["key"]
    G.update(depth=a.depth_root, out=a.out, par=cnlib.read_par(a.par))
    with ProcessPoolExecutor(a.threads) as ex:
        loaded = list(ex.map(load, [g["key"] for g in genomes]))
    G["fine"] = {k: f for k, f, c, m in loaded}
    coarse = {k: (c, m) for k, f, c, m in loaded}
    G["pon"] = {}
    for stratum in sorted({g["stratum"] for g in genomes}):
        members = [g for g in genomes if g["stratum"] == stratum and g["role"] == "panel" and g["flagged"] != "True"]
        pon = cnlib.PanelOfNormals([g["key"] for g in members], coarse, {g["key"]: g["sex"] for g in members},
                                   {g["key"]: g["family"] for g in members}, G["par"])
        G["pon"][stratum] = pon
        np.savez_compressed(f"{a.out}/pon/{stratum.replace('/', '_')}.npz", genomes=np.array(pon.genomes),
                            sex=pon.sex, family=pon.family, **{c: pon.mat[c] for c in cnlib.CONTIGS})
        print(f"PoN {stratum}: {len(members)} genomes ({sum(g['sex'] == 'male' for g in members)} male)", file=sys.stderr)
    with ProcessPoolExecutor(a.threads) as ex:          # fork: workers inherit G
        res = list(ex.map(run_one, genomes))
    with open(f"{a.out}/segments_all.tsv", "w") as fh, open(f"{a.out}/genome_stats.tsv", "w") as fs:
        fh.write("key\tfamily\tstratum\tsex\t" + "\t".join(SEG_COLS) + "\n")
        cols = ["key", "stratum", "sex", "role", "flagged", "genome_median", "sigma", "pon_genomes", "masked_autosomal",
                "segments", "wall_s"]
        fs.write("\t".join(cols) + "\n")
        gi = {g["key"]: g for g in genomes}
        for key, segs, st in res:
            g = gi[key]
            for s in segs:
                fh.write(f"{key}\t{g['family']}\t{g['stratum']}\t{g['sex']}\t" + "\t".join(
                    str(round(s.get(c, 0), 4)) if isinstance(s.get(c), float) else str(s.get(c, False if c == "nested" else 0)) for c in SEG_COLS) + "\n")
            st.update(stratum=g["stratum"], sex=g["sex"], role=g["role"], flagged=g["flagged"])
            fs.write("\t".join(str(st[c]) for c in cols) + "\n")


if __name__ == "__main__":
    main()
