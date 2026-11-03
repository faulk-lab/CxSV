#!/usr/bin/env python3
"""Panel_Callable of every panel genome's own CxSV loci, with and without leave-out.

Panel_Callable (classify_patient_cxsv.py) is the number of panel individuals
callable at every junction of a locus, using the release's callability masks
and the rules of workflow/lib/cxsv_callability.py (Mask.breakpoints_callable:
a `flank`-bp window on at least one side of each breakpoint inside one
CALLABLE interval; individual callable = any of its samples callable).

Running classify_patient_cxsv.py once per panel genome would load all 1,620
masks 1,619 times. Here every mask is loaded once and tested against the
junctions of ALL panel loci at once, and three counts are kept per locus:

  callable_NONE   all panel individuals (what the scorer sees without leave-out)
  callable_LOIO   the locus owner's individual removed
  callable_LOFO   the owner's whole panel family removed (families.py)

Output (--out, gzip TSV): cohort sample_id individual_id locus_id n_junctions
callable_NONE callable_LOIO callable_LOFO  (NA for loci without junctions).
"""
import argparse
import gzip
import multiprocessing as mp
import os
import sys
import time
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gp_common import load_families, read_tsv  # noqa: E402
import cxsv_callability as cc  # noqa: E402

G = {}


def read_sample(args):
    """(locus ids, junction arrays, membership) of one panel sample."""
    cohort, sample, junctions_path, master_path = args
    loci = []
    with open(master_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > 2:
                loci.append(f"{f[ix['Chr']]}:{f[ix['Start']]}-{f[ix['End']]}")
    lix = {L: k for k, L in enumerate(loci)}
    c1, p1, c2, p2, memb_l, memb_j = [], [], [], [], [], []
    with open(junctions_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            lid = f[ix["locus_id"]]
            if lid == "." or not lid:
                continue
            k = len(c1)
            c1.append(f[ix["chrom1"]]); p1.append(int(f[ix["pos1"]]))
            c2.append(f[ix["chrom2"]]); p2.append(int(f[ix["pos2"]]))
            for L in lid.split(","):
                if L in lix:
                    memb_l.append(lix[L]); memb_j.append(k)
    return cohort, sample, loci, c1, p1, c2, p2, memb_l, memb_j


def init_worker():
    pass


def individual_counts(task):
    """Callable-at-every-junction per locus for one individual (OR over its samples)."""
    ind, masks, fam_code, ind_code = task
    J = G["nJ"]
    hit = np.zeros(J, dtype=bool)
    for m in masks:
        mask = cc.Mask.from_bed(os.path.join(G["cdir"], m))
        e1 = np.zeros(J, dtype=bool)
        e2 = np.zeros(J, dtype=bool)
        for chrom, idx, pos in G["g1"]:
            e1[idx] = mask.breakpoints_callable(chrom, pos, G["flank"])
        for chrom, idx, pos in G["g2"]:
            e2[idx] = mask.breakpoints_callable(chrom, pos, G["flank"])
        hit |= e1 & e2
    allc = np.logical_and.reduceat(hit[G["memb_j"]], G["starts"])     # per locus with junctions
    none = allc.astype(np.int32)
    loio = (allc & (G["owner_ind"] != ind_code)).astype(np.int32)
    lofo = (allc & (G["owner_fam"] != fam_code)).astype(np.int32)
    return none, loio, lofo


def worker_chunk(tasks):
    n = len(G["starts"])
    acc = [np.zeros(n, dtype=np.int32) for _ in range(3)]
    for t in tasks:
        for a, b in zip(acc, individual_counts(t)):
            a += b
    return acc


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--sample-manifest", required=True, help="panel sample_manifest.tsv")
    ap.add_argument("--samples-root", required=True, help="results/panel/samples")
    ap.add_argument("--callability-dir", required=True)
    ap.add_argument("--flank", type=int, required=True)
    ap.add_argument("--families", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    t0 = time.time()
    fams = load_families(a.families)
    man = read_tsv(a.sample_manifest)
    jobs = [(r["cohort"], r["sample_id"], f"{a.samples_root}/{r['cohort']}/{r['sample_id']}/junctions/junctions.tsv",
             f"{a.samples_root}/{r['cohort']}/{r['sample_id']}/cxsv/cxsv_master_table.tsv") for r in man]
    ind_of = {(r["cohort"], r["sample_id"]): r["individual_id"] for r in man}
    with mp.Pool(a.threads) as pool:
        per = pool.map(read_sample, jobs, chunksize=8)
    print(f"read {len(per)} samples in {time.time() - t0:.0f} s", file=sys.stderr)

    ind_ids = sorted(fams)
    ind_code = {i: k for k, i in enumerate(ind_ids)}
    fam_ids = sorted({fams[i]["family_id"] for i in ind_ids})
    fam_code_of = {f: k for k, f in enumerate(fam_ids)}
    fam_code = {i: fam_code_of[fams[i]["family_id"]] for i in ind_ids}

    rows_meta, c1, p1, c2, p2, ml, mj = [], [], [], [], [], [], []
    owner_ind, owner_fam = [], []
    for cohort, sample, loci, a1, b1, a2, b2, l, j in per:
        ind = ind_of[(cohort, sample)]
        lbase, jbase = len(rows_meta), len(c1)
        for L in loci:
            rows_meta.append((cohort, sample, ind, L))
            owner_ind.append(ind_code[ind]); owner_fam.append(fam_code[ind])
        c1 += a1; p1 += b1; c2 += a2; p2 += b2
        ml += [lbase + x for x in l]; mj += [jbase + x for x in j]
    del per
    nJ, nL = len(c1), len(rows_meta)
    ml = np.asarray(ml, dtype=np.int64); mj = np.asarray(mj, dtype=np.int64)
    order = np.lexsort((mj, ml))
    ml, mj = ml[order], mj[order]
    loci_with = np.unique(ml)
    starts = np.searchsorted(ml, loci_with)
    njunc = np.diff(np.append(starts, len(ml)))
    G.update(nJ=nJ, memb_j=mj, starts=starts, cdir=a.callability_dir, flank=a.flank,
             owner_ind=np.asarray(owner_ind, dtype=np.int32)[loci_with],
             owner_fam=np.asarray(owner_fam, dtype=np.int32)[loci_with])
    for end, (cs, ps) in enumerate(((c1, p1), (c2, p2))):
        by = defaultdict(list)
        for k, ch in enumerate(cs):
            by[ch].append(k)
        pa = np.asarray(ps, dtype=np.int64)
        G[f"g{end + 1}"] = [(ch, np.asarray(idx, dtype=np.int64), pa[np.asarray(idx, dtype=np.int64)])
                            for ch, idx in by.items()]
    del c1, p1, c2, p2
    print(f"{nL} loci ({len(loci_with)} with junctions), {nJ} locus junctions, "
          f"{len(ind_ids)} individuals; setup {time.time() - t0:.0f} s", file=sys.stderr)

    samples = cc.read_samples(a.callability_dir)
    masks = defaultdict(list)
    for s in samples:
        masks[s["individual_id"]].append(s["mask"])
    missing = [i for i in ind_ids if i not in masks]
    if missing:
        sys.exit(f"individuals without a callability mask: {missing[:10]}")
    tasks = [(ind_code[i], masks[i], fam_code[i], ind_code[i]) for i in ind_ids]
    chunks = [tasks[k::a.threads * 4] for k in range(a.threads * 4)]
    tot = [np.zeros(len(loci_with), dtype=np.int32) for _ in range(3)]
    with mp.get_context("fork").Pool(a.threads) as pool:
        for k, acc in enumerate(pool.imap_unordered(worker_chunk, chunks)):
            for t_, x in zip(tot, acc):
                t_ += x
            print(f"  chunk {k + 1}/{len(chunks)} {time.time() - t0:.0f} s", file=sys.stderr)
    full = {int(L): (int(n), int(x), int(y), int(z)) for L, n, x, y, z in zip(loci_with, njunc, *tot)}
    tmp = a.out + ".tmp"
    with gzip.open(tmp, "wt") as fh:
        fh.write("cohort\tsample_id\tindividual_id\tlocus_id\tn_junctions\tcallable_NONE\tcallable_LOIO\tcallable_LOFO\n")
        for k, (cohort, sample, ind, L) in enumerate(rows_meta):
            v = full.get(k)
            vals = [str(x) for x in v] if v else ["0", "NA", "NA", "NA"]
            fh.write("\t".join([cohort, sample, ind, L] + vals) + "\n")
    os.replace(tmp, a.out)
    print(f"done {time.time() - t0:.0f} s", file=sys.stderr)


if __name__ == "__main__":
    main()
