#!/usr/bin/env python3
"""Panel-of-normals of the copy-number track, one per platform stratum (rule
panel_cn_pon, workflow/rules/cn.smk; docs/CN_LINKAGE.md).

  cn_pon.py --members members.tsv [--par par.bed] --exclude flagged.tsv --params cn_params.json (rule cn_params_json; its cn_params.segmentation block) \\
      --out-dir DIR [--threads N]

members.tsv: key (cohort/sample), individual_id, stratum, depth (1 kb regions.bed.gz).
Every member is read once: genome median of its 10 kb bins and sex from the chrX ratio
(cxsv_cn.infer_sex). Libraries listed in --exclude (cohort, sample_id; the library-QC
flag list) do not enter the PoN. Per stratum, DIR/<stratum>.npz holds the per-bin
median, relative MAD and genome count of normalised depth (cxsv_cn.Pon; chrX pools
both sexes, chrY males only). DIR/pon_members.tsv lists every library with its sex,
chrX ratio, genome median and whether it is in the PoN (and why not).
"""
import argparse
import csv
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import cxsv_cn as cn  # noqa: E402

PAR = {}
P = {}


def load(row):
    fine = cn.read_bins(row["depth"])
    coarse = cn.to_coarse(fine)
    gmed = cn.genome_median(coarse)
    sex, ratio = cn.infer_sex(coarse, gmed, P)
    # no chrX track (e.g. a single-contig test reference): every contig is diploid
    return row["key"], cn.normalised(coarse, gmed, sex if sex != "unknown" else "female", PAR), sex, ratio, gmed


def init(par, params):
    PAR.update(par)
    P.update(params)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--members", required=True)
    ap.add_argument("--par", default=None, help="PAR BED; default cn_params.par of --params")
    ap.add_argument("--exclude", default="")
    ap.add_argument("--params", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()
    cnp = json.load(open(a.params))["cn_params"]
    params = {**cn.DEFAULTS, **cnp["segmentation"]}
    par = cn.read_par(a.par) if a.par else cn.par_from_params(cnp.get("par"))
    members = list(csv.DictReader(open(a.members), delimiter="\t"))
    excl = set()
    if a.exclude and a.exclude != "none":
        for r in csv.DictReader([l for l in open(a.exclude) if not l.startswith("#")], delimiter="\t"):
            excl.add(f"{r['cohort']}/{r['sample_id']}")
    os.makedirs(a.out_dir, exist_ok=True)
    by_stratum = {}
    info = {}
    by_key = {r["key"]: r for r in members}
    if a.threads > 1:
        ex = ProcessPoolExecutor(a.threads, initializer=init, initargs=(par, params))
        results = ex.map(load, members, chunksize=4)
    else:                                     # serial (tests; hosts without POSIX semaphores)
        ex = None
        init(par, params)
        results = map(load, members)
    try:
        for key, norm, sex, ratio, gmed in results:
            info[key] = (sex, ratio, gmed)
            m = by_key[key]
            if key in excl:
                continue
            if not gmed > 0:
                continue
            by_stratum.setdefault(m["stratum"], []).append((norm, sex))
    finally:
        if ex is not None:
            ex.shutdown()
    strata = sorted({m["stratum"] for m in members})
    for s in strata:
        mem = by_stratum.get(s, [])
        pon = cn.Pon.build(mem, par, meta={"stratum": s})
        pon.save(f"{a.out_dir}/{s}.npz")
        print(f"{s}: {pon.meta['genomes']} genomes ({pon.meta['males']} male) in the PoN", file=sys.stderr)
    with open(f"{a.out_dir}/pon_members.tsv", "w") as fh:
        fh.write("key\tindividual_id\tstratum\tsex\tx_ratio\tgenome_median\tin_pon\n")
        for m in members:
            sex, ratio, gmed = info[m["key"]]
            why = ("no: library-QC flag" if m["key"] in excl else
                   "no: no depth" if not gmed > 0 else "yes")
            fh.write(f"{m['key']}\t{m['individual_id']}\t{m['stratum']}\t{sex}\t{ratio}\t{round(gmed, 3)}\t{why}\n")


if __name__ == "__main__":
    main()
