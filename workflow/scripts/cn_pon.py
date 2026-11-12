#!/usr/bin/env python3
"""Panel-of-normals of the copy-number track, one per platform stratum (rule
panel_cn_pon, workflow/rules/cn.smk; docs/CN_LINKAGE.md).

  cn_pon.py --members members.tsv [--par par.bed] [--exclude flagged.tsv] --params cn_params.json (rule cn_params_json; its cn_params.segmentation block) \\
      --out-dir DIR [--threads N]

members.tsv: key (cohort/sample), individual_id, stratum, depth (1 kb regions.bed.gz).
Every member is read once: genome median of its 10 kb bins and its sex call from the
chrX and chrY ratios (cxsv_cn.call_sex: male / female / ambiguous; ploidy follows the
call, an ambiguous library gets the legacy chrX rule). Libraries in cn_params.pon_exclude (cohort/sample_id; the library-QC
flag list) or in --exclude (TSV cohort, sample_id) do not enter the PoN. Per stratum, DIR/<stratum>.npz holds the per-bin
median, relative MAD and genome count of normalised depth (cxsv_cn.Pon; autosomes from
every member; chrX one reference per sex, chrY males only; ambiguous calls enter the
autosomes only). DIR/pon_members.tsv lists every library with its sex call, ploidy sex,
chrX and chrY ratios, genome median, whether it is in the PoN (and why not), whether
it is in its sex's chrX/chrY reference, and whether its chrX/chrY segments count towards
the segment-level panel frequency (sexchrom_public: 'no: discordant libraries of <ID>: ...'
when another library of the same individual disagrees on sex; cxsv_cnlink.
sexchrom_discordant, cn_params.linkage sexchrom_discordance_x / _y; panel-v1.3 D2). That
column does not change the PoN.
"""
import argparse
import csv
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import cxsv_cn as cn  # noqa: E402
import cxsv_cnlink as cl  # noqa: E402

PAR = {}
P = {}


def load(row):
    fine = cn.read_bins(row["depth"])
    coarse = cn.to_coarse(fine)
    gmed = cn.genome_median(coarse)
    sx = cn.call_sex(coarse, gmed, PAR, P)
    # no chrX track (e.g. a single-contig test reference): every contig is diploid (ploidy_sex female)
    return row["key"], cn.normalised(coarse, gmed, sx["ploidy_sex"], PAR), sx, gmed


def init(par, params):
    PAR.update(par)
    P.update(params)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--members", required=True)
    ap.add_argument("--par", default=None, help="PAR BED; default cn_params.par of --params")
    ap.add_argument("--exclude", default="", help="TSV cohort, sample_id; default cn_params.pon_exclude of --params")
    ap.add_argument("--params", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()
    cnp = json.load(open(a.params))["cn_params"]
    params = {**cn.DEFAULTS, **cnp["segmentation"]}
    par = cn.read_par(a.par) if a.par else cn.par_from_params(cnp.get("par"))
    members = list(csv.DictReader(open(a.members), delimiter="\t"))
    excl = set(cnp.get("pon_exclude") or [])
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
        for key, norm, sx, gmed in results:
            info[key] = (sx, gmed)
            m = by_key[key]
            if key in excl:
                continue
            if not gmed > 0:
                continue
            by_stratum.setdefault(m["stratum"], []).append((norm, sx["sex_call"]))
    finally:
        if ex is not None:
            ex.shutdown()
    strata = sorted({m["stratum"] for m in members})
    for s in strata:
        mem = by_stratum.get(s, [])
        pon = cn.Pon.build(mem, par, meta={"stratum": s})
        pon.save(f"{a.out_dir}/{s}.npz")
        print(f"{s}: {pon.meta['genomes']} genomes in the PoN ({pon.meta['males']} male, "
              f"{pon.meta['females']} female, {pon.meta['sex_ambiguous']} ambiguous: autosomes only)", file=sys.stderr)
    link = {**cl.CN_DEFAULTS, **(cnp.get("linkage") or {})}
    disc = cl.sexchrom_discordant([dict(key=m["key"], individual_id=m["individual_id"], **info[m["key"]][0])
                                   for m in members], link)
    if disc:
        print(f"{len(disc)} libraries left out of chrX/chrY public counting (sex-discordant): "
              + ", ".join(sorted(disc)), file=sys.stderr)
    with open(f"{a.out_dir}/pon_members.tsv", "w") as fh:
        fh.write("key\tindividual_id\tstratum\tsex\tsex_call\tx_ratio\ty_ratio\tgenome_median\tin_pon\t"
                 "in_sex_reference\tsexchrom_public\n")
        for m in members:
            sx, gmed = info[m["key"]]
            why = ("no: library-QC flag" if m["key"] in excl else
                   "no: no depth" if not gmed > 0 else "yes")
            sexref = "yes" if why == "yes" and sx["sex_call"] in ("male", "female") else \
                     ("no: sex " + sx["sex_call"] if why == "yes" else "no")
            fh.write(f"{m['key']}\t{m['individual_id']}\t{m['stratum']}\t{sx['ploidy_sex']}\t{sx['sex_call']}\t"
                     f"{sx['x_ratio']}\t{sx['y_ratio']}\t{round(gmed, 3)}\t{why}\t{sexref}\t"
                     f"{('no: ' + disc[m['key']]) if m['key'] in disc else 'yes'}\n")


if __name__ == "__main__":
    main()
