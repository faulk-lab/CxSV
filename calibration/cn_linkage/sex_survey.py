#!/usr/bin/env python3
"""Sex-chromosome depth survey of the panel (roadmap item 15, panel-v1.2 sex-stratified
chrX reference; docs/CN_LINKAGE.md "Sex"). Calibration only, not part of a run.

For every library of a members table (key, individual_id, stratum, depth), from its 1 kb
depth track (cn_depth_v1): the genome median of the 10 kb autosomal bins and, relative to
it, the chrX and chrY depth outside the PARs:

  x_ratio        median chrX non-PAR 10 kb bin / genome median (female ~1, male ~0.5)
  y_ratio        median chrY non-PAR 10 kb bin, over bins with depth > 0, / genome median
  y_frac         fraction of chrY non-PAR 10 kb bins with depth >= 0.1 x genome median
  y_mean         mean chrY non-PAR depth / genome median

plus the registry sex (--individuals: the release's panel_individuals.tsv, columns
individual_id, sex; from the sample registry).
With --classify, the pipeline's sex call (workflow/lib/cxsv_cn.call_sex, thresholds from
cn_params.segmentation of --par-yaml; its chrY ratio is y_mean) is added too, so the
thresholds can be checked against the same table (columns sex_call, ploidy_sex).

  sex_survey.py --members pon_input.tsv --individuals panel_individuals.tsv --par-yaml config/calling.yaml \
      --out sex_survey.tsv [--threads 12] [--depth-root /mnt/sargasso1/CxSV/build/cn_depth_v1]
"""
import argparse
import csv
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "workflow", "lib"))
import cxsv_cn as cn  # noqa: E402

PAR = {}
P = {}


def nonpar_mask(contig, nbins, par, binsize=cn.COARSE * 1000):
    m = np.ones(nbins, bool)
    for s, e in par.get(contig, []):
        m[s // binsize:int(np.ceil(e / binsize))] = False
    return m


def survey(row):
    try:
        fine = cn.read_bins(row["depth"])
    except Exception as e:                       # unreadable track: reported, not dropped
        return dict(key=row["key"], error=str(e)[:200])
    coarse = cn.to_coarse(fine)
    gmed = cn.genome_median(coarse)
    out = dict(key=row["key"], genome_median=round(gmed, 3))
    for c, tag in (("chrX", "x"), ("chrY", "y")):
        a = coarse.get(c)
        if a is None or not gmed > 0:
            out.update({f"{tag}_ratio": float("nan")})
            continue
        a = a[nonpar_mask(c, len(a), PAR)]
        a = a[np.isfinite(a)]
        if tag == "x":
            v = a[a > 0]
            out["x_ratio"] = round(float(np.median(v)) / gmed, 4) if len(v) else float("nan")
        else:
            v = a[a > 0]
            out["y_ratio"] = round(float(np.median(v)) / gmed, 4) if len(v) else 0.0
            out["y_frac"] = round(float(np.mean(a >= 0.1 * gmed)), 4) if len(a) else 0.0
            out["y_mean"] = round(float(np.mean(a)) / gmed, 4) if len(a) else 0.0
    if P:          # the pipeline's own call (its chrY ratio is the non-PAR mean, column y_mean here)
        sx = cn.call_sex(coarse, gmed, PAR, P)
        out["sex_call"], out["ploidy_sex"] = sx["sex_call"], sx["ploidy_sex"]
    return out


def init(par, params):
    PAR.update(par)
    P.update(params)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--members", required=True)
    ap.add_argument("--individuals", default=None, help="panel_individuals.tsv (individual_id, sex)")
    ap.add_argument("--par-yaml", required=True, help="config/calling.yaml (cn_params.par)")
    ap.add_argument("--classify", action="store_true", help="apply cxsv_cn.classify_sex with cn_params.segmentation")
    ap.add_argument("--depth-root", default=None, help="read <root>/<key>/depth.regions.bed.gz instead of the depth column")
    ap.add_argument("--extra", nargs="*", default=[], help="key=path of extra depth tracks (e.g. the positive control)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    cnp = yaml.safe_load(open(a.par_yaml))["cn_params"]
    par = cn.par_from_params(cnp.get("par"))
    params = {**cn.DEFAULTS, **cnp["segmentation"]} if a.classify else {}
    rows = list(csv.DictReader(open(a.members), delimiter="\t"))
    if a.depth_root:
        for r in rows:
            r["depth"] = f"{a.depth_root}/{r['key']}/depth.regions.bed.gz"
    for x in a.extra:
        k, p = x.split("=", 1)
        rows.append(dict(key=k, individual_id=k.split("/")[-1], stratum="extra", depth=p))
    reg = {}
    if a.individuals:
        for r in csv.DictReader(open(a.individuals), delimiter="\t"):
            reg[r["individual_id"]] = (r.get("sex") or "").strip() or "unknown"
    by = {r["key"]: r for r in rows}
    with ProcessPoolExecutor(a.threads, initializer=init, initargs=(par, params)) as ex:
        res = list(ex.map(survey, rows, chunksize=2))
    cols = ["key", "individual_id", "stratum", "registry_sex", "genome_median", "x_ratio", "y_ratio",
            "y_frac", "y_mean"] + (["sex_call", "ploidy_sex"] if params else []) + ["error"]
    with open(a.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for o in res:
            r = by[o["key"]]
            o.update(individual_id=r["individual_id"], stratum=r["stratum"], registry_sex=reg.get(r["individual_id"], "unknown"))
            fh.write("\t".join(str(o.get(c, "")) for c in cols) + "\n")
    print(f"{len(res)} libraries -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
