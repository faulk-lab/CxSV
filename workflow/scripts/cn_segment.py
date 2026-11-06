#!/usr/bin/env python3
"""Copy-number segments of one sample against its stratum's panel-of-normals (rule
cn_segment, workflow/rules/cn.smk; docs/CN_LINKAGE.md).

  cn_segment.py --depth depth.regions.bed.gz --pon <stratum>.npz [--par par.bed] \\
      --params cn_params.json (rule cn_params_json; its cn_params.segmentation block) --stratum S --out segments.tsv --stats stats.json

Sex is inferred from the sample's chrX depth ratio (cxsv_cn.infer_sex). Segments:
cxsv_cn.segment_genome (ploidy-aware HMM on 10 kb bins, 1 kb boundary refinement,
nested-level scan, masked-bin extents at each boundary).
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import cxsv_cn as cn  # noqa: E402

COLS = ["chrom", "start", "end", "state", "ploidy", "cn", "mean_ratio", "n_bins", "z", "nested",
        "start_mask_before", "start_mask_after", "end_mask_before", "end_mask_after"]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    for k in ("depth", "pon", "params", "stratum", "out", "stats"):
        ap.add_argument(f"--{k}", required=True)
    ap.add_argument("--par", default=None, help="PAR BED; default cn_params.par of --params")
    a = ap.parse_args()
    cnp = json.load(open(a.params))["cn_params"]
    P = {**cn.DEFAULTS, **cnp["segmentation"]}
    par = cn.read_par(a.par) if a.par else cn.par_from_params(cnp.get("par"))
    fine = cn.read_bins(a.depth)
    coarse = cn.to_coarse(fine)
    gmed = cn.genome_median(coarse)
    sex, ratio = cn.infer_sex(coarse, gmed, P)
    pon = cn.Pon.load(a.pon)
    segs, st = cn.segment_genome(fine, sex if sex != "unknown" else "female", par, pon, P)
    with open(a.out, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for s in segs:
            vals = []
            for c in COLS:
                v = s.get(c, False if c == "nested" else 0)
                vals.append(str(round(v, 4)) if isinstance(v, float) else str(v))
            fh.write("\t".join(vals) + "\n")
    stats = dict(stratum=a.stratum, sex=sex, x_ratio=ratio, pon_genomes=pon.meta.get("genomes"),
                 genome_median=round(st["genome_median"], 3), sigma=round(st["sigma"], 4) if st["sigma"] == st["sigma"] else None,
                 masked_fraction=round(st["masked_fraction"], 4) if st["masked_fraction"] == st["masked_fraction"] else None,
                 segments=len(segs))
    json.dump(stats, open(a.stats, "w"), indent=1, sort_keys=True)
    print(json.dumps(stats), file=sys.stderr)


if __name__ == "__main__":
    main()
