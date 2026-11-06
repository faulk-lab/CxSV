#!/usr/bin/env python3
"""Copy-number segments of one sample against its stratum's panel-of-normals (rule
cn_segment, workflow/rules/cn.smk; docs/CN_LINKAGE.md).

  cn_segment.py --depth depth.regions.bed.gz --pon <stratum>.npz|<pon dir> [--par par.bed] \\
      --params cn_params.json (rule cn_params_json; its cn_params.segmentation block) \\
      --stratum S|auto [--library-qc qc/library_qc.tsv --default-stratum D] --out segments.tsv --stats stats.json

Stratum (patient runs, --pon = the release's cn/pon/): a configured --stratum must have a
PoN. With --stratum auto, the library-QC group of the sample (qc_group of
qc/library_qc.tsv, detected from the alignment; lib/cxsv_libgroup.py) is used when the
release has its PoN, else --default-stratum (platform default). stats.json records the
stratum and its source.

Sex is inferred from the sample's chrX depth ratio (cxsv_cn.infer_sex). Segments:
cxsv_cn.segment_genome (ploidy-aware HMM on 10 kb bins, 1 kb boundary refinement,
nested-level scan, masked-bin extents at each boundary).
"""
import argparse
import csv
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import cxsv_cn as cn  # noqa: E402

COLS = ["chrom", "start", "end", "state", "ploidy", "cn", "mean_ratio", "n_bins", "z", "nested",
        "start_mask_before", "start_mask_after", "end_mask_before", "end_mask_after"]


def resolve_stratum(a):
    """(stratum, source): --stratum unless 'auto'; else the library-QC group of the
    sample (qc/library_qc.tsv qc_group, detected from the alignment); else
    --default-stratum (the platform default)."""
    if a.stratum != "auto":
        return a.stratum, "configured"
    if a.library_qc and os.path.exists(a.library_qc):
        with open(a.library_qc, newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        g = (rows[0].get("qc_group") or "").strip() if rows else ""
        if g and g != "NA":
            src = (rows[0].get("qc_group_source") or "").strip()
            return g, "library_qc" + (f" ({src})" if src else "")
    if not a.default_stratum:
        sys.exit("--stratum auto: no library-QC group and no --default-stratum")
    return a.default_stratum, "platform default"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    for k in ("depth", "pon", "params", "stratum", "out", "stats"):
        ap.add_argument(f"--{k}", required=True)
    ap.add_argument("--library-qc", default=None, help="qc/library_qc.tsv: its qc_group is the stratum "
                                                         "when --stratum is auto")
    ap.add_argument("--default-stratum", default=None, help="stratum when --stratum auto finds no group")
    ap.add_argument("--par", default=None, help="PAR BED; default cn_params.par of --params")
    a = ap.parse_args()
    cnp = json.load(open(a.params))["cn_params"]
    P = {**cn.DEFAULTS, **cnp["segmentation"]}
    par = cn.read_par(a.par) if a.par else cn.par_from_params(cnp.get("par"))
    fine = cn.read_bins(a.depth)
    coarse = cn.to_coarse(fine)
    gmed = cn.genome_median(coarse)
    sex, ratio = cn.infer_sex(coarse, gmed, P)
    stratum, source = resolve_stratum(a)
    if os.path.isdir(a.pon) and source.startswith("library_qc") and \
            not os.path.exists(os.path.join(a.pon, f"{stratum}.npz")) and a.default_stratum:
        print(f"no PoN for the library-QC group {stratum!r}; using the platform default "
              f"{a.default_stratum!r}", file=sys.stderr)
        stratum, source = a.default_stratum, f"platform default (no PoN for library-QC group {stratum})"
    pon_path = os.path.join(a.pon, f"{stratum}.npz") if os.path.isdir(a.pon) else a.pon
    if not os.path.exists(pon_path):
        have = sorted(x[:-4] for x in os.listdir(a.pon) if x.endswith(".npz")) if os.path.isdir(a.pon) else []
        sys.exit(f"no PoN for stratum {stratum!r} ({source}) in {a.pon}; available: {have}. "
                 "Set the sample's cn_stratum (or qc_group).")
    pon = cn.Pon.load(pon_path)
    segs, st = cn.segment_genome(fine, sex if sex != "unknown" else "female", par, pon, P)
    with open(a.out, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for s in segs:
            vals = []
            for c in COLS:
                v = s.get(c, False if c == "nested" else 0)
                vals.append(str(round(v, 4)) if isinstance(v, float) else str(v))
            fh.write("\t".join(vals) + "\n")
    stats = dict(stratum=stratum, stratum_source=source, sex=sex, x_ratio=ratio, pon_genomes=pon.meta.get("genomes"),
                 genome_median=round(st["genome_median"], 3), sigma=round(st["sigma"], 4) if st["sigma"] == st["sigma"] else None,
                 masked_fraction=round(st["masked_fraction"], 4) if st["masked_fraction"] == st["masked_fraction"] else None,
                 segments=len(segs))
    json.dump(stats, open(a.stats, "w"), indent=1, sort_keys=True)
    print(json.dumps(stats), file=sys.stderr)


if __name__ == "__main__":
    main()
