#!/usr/bin/env python3
"""Library QC: chimeric-read and single-caller-breakend metrics of one genome
(docs/LIBRARY_QC.md). Independent of CxSV discovery and scoring; NOT part of the
calling fingerprint.

    library_qc.py --aln ALN --ref REF --windows BED --merged-vcf merged_sv.vcf --out library_qc.tsv
                  [--sample NAME] [--min-mapq 20] [--thresholds config/library_qc.yaml --group GROUP]

Read level, over the fixed windows (config/library_qc_windows.hs1.bed: 20 x 1 Mb of
unique autosomal sequence), primary mapped reads with MAPQ >= --min-mapq that START
in a window (so no read is counted twice):
  frac_sa             fraction carrying a supplementary alignment (SA tag)
  frac_sa_interchrom  fraction whose SA tag names another chromosome
  coverage_est        aligned bases of those reads / window length
Call level, from the merged SURVIVOR VCF (SUPP_VEC order cuteSV, SVIM, Sniffles2):
  single_caller_bnd   BND/TRA records called by exactly one caller
  svim_only_bnd       BND/TRA records with SUPP_VEC 010
  single_caller_any   records of any type called by exactly one caller
  svim_only_any       records of any type with SUPP_VEC 010
  *_per_gb            divided by the sequence estimated from the windows
                      (coverage_est x length of chr1-22 + chrX, in Gb)
Flag (--thresholds config/library_qc.yaml, --group): for each flag metric of the
group, a coverage-adjusted robust z,
  z = (log10(metric + eps) - (a + b * log10(coverage_est)) - r0) / mad
(a, b: robust fit over the healthy panel genomes of the group; r0, mad: median and
1.4826 x MAD of its residuals; calibration/library_qc/calibrate.py). qc_flag is
FLAG when any z > k, PASS otherwise, NA when the group has no calibration or no
QC window lies on the alignment's reference.
"""
import argparse
import csv
import math
import statistics
import sys

METRICS = ("frac_sa", "frac_sa_interchrom", "single_caller_bnd_per_gb", "svim_only_bnd_per_gb",
           "single_caller_any_per_gb", "svim_only_any_per_gb")


def read_windows(path):
    out = []
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        out.append((f[0], int(f[1]), int(f[2])))
    return out


def read_metrics(aln, ref, windows, min_mapq):
    """Read-level metrics over the windows that lie on the alignment's contigs
    (a window off the header, e.g. on a test reference, is skipped)."""
    import pysam  # only this part needs pysam
    f = pysam.AlignmentFile(aln, reference_filename=ref)
    lengths = dict(zip(f.references, f.lengths))
    windows = [(c, s, e) for c, s, e in windows if c in lengths and e <= lengths[c]]
    n = sa = inter = sa2 = bases = 0
    lens = []
    for c, s, e in windows:
        for r in f.fetch(c, s, e):
            if r.reference_start < s or r.is_secondary or r.is_supplementary or r.is_unmapped:
                continue
            if r.mapping_quality < min_mapq:
                continue
            n += 1
            bases += r.query_alignment_length or 0
            lens.append(r.query_length or r.infer_read_length() or 0)
            if r.has_tag("SA"):
                sa += 1
                parts = [p for p in r.get_tag("SA").split(";") if p]
                sa2 += len(parts) >= 2
                inter += any(p.split(",", 1)[0] != r.reference_name for p in parts)
    span = sum(e - s for _, s, e in windows)
    return {"windows_used": len(windows), "reads_sampled": n, "frac_sa": sa / n if n else float("nan"),
            "frac_sa_2plus": sa2 / n if n else float("nan"),
            "frac_sa_interchrom": inter / n if n else float("nan"),
            "coverage_est": bases / span if span else float("nan"),
            "median_read_len": int(statistics.median(lens)) if lens else 0}


def info_get(info, key):
    for kv in info.split(";"):
        if kv.startswith(key + "="):
            return kv.split("=", 1)[1]
    return None


def call_metrics(vcf):
    n = bnd = single = svim = single_any = svim_any = 0
    with open(vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            n += 1
            info = line.split("\t", 8)[7]
            sv = info_get(info, "SUPP_VEC") or ""
            single_any += sv.count("1") == 1
            svim_any += sv == "010"
            if info_get(info, "SVTYPE") not in ("BND", "TRA"):
                continue
            bnd += 1
            single += sv.count("1") == 1
            svim += sv == "010"
    return {"merged_records": n, "bnd_records": bnd, "single_caller_bnd": single, "svim_only_bnd": svim,
            "single_caller_any": single_any, "svim_only_any": svim_any}


def genome_gb(ref):
    keep = {f"chr{i}" for i in range(1, 23)} | {"chrX"}
    return sum(int(l.split("\t")[1]) for l in open(ref + ".fai") if l.split("\t")[0] in keep) / 1e9


def zscore(value, coverage, c, eps):
    """Coverage-adjusted robust z of one metric (calibration c: a, b, r0, mad)."""
    if not (value == value and coverage and coverage > 0):
        return float("nan")
    return (math.log10(value + eps) - (c["a"] + c["b"] * math.log10(coverage)) - c["r0"]) / c["mad"]


def flag(row, cal, k, eps):
    """(qc_flag, qc_reason, {metric: z}) for a row and its group's calibration."""
    if not cal:
        return "NA", "no calibration for this group", {}
    zs = {m: zscore(row[m], row["coverage_est"], c, eps) for m, c in cal.items()}
    if all(z != z for z in zs.values()):
        return "NA", "not assessable (no coverage estimate from the QC windows)", zs
    over = [f"{m} {row[m]:.3g} (z {z:.1f} > {k:g})" for m, z in zs.items() if z == z and z > k]
    return ("FLAG", "; ".join(over), zs) if over else ("PASS", f"every z <= {k:g}", zs)


def load_thresholds(path):
    import yaml
    return yaml.safe_load(open(path))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("aln", "ref", "windows", "merged-vcf", "out"):
        ap.add_argument(f"--{k}", required=True)
    ap.add_argument("--sample", default="sample")
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--thresholds", default=None)
    ap.add_argument("--group", default=None)
    a = ap.parse_args()
    row = {"sample": a.sample}
    row.update(read_metrics(a.aln, a.ref, read_windows(a.windows), a.min_mapq))
    row.update(call_metrics(a.merged_vcf))
    gb = row["coverage_est"] * genome_gb(a.ref) if row["coverage_est"] == row["coverage_est"] else 0.0
    row["sequence_gb_est"] = gb
    for m in ("single_caller_bnd", "svim_only_bnd", "single_caller_any", "svim_only_any"):
        row[f"{m}_per_gb"] = row[m] / gb if gb else float("nan")
    if a.thresholds:
        t = load_thresholds(a.thresholds)
        cal = (t.get("groups") or {}).get(a.group or "", {})
        row["qc_group"] = a.group or "NA"
        row["qc_flag"], row["qc_reason"], zs = flag(row, cal, float(t["k"]), float(t["eps"]))
        for m in cal:
            row[f"z_{m}"] = zs[m]
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(list(row))
        w.writerow([f"{v:.6g}" if isinstance(v, float) else v for v in row.values()])
    print("\t".join(f"{k}={v}" for k, v in row.items()), file=sys.stderr)


if __name__ == "__main__":
    main()
