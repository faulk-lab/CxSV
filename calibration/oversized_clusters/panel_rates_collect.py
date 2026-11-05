#!/usr/bin/env python3
"""Collect healthy-genome Tier rates: panel-v1.0 (item-10 calibration, LOFO,
graded privacy) vs a rebuilt panel (panel_rates.sh), per genome and per stratum.

    panel_rates_collect.py --old-scores OLD_CAL/loci_scores_graded.tsv.gz \
        --manifest NEW/results/panel/panel/sample_manifest.tsv --cal CAL_ROOT --out-prefix PREFIX

Strata as in the item-10 report (graded_privacy_results/strata_before_after.tsv):
all, platform (ONT 1kGP-1019, ONT LRSC R10, ONT LRSC R9, PacBio HiFi), coverage
bins. Writes PREFIX_per_genome.tsv, PREFIX_strata.tsv, PREFIX_tier1_loci.tsv.
The scorer's points did not change between the item-10 code (dc7f582) and the
rebuild commit (65ef970 adds literature annotation columns only), so the
item-10 LOFO scores are the panel-v1.0 baseline.
"""
import argparse
import csv
import gzip
import os
from collections import defaultdict

T1, T2 = "Tier1_high_priority", "Tier2"
COV_BINS = [(0, 5, "coverage <5x"), (5, 10, "coverage 5-10x"), (10, 20, "coverage 10-20x"),
            (20, 30, "coverage 20-30x"), (30, 1e9, "coverage >=30x")]


def group_of(cohort, chemistry):
    if cohort == "1kGP_ONT_1019_CHM13":
        return "ONT 1kGP-1019"
    if cohort == "1kGP_PacBio_CHM13":
        return "PacBio HiFi"
    return f"ONT LRSC {chemistry}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--old-scores", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--cal", required=True)
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()

    old = defaultdict(lambda: [0, 0, 0])
    chem = {}
    with gzip.open(a.old_scores, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["mode"] != "LOFO":
                continue
            k = (r["cohort"], r["sample_id"])
            chem[k] = r["chemistry"]
            o = old[k]
            o[2] += 1
            o[0] += r["Tier"] == T1
            o[1] += r["Tier"] == T2
    per, t1loci = [], []
    missing = []
    for m in csv.DictReader(open(a.manifest), delimiter="\t"):
        k = (m["cohort"], m["sample_id"])
        p = f"{a.cal}/persample/{k[0]}/{k[1]}/scores_LOFO.tsv"
        if not os.path.exists(p):
            missing.append(k)
            continue
        n1 = n2 = n = 0
        for r in csv.DictReader(open(p), delimiter="\t"):
            n += 1
            n1 += r["Tier"] == T1
            n2 += r["Tier"] == T2
            if r["Tier"] == T1:
                t1loci.append({"cohort": k[0], "sample_id": k[1], "Locus_ID": r["Locus_ID"], "Size": r["Size"],
                               "Cmplx_Class": r["Cmplx_Class"], "Panel_Status": r["Panel_Status"],
                               "Total_Score": r["Total_Score"], "Top_Features": r["Top_Features"]})
        o = old[k]
        per.append({"cohort": k[0], "sample_id": k[1], "individual_id": m["individual_id"],
                    "group": group_of(k[0], chem.get(k, "NA")), "coverage_x": m["coverage_x"],
                    "loci_v1.0": o[2], "T1_v1.0": o[0], "T2_v1.0": o[1],
                    "loci_v1.1": n, "T1_v1.1": n1, "T2_v1.1": n2})
    with open(a.out_prefix + "_per_genome.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(per[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(per)
    strata = [("all", lambda r: True)]
    for g in ("ONT 1kGP-1019", "ONT LRSC R10", "ONT LRSC R9", "PacBio HiFi"):
        strata.append((f"platform: {g}", lambda r, g=g: r["group"] == g))
    for lo, hi, name in COV_BINS:
        strata.append((name, lambda r, lo=lo, hi=hi: lo <= float(r["coverage_x"]) < hi))
    cols = ["stratum", "n", "T1_per_genome_v1.0", "T1_per_genome_v1.1", "T1_ratio", "T2_per_genome_v1.0",
            "T2_per_genome_v1.1", "share_genomes_T1_v1.0", "share_genomes_T1_v1.1", "loci_per_genome_v1.0",
            "loci_per_genome_v1.1"]
    with open(a.out_prefix + "_strata.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for name, f in strata:
            rs = [r for r in per if f(r)]
            n = len(rs)
            if not n:
                continue
            s = lambda c: sum(int(r[c]) for r in rs)  # noqa: E731
            b1, a1 = s("T1_v1.0") / n, s("T1_v1.1") / n
            w.writerow([name, n, f"{b1:.4f}", f"{a1:.4f}", f"{a1 / b1:.3f}" if b1 else "NA",
                        f"{s('T2_v1.0') / n:.4f}", f"{s('T2_v1.1') / n:.4f}",
                        f"{sum(int(r['T1_v1.0']) > 0 for r in rs) / n:.4f}",
                        f"{sum(int(r['T1_v1.1']) > 0 for r in rs) / n:.4f}",
                        f"{s('loci_v1.0') / n:.1f}", f"{s('loci_v1.1') / n:.1f}"])
    with open(a.out_prefix + "_tier1_loci.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(t1loci[0]) if t1loci else ["cohort"], delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(t1loci)
    print(f"genomes {len(per)}, missing scores {len(missing)} {missing[:5]}")
    print(open(a.out_prefix + "_strata.tsv").read())


if __name__ == "__main__":
    main()
