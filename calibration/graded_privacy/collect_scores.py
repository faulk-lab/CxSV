#!/usr/bin/env python3
"""One table of every scored panel locus (item 10 calibration).

Joins, per panel sample and leave-out mode, the scorer output
(persample/<c>/<s>/scores<TAG>_<MODE>.tsv), the leave-out status
(leaveout/<c>/<s>/status_<MODE>.tsv: junctions, junctions seen, matched
fraction, best individual fraction) and the catalogue hits (cats.tsv: per
catalogue k/n, max AF, contains), with sample metadata (platform, chemistry,
coverage, superpopulation, family).

nonprivacy_points = Total_Score - panel_privacy_points - catalogue_frequency_points
is the part of the score that the privacy rule does not touch; sweep.py
re-adds candidate privacy points to it.
"""
import argparse
import csv
import gzip
import multiprocessing as mp
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gp_common import load_families, read_tsv  # noqa: E402

CATS = ("CoLoRSdb", "HGSVC3", "HPRC", "gnomAD-SV")
COLS = (["mode", "cohort", "sample_id", "individual_id", "family_id", "n_panel_family", "platform", "chemistry",
         "coverage_x", "superpopulation", "Locus_ID", "Size", "Cmplx_Class", "N_Breakends", "Panel_Status",
         "Panel_Junctions", "Panel_Junctions_Seen", "Panel_Matched_Fraction", "Panel_Best_Junction_Fraction",
         "Panel_Carriers", "Panel_Callable", "Panel_Individuals", "Leaveout_Excluded", "Total_Score", "Tier",
         "panel_privacy_points", "catalogue_frequency_points", "nonprivacy_points", "gene_dosage_points",
         "no_genes_points", "Top_Features"]
        + [f"{c}_{k}" for c in CATS for k in ("k", "n", "Max_AF", "Contains")] + ["Catalogue_Max_AF"])
G = {}


def one(task):
    cohort, sample, meta = task
    a = G["a"]
    P = f"{a.cal_root}/persample/{cohort}/{sample}"
    L = f"{a.cal_root}/leaveout/{cohort}/{sample}"
    cats = {r["Locus_ID"]: r for r in read_tsv(f"{P}/cats.tsv")}
    out = []
    for mode in a.modes:
        st = {f"{r['Chr']}:{r['Start']}-{r['End']}": r for r in read_tsv(f"{L}/status_{mode}.tsv")}
        for r in read_tsv(f"{P}/scores{a.tag}_{mode}.tsv"):
            s, c = st[r["Locus_ID"]], cats.get(r["Locus_ID"], {})
            pp, cp = float(r["panel_privacy_points"]), float(r["catalogue_frequency_points"])
            row = dict(meta, mode=mode, Locus_ID=r["Locus_ID"], Size=r["Size"], Cmplx_Class=r["Cmplx_Class"],
                       N_Breakends=r["N_Breakends"], Panel_Status=r["Panel_Status"],
                       Panel_Junctions=s["Panel_Junctions"], Panel_Junctions_Seen=s["Panel_Junctions_Seen"],
                       Panel_Matched_Fraction=s["Panel_Matched_Fraction"],
                       Panel_Best_Junction_Fraction=s["Panel_Best_Junction_Fraction"],
                       Panel_Carriers=s["Panel_Carriers"], Panel_Callable=s["Panel_Callable"],
                       Panel_Individuals=s["Panel_Individuals"], Leaveout_Excluded=s["Leaveout_Excluded"],
                       Total_Score=r["Total_Score"], Tier=r["Tier"], panel_privacy_points=pp,
                       catalogue_frequency_points=cp,
                       nonprivacy_points=round(float(r["Total_Score"]) - pp - cp, 4),
                       gene_dosage_points=r["gene_dosage_points"], no_genes_points=r["no_genes_points"],
                       Top_Features=r["Top_Features"], Catalogue_Max_AF=c.get("Catalogue_Max_AF", "."))
            for cn in CATS:
                k, n = (c.get(f"{cn}_Junctions_Matched") or "0/0").split("/")
                row.update({f"{cn}_k": k, f"{cn}_n": n, f"{cn}_Max_AF": c.get(f"{cn}_Max_AF", "."),
                            f"{cn}_Contains": c.get(f"{cn}_Contains", "no")})
            out.append(row)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cal-root", required=True)
    ap.add_argument("--sample-manifest", required=True)
    ap.add_argument("--registry", required=True, help="registry sample_registry.tsv (chemistry)")
    ap.add_argument("--modes", default="LOIO,LOFO")
    ap.add_argument("--tag", default="")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    a.modes = a.modes.split(",")
    G["a"] = a
    fams = load_families(f"{a.cal_root}/families.tsv")
    chem = {}
    for r in read_tsv(a.registry):
        if r.get("status") == "verified":
            chem[(r["cohort"], r["individual_id"])] = r.get("chemistry") or ""
    tasks = []
    for r in read_tsv(a.sample_manifest):
        f = fams[r["individual_id"]]
        tasks.append((r["cohort"], r["sample_id"], {
            "cohort": r["cohort"], "sample_id": r["sample_id"], "individual_id": r["individual_id"],
            "family_id": f["family_id"], "n_panel_family": f["n_panel_family"], "platform": r["platform"],
            "chemistry": chem.get((r["cohort"], r["individual_id"]), ""), "coverage_x": r["coverage_x"],
            "superpopulation": f["superpopulation"] or "unlabelled"}))
    n = 0
    tmp = a.out + ".tmp"
    with gzip.open(tmp, "wt") as fh, mp.Pool(a.threads) as pool:
        w = csv.DictWriter(fh, COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for rows in pool.imap(one, tasks, chunksize=4):
            for r in rows:
                w.writerow(r)
            n += len(rows)
    os.replace(tmp, a.out)
    print(f"{len(tasks)} samples, {n} locus rows -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
