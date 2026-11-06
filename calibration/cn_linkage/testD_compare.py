#!/usr/bin/env python3
"""Test D comparison (roadmap item 15, Phase 1): per-genome LOFO scores of the
augmented subset (cn_lofo.sh) against the panel-v1.1-r2 LOFO run on the same genomes.

  testD_compare.py NEW_CAL OLD_CAL SUBSET_MANIFEST OUT.tsv

OUT.tsv: per genome loci / Tier 1 / Tier 2 before and after, loci whose score or tier
changed (must be 0: the CN track only adds loci), new loci and new Tier 1.
OUT_loci.tsv: every added locus with class, panel status, score, tier, top features.
"""
import csv, glob, os, sys
new, old, man, out = sys.argv[1:5]
rows = []
keep = [(r["cohort"], r["sample_id"]) for r in csv.DictReader(open(man), delimiter="\t")]
def read(p):
    return list(csv.DictReader(open(p), delimiter="\t")) if os.path.exists(p) else None
with open(out, "w") as fo, open(out.replace(".tsv", "_loci.tsv"), "w") as fl:
    fo.write("cohort\tsample\tloci_old\tloci_new\ttier1_old\ttier1_new\ttier2_old\ttier2_new\tchanged_scores\tnew_loci\tnew_tier1\n")
    fl.write("cohort\tsample\tlocus\tclass\tstatus\tscore\ttier\tcriteria\ttop\n")
    for c, s in keep:
        a = read(f"{old}/persample/{c}/{s}/scores_LOFO.tsv"); b = read(f"{new}/persample/{c}/{s}/scores_LOFO.tsv")
        if a is None or b is None:
            fo.write(f"{c}\t{s}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n"); continue
        A = {r["Locus_ID"]: r for r in a}; Bm = {r["Locus_ID"]: r for r in b}
        changed = sum(1 for k in A if k in Bm and (A[k]["Total_Score"], A[k]["Tier"]) != (Bm[k]["Total_Score"], Bm[k]["Tier"]))
        newk = [k for k in Bm if k not in A]
        t1 = lambda R: sum(r["Tier"].startswith("Tier1") for r in R)
        t2 = lambda R: sum(r["Tier"] == "Tier2" for r in R)
        nt1 = sum(Bm[k]["Tier"].startswith("Tier1") for k in newk)
        fo.write(f"{c}\t{s}\t{len(a)}\t{len(b)}\t{t1(a)}\t{t1(b)}\t{t2(a)}\t{t2(b)}\t{changed}\t{len(newk)}\t{nt1}\n")
        for k in newk:
            r = Bm[k]
            fl.write("\t".join([c, s, k, r.get("Cmplx_Class", ""), r.get("Panel_Status", ""), r["Total_Score"], r["Tier"],
                                r.get("CxSV_Criteria", "")[:60], r.get("Top_Features", "")]) + "\n")
print("ok")
