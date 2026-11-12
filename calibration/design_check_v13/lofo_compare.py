#!/usr/bin/env python3
"""Design check step 4 (subset): LOFO Tier 1 / Tier 2 per genome, v1.3 integration code vs the
panel-v1.2-r2 LOFO, on the same genomes. Reported tier = column Tier (after tier ceiling,
library-QC and depth-noise gates); Tier 1 = Tier1_high_priority.
  lofo_compare.py NEW_CAL OLD_CAL MANIFEST OUT_DIR"""
import csv, os, sys
NEW, OLD, MAN, OUT = sys.argv[1:5]
os.makedirs(OUT, exist_ok=True)
T1, T2 = "Tier1_high_priority", "Tier2"
rd = lambda p: list(csv.DictReader(open(p), delimiter="\t"))
per, loci = [], []
for m in rd(MAN):
    c, s = m["cohort"], m["sample_id"]
    row = dict(cohort=c, sample_id=s, platform=m["platform"])
    t1 = {}
    for tag, root in (("v12", OLD), ("v13", NEW)):
        sc = rd(f"{root}/persample/{c}/{s}/scores_LOFO.tsv")
        row[f"loci_{tag}"] = len(sc)
        row[f"T1_{tag}"] = sum(r["Tier"] == T1 for r in sc)
        row[f"T2_{tag}"] = sum(r["Tier"] == T2 for r in sc)
        row[f"T1_public_{tag}"] = sum(r["Tier"] == T1 and r.get("Panel_Status") == "public" for r in sc)
        row[f"T1_by_points_public_{tag}"] = sum(r.get("Score_Tier", r["Tier"]) == T1 and r.get("Panel_Status") == "public" for r in sc)
        if tag == "v13":
            row["depth_noise_gate"] = ";".join(sorted({r.get("Depth_Noise_Gate", "") for r in sc}))
            row["capped_public"] = sum(r.get("Tier_Ceiling", "not_applied") not in ("not_applied", "") for r in sc)
        t1[tag] = {f"{r['Chr']}:{r['Start']}-{r['End']}": r for r in sc if r["Tier"] == T1}
    for lid in sorted(set(t1["v12"]) | set(t1["v13"])):
        a, b = t1["v12"].get(lid), t1["v13"].get(lid)
        if a and b:
            continue
        r = b or a
        loci.append(dict(cohort=c, sample_id=s, locus=lid, change="gained" if b else "lost",
                         cls=r.get("Cmplx_Class", ""), status=r.get("Panel_Status", ""), score=r.get("Total_Score", ""),
                         criteria=r.get("CxSV_Criteria", "")))
    per.append(row)
def wr(name, rows):
    if rows:
        w = csv.DictWriter(open(f"{OUT}/{name}", "w"), fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
wr("lofo_subset_per_genome.tsv", per); wr("lofo_subset_tier1_changes.tsv", loci)
n = len(per)
summ = dict(genomes=n, **{f"{k}_per_genome": round(sum(r[k] for r in per) / n, 4) for k in ("T1_v12", "T1_v13", "T2_v12", "T2_v13")},
            T1_v12=sum(r["T1_v12"] for r in per), T1_v13=sum(r["T1_v13"] for r in per),
            T1_public_v13=sum(r["T1_public_v13"] for r in per), T1_public_v12=sum(r["T1_public_v12"] for r in per),
            T1_by_points_public_v13=sum(r["T1_by_points_public_v13"] for r in per),
            genomes_any_T1_v12=sum(r["T1_v12"] > 0 for r in per), genomes_any_T1_v13=sum(r["T1_v13"] > 0 for r in per),
            T1_gained=sum(x["change"] == "gained" for x in loci), T1_lost=sum(x["change"] == "lost" for x in loci))
wr("lofo_subset_summary.tsv", [summ])
for k, v in summ.items():
    print(k, v)
