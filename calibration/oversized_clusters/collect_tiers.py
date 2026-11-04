#!/usr/bin/env python3
"""Per-sample Tier counts of the eval_scores.sh arms, plus every Tier-1/Tier-2 locus with
the length of its longest junction (a junction > 500 kb joins a locus only under the
breakends rule) and whether its id existed before item 13.
Usage: collect_tiers.py SUBSET.tsv OUT_PREFIX ARM=CAL_ROOT:SAMPLES_ROOT ..."""
import csv, os, sys
from collections import Counter

sub, pre = sys.argv[1], sys.argv[2]
arms = [a.split("=", 1) for a in sys.argv[3:]]
samples = [(r["cohort"], r["sample_id"], r["group"], r["coverage_x"]) for r in csv.DictReader(open(sub), delimiter="\t")]


def longest(sroot, c, s):
    out = {}
    with open(f"{sroot}/{c}/{s}/junctions/junctions.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["locus_id"] == ".":
                continue
            L = abs(int(r["pos2"]) - int(r["pos1"])) if r["chrom1"] == r["chrom2"] else 0
            for loc in r["locus_id"].split(","):
                out[loc] = max(out.get(loc, 0), L)
    return out


with open(pre + "_per_sample.tsv", "w") as fs, open(pre + "_tier12_loci.tsv", "w") as fl:
    fs.write("arm\tcohort\tsample\tgroup\tcoverage_x\tscored\tloci\ttier1\ttier2\ttier3\tbenign\n")
    fl.write("arm\tcohort\tsample\tgroup\tlocus\ttier\ttotal\tclass\tpanel_status\tlongest_junction\ttop_features\n")
    for arm, spec in arms:
        cal, sroot = spec.split(":")
        for c, s, g, cov in samples:
            p = f"{cal}/persample/{c}/{s}/scores_LOFO.tsv"
            if not os.path.exists(p):
                fs.write(f"{arm}\t{c}\t{s}\t{g}\t{cov}\t0\tNA\tNA\tNA\tNA\tNA\n")
                continue
            rows = list(csv.DictReader(open(p), delimiter="\t"))
            t = Counter(r["Tier"] for r in rows)
            fs.write(f"{arm}\t{c}\t{s}\t{g}\t{cov}\t1\t{len(rows)}\t{t['Tier1_high_priority']}\t{t['Tier2']}\t"
                     f"{t['Tier3']}\t{t['Benign_likely']}\n")
            lj = longest(sroot, c, s)
            for r in rows:
                if r["Tier"] in ("Tier1_high_priority", "Tier2"):
                    fl.write("\t".join([arm, c, s, g, r["Locus_ID"], r["Tier"], r["Total_Score"], r["Cmplx_Class"],
                                        r["Panel_Status"], str(lj.get(r["Locus_ID"], 0)), r["Top_Features"]]) + "\n")
