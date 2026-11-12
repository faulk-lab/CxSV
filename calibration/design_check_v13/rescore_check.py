#!/usr/bin/env python3
"""Design check step 4: per-sample criteria on the v1.3 downstream re-scores vs the original runs.
  rescore_check.py OUT_ROOT TASKS.tsv OUT_DIR
Per sample: Tier 1 / Tier 2 counts (column Tier, reported) before -> after, Tier 1 loci gained /
lost, public Tier 1, malformed locus ids (End <= Start, Start < 0, or length <= 2 bp), the sex
call and depth-noise lines of the score report, the depth-noise gate state, and the two named
loci (BOC_2057 BRSK2, BOC_4433 NF1)."""
import csv, os, re, sys
ROOT, TASKS, OUT = sys.argv[1:4]
os.makedirs(OUT, exist_ok=True)
T1, T2 = "Tier1_high_priority", "Tier2"
rd = lambda p: list(csv.DictReader(open(p), delimiter="\t"))
NAMED = {"BOC_2057": ("BRSK2", "chr11", 1_465_980, 1_509_727), "BOC_4433": ("NF1", "chr17", 31_902_782, 32_067_453)}
per, changes, named = [], [], []
for line in open(TASKS):
    name, src, label, sample = line.rstrip("\n").split("\t")[:4]
    S = f"{src}/results/patient/samples/{label}/{sample}"
    old, new = rd(f"{S}/scoring/cxsv_scores.tsv"), rd(f"{ROOT}/{name}/scoring/cxsv_scores.tsv")
    lid = lambda r: f"{r['Chr']}:{r['Start']}-{r['End']}"
    bad = [lid(r) for r in new if int(r["End"]) <= int(r["Start"]) or int(r["Start"]) < 0 or int(r["End"]) - int(r["Start"]) <= 2]
    rep = open(f"{ROOT}/{name}/scoring/cxsv_score_report.md").read()
    sexl = [l.strip() for l in rep.splitlines() if re.search(r"\bsex\b", l, re.I)]
    dnl = [l.strip() for l in rep.splitlines() if "depth noise" in l.lower() or "depth-noise" in l.lower()]
    o1 = {lid(r) for r in old if r["Tier"] == T1}; n1 = {lid(r) for r in new if r["Tier"] == T1}
    gate = sorted({r.get("Depth_Noise_Gate", "") for r in new})
    per.append(dict(sample=name, loci_before=len(old), loci_after=len(new),
                    T1_before=len(o1), T1_after=len(n1), T2_before=sum(r["Tier"] == T2 for r in old),
                    T2_after=sum(r["Tier"] == T2 for r in new),
                    T1_public_after=sum(r["Tier"] == T1 and r.get("Panel_Status") == "public" for r in new),
                    capped_public_after=sum(r.get("Tier_Ceiling", "not_applied") not in ("", "not_applied") for r in new),
                    malformed_ids=";".join(bad) or "none",
                    identical_scores_file=open(f"{S}/scoring/cxsv_scores.tsv").read() == open(f"{ROOT}/{name}/scoring/cxsv_scores.tsv").read(),
                    sex_line=" | ".join(sexl[:2]) or "MISSING", depth_noise_line=" | ".join(dnl[:1]) or "MISSING",
                    depth_noise_gate=";".join(gate)))
    for l in sorted(o1 ^ n1):
        r = next(x for x in (new if l in n1 else old) if lid(x) == l)
        changes.append(dict(sample=name, locus=l, change="gained" if l in n1 else "lost", cls=r.get("Cmplx_Class", ""),
                            status=r.get("Panel_Status", ""), score=r.get("Total_Score", ""), criteria=r.get("CxSV_Criteria", "")))
    if sample in NAMED and not name.endswith("maincode"):
        g, ch, a, b = NAMED[sample]
        for tag, rows in (("before", old), ("after", new)):
            for r in rows:
                if r["Chr"] == ch and int(r["Start"]) == a and int(r["End"]) == b:
                    named.append(dict(sample=name, gene=g, locus=lid(r), version=tag, Tier=r["Tier"], Score_Tier=r.get("Score_Tier", r["Tier"]),
                                      Total_Score=r["Total_Score"], Panel_Status=r.get("Panel_Status", ""),
                                      Tier_Ceiling=r.get("Tier_Ceiling", ""), caller_agreement=r.get("caller_agreement_value", ""),
                                      caller_points=r.get("caller_agreement_points", ""), Depth_Noise_Gate=r.get("Depth_Noise_Gate", "")))
def wr(n, rows):
    if rows:
        w = csv.DictWriter(open(f"{OUT}/{n}", "w"), fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
wr("rescore_per_sample.tsv", per); wr("rescore_tier1_changes.tsv", changes); wr("rescore_named_loci.tsv", named)
for r in per:
    print({k: r[k] for k in ("sample", "T1_before", "T1_after", "T2_before", "T2_after", "T1_public_after", "capped_public_after", "malformed_ids", "identical_scores_file", "depth_noise_gate")})
    print("   sex:", r["sex_line"][:160]); print("   dn :", r["depth_noise_line"][:160])
for r in named: print(r)
for r in changes: print(r)
