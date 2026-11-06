#!/usr/bin/env python3
"""Mechanics check of the CN track on the positive control: segments taken from the
truth copy numbers (no depth), then cnlink + augment (+ patient_rescore separately).
Not a test result; it checks that a correctly segmented L01 becomes a scored locus."""
import csv, json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import augment, cnlink
run, out, code, fai = sys.argv[1:5]
truths = sys.argv[5:]
segs = []
for t in truths:
    for r in csv.DictReader(open(t), delimiter="\t"):
        cp = int(r["copies_on_derivative"])
        if cp == 1:
            continue
        pl = 2
        segs.append(dict(chrom=r["chrom"], start=int(r["start"]), end=int(r["end"]), state=pl - 1 + cp, ploidy=pl,
                         cn=pl - 1 + cp, mean_ratio=(pl - 1 + cp) / pl, n_bins=10, z=10.0))
segs = [s for s in segs if s["end"] - s["start"] >= 30000]
print(json.dumps(segs[:20]))
recs = cnlink.read_records(f"{run}/merged_sv.vcf")
cands = cnlink.candidates(segs, recs, [])
for c in cands:
    print(c["decision"], c["cls"], c["parts"], [ (r["id"], r["svtype"]) for r in c["records"]])
os.makedirs(f"{out}/callers", exist_ok=True)
for f in ("sniffles_filtered.vcf", "cutesv_filtered.vcf", "svim_filtered.vcf"):
    if not os.path.exists(f"{out}/callers/{f}"):
        os.link(f"{run}/callers/{f}", f"{out}/callers/{f}")
rep = augment.augment(run, out, f"{run}/merged_sv.vcf", cands, code, fai, f"{out}/augment.log")
print(json.dumps(rep, indent=1))
