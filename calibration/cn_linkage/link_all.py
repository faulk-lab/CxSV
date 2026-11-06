#!/usr/bin/env python3
"""Candidates, junction association, structure model and augmented discovery trees
for every genome (roadmap item 15, Phase 1). Run in the pipeline image.

  link_all.py --genomes genomes.tsv --seg OUT_OF_segment_all --code CODE --fai hs1.fa.fai \
      --out OUT [--threads 8]

genomes.tsv as for segment_all.py plus column src: the genome's per-sample directory
(merged_sv.vcf, cxsv/, junctions/, callers/; panel build or patient run).

For each genome:
  - records: every merged_sv.vcf record; cxsv/single_caller_filtered.tsv ids are
    marked (svim_lowq, set aside by discovery)
  - reference genomes for the segment-level panel frequency: every role=panel genome
    of ANOTHER family (all strata; for a patient: every panel genome)
  - cnlink.candidates -> OUT/candidates/<key>.tsv
  - augment.augment -> OUT/samples/<cohort>/<sample>/{cxsv,junctions,callers}
    (callers/ hard-linked from src; the scorer reads them)
Summary tables: OUT/candidates_all.tsv, OUT/loci_all.tsv.
"""
import argparse
import csv
import os
import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import augment  # noqa: E402
import cnlink  # noqa: E402

G = {}
CAND_COLS = ["key", "candidate", "decision", "cls", "nested", "inverted", "parts", "n_segments", "segments",
             "n_records", "record_types", "n_long", "n_single_filtered", "public_counts", "span"]


def read_tsv(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader([l for l in fh if not l.startswith("#")], delimiter="\t"))


def read_segs(path):
    out = []
    for r in read_tsv(path):
        out.append(dict(chrom=r["chrom"], start=int(r["start"]), end=int(r["end"]), state=int(r["state"]),
                        ploidy=int(r["ploidy"]), cn=float(r["cn"]), mean_ratio=float(r["mean_ratio"]),
                        n_bins=int(r["n_bins"]), z=float(r["z"]), nested=r.get("nested") == "True",
                        **{k: int(float(r.get(k) or 0)) for k in ("start_mask_before", "start_mask_after",
                                                                   "end_mask_before", "end_mask_after")}))
    return out


def one(g):
    t0 = time.time()
    key, src = g["key"], g["src"]
    segs = G["segs"][key]
    others = [(k, G["fam"][k], s) for k, s in G["segs"].items()
              if G["role"][k] == "panel" and G["fam"][k] != G["fam"][key] and k != key]
    recs = cnlink.read_records(f"{src}/merged_sv.vcf")
    single = set()
    sf = f"{src}/cxsv/single_caller_filtered.tsv"
    if os.path.exists(sf):
        single = {r["ID"] for r in read_tsv(sf)}
    cands = cnlink.candidates(segs, recs, others, single_filtered=single)
    rows = []
    for c in cands:
        rows.append(dict(key=key, candidate=c["id"], decision=c["decision"], cls=c["cls"], nested=c["nested"],
                         inverted=c["inverted"], parts=";".join(f"{a}:{b}-{e}" for a, b, e in c["parts"]),
                         n_segments=len(c["segs"]),
                         segments=";".join(f"{s['chrom']}:{s['start']}-{s['end']}:{s['state']}/{s['ploidy']}" for s in c["segs"]),
                         n_records=len(c["records"]), record_types=",".join(sorted(r["svtype"] for r in c["records"])),
                         n_long=c["n_long"], n_single_filtered=c["n_single_filtered"],
                         public_counts=",".join(map(str, c["public_counts"])), span=c["span"]))
    cohort, sample = key.split("/", 1) if g["role"] == "panel" else (g["cohort"], g["sample_id"])
    out_dir = f"{G['out']}/samples/{cohort}/{sample}"
    os.makedirs(f"{out_dir}/callers", exist_ok=True)
    os.makedirs(f"{out_dir}/logs", exist_ok=True)
    for f in ("sniffles_filtered.vcf", "cutesv_filtered.vcf", "svim_filtered.vcf"):
        d = f"{out_dir}/callers/{f}"
        if not os.path.exists(d):
            try:
                os.link(f"{src}/callers/{f}", d)
            except OSError:
                import shutil
                shutil.copy(f"{src}/callers/{f}", d)
    log = f"{out_dir}/logs/augment.log"
    open(log, "w").close()
    rep = augment.augment(src, out_dir, f"{src}/merged_sv.vcf", cands, G["code"], G["fai"], log)
    for r in rep:
        r["key"] = key
    return key, rows, rep, round(time.time() - t0, 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genomes", required=True)
    ap.add_argument("--seg", required=True)
    ap.add_argument("--code", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--families", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--only", default="", help="comma-separated keys to process (default all)")
    a = ap.parse_args()
    genomes = read_tsv(a.genomes)
    fam = {r["individual_id"]: r["family_id"] for r in read_tsv(a.families)}
    G.update(out=a.out, code=a.code, fai=a.fai, segs={}, fam={}, role={})
    for g in genomes:
        G["segs"][g["key"]] = read_segs(f"{a.seg}/segments/{g['key'].replace('/', '__')}.tsv")
        G["fam"][g["key"]] = fam.get(g["individual_id"], g["individual_id"]) if g["role"] == "panel" else "patient:" + g["key"]
        G["role"][g["key"]] = g["role"]
    todo = [g for g in genomes if not a.only or g["key"] in a.only.split(",")]
    os.makedirs(f"{a.out}/candidates", exist_ok=True)
    with ProcessPoolExecutor(a.threads) as ex:
        res = list(ex.map(one, todo))
    with open(f"{a.out}/candidates_all.tsv", "w", newline="") as fc, open(f"{a.out}/loci_all.tsv", "w", newline="") as fl:
        wc = csv.DictWriter(fc, CAND_COLS, delimiter="\t", lineterminator="\n")
        wc.writeheader()
        lcols = ["key", "candidate", "locus_id", "cls", "classifier_class", "classifier_iscxsv", "n_records", "n_long",
                 "n_single_filtered", "written"]
        wl = csv.DictWriter(fl, lcols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        wl.writeheader()
        for key, rows, rep, wall in res:
            wc.writerows(rows)
            wl.writerows(rep)
            print(f"{key}\tcandidates {len(rows)}\tloci {sum(r['written'] == 'yes' for r in rep)}\twall_s {wall}", file=sys.stderr)


if __name__ == "__main__":
    main()
