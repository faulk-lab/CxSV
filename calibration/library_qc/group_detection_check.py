#!/usr/bin/env python3
"""Validate cxsv_libgroup (library-QC calibration group detection, docs/LIBRARY_QC.md 4a) on every panel alignment: header + up to N sampled reads in the QC windows,
compared with the calibration group of the genome's cohort/chemistry. Usage:
panel_group_check.py REGISTRY DATASETS WINDOWS LIBDIR OUT.tsv [N_READS] [JOBS]"""
import csv, sys, os
from concurrent.futures import ProcessPoolExecutor
reg, ds, win, lib, out = sys.argv[1:6]
N = int(sys.argv[6]) if len(sys.argv) > 6 else 3000
J = int(sys.argv[7]) if len(sys.argv) > 7 else 12
sys.path.insert(0, lib)
import cxsv_libgroup as g
REF = "/mnt/sargasso1/CxSV/refs/hs1.fa"
W = [(l.split()[0], int(l.split()[1]), int(l.split()[2])) for l in open(win) if l.strip() and not l.startswith("#")]

def expected(r):
    if r["cohort"] == "1kGP_ONT_1019_CHM13": return g.ONT1019
    if r["cohort"] == "1kGP_PacBio_CHM13": return g.PACBIO
    return {"R10": g.R10, "R9": g.R9}.get(r["chemistry"], "NA")

def one(r):
    import pysam
    p = os.path.join(ds, r["rel_path"])
    try:
        f = pysam.AlignmentFile(p, reference_filename=REF)
        h = g.classify_text(g.header_text(f.header.to_dict()), "header")
        c = g.ReadCues(N)
        for ch, s, e in W:
            for x in f.fetch(ch, s, e):
                if x.reference_start < s or x.is_secondary or x.is_supplementary or x.mapping_quality < 20:
                    continue
                c.add(x.query_name, x.query_qualities, x.get_tag('RG') if x.has_tag('RG') else None, x.has_tag('rq'))
                if c.n >= N: break
            if c.n >= N: break
        s_ = c.summary(); rr = g.classify_reads(s_, {})  # DEFAULTS = config/library_qc.yaml group_detection
        res = g.resolve(None, h, rr, r["platform"], {"ONT": g.ONT1019, "PACBIO": g.PACBIO})
        return [r["file_id"], r["cohort"], r["chemistry"], r["basecaller"], expected(r), res["qc_group"],
                res["qc_group_source"], h[0] or "", rr[0] or "", s_["reads"], f"{s_['median_mean_q']:.2f}",
                f"{s_['uniform_q_fraction']:.3f}", f"{s_['pacbio_fraction']:.3f}", res["qc_group_evidence"]]
    except Exception as e:
        return [r["file_id"], r["cohort"], r["chemistry"], r["basecaller"], expected(r), "ERROR", str(e)[:80], "", "", 0, "", "", "", ""]

rows = [r for r in csv.DictReader(open(reg), delimiter="\t") if r["status"] == "verified"]
with ProcessPoolExecutor(J) as ex:
    res = list(ex.map(one, rows, chunksize=4))
with open(out, "w") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["file_id", "cohort", "chemistry", "basecaller", "expected", "detected", "source", "header_group",
                "read_group", "reads", "median_mean_q", "uniform_q_fraction", "pacbio_fraction", "evidence"])
    w.writerows(res)
from collections import Counter
print(Counter((x[4], x[5], x[6]) for x in res))
