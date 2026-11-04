#!/usr/bin/env python3
"""Genome span inside FILTERED_oversized_cluster clusters (classify_cxsv: cluster > max_cluster_size
is dropped whole, so any CxSV inside it is not called), per panel genome and for given extra samples.
Also: for each positive-control truth event span, the number of panel genomes in which it lies inside
an oversized cluster (a planted event there would be lost in that genome).
Usage: oversized_survey.py OUT.tsv EVENTS_OUT.tsv [label=summary.tsv ...]"""
import csv, glob, os, sys
from collections import defaultdict
BUILD = "/mnt/sargasso1/CxSV/build/panel-v1.0/results/panel/samples"
T = "/mnt/sargasso1/CxSV/datasets/pos_control"
events = []
for f in (f"{T}/release/truth/truth_events.bed", f"{T}/addons/E13/truth/truth_events.bed"):
    for l in open(f):
        c, s, e, n = l.split("\t")[:4]; events.append((n.strip(), c, int(s), int(e)))

def merged(iv):
    out = []
    for c, s, e in sorted(iv):
        if out and out[-1][0] == c and s <= out[-1][2]: out[-1][2] = max(out[-1][2], e)
        else: out.append([c, s, e])
    return out

def survey(path):
    iv, n = [], 0
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["CxSV_Criteria"] == "FILTERED_oversized_cluster":
                iv.append((r["Chr"], int(r["Start"]), int(r["End"]))); n += 1
    m = merged(iv)
    hit = [ev for ev, c, s, e in events if any(c == mc and s < me and e > ms for mc, ms, me in m)]
    return n, sum(e - s for _, s, e in m), max([e - s for _, s, e in m], default=0), hit

rows = []
extra = dict(a.split("=", 1) for a in sys.argv[3:])
paths = [(p.split("/")[-4], p.split("/")[-3], p) for p in sorted(glob.glob(f"{BUILD}/*/*/cxsv/cxsv_summary.tsv"))]
paths += [("extra", k, v) for k, v in extra.items()]
ev_hits = defaultdict(int)
for cohort, sample, p in paths:
    n, bp, mx, hit = survey(p)
    rows.append((cohort, sample, n, bp, mx, ",".join(hit) or "."))
    if cohort != "extra":
        for h in hit: ev_hits[h] += 1
with open(sys.argv[1], "w") as fh:
    fh.write("cohort\tsample\tn_oversized_clusters\toversized_span_bp\tlargest_bp\ttruth_events_inside\n")
    for r in rows: fh.write("\t".join(map(str, r)) + "\n")
npanel = sum(1 for r in rows if r[0] != "extra")
with open(sys.argv[2], "w") as fh:
    fh.write("event\tchrom\tstart\tend\tpanel_genomes_with_event_span_in_oversized_cluster\tpanel_genomes\n")
    for ev, c, s, e in events: fh.write(f"{ev}\t{c}\t{s}\t{e}\t{ev_hits[ev]}\t{npanel}\n")
import statistics as st
p = [r for r in rows if r[0] != "extra"]
print("panel genomes", len(p), "median span Mb", st.median(r[3] for r in p) / 1e6,
      "mean", sum(r[3] for r in p) / len(p) / 1e6, "genomes with any", sum(r[2] > 0 for r in p))
for r in rows:
    if r[0] == "extra": print(r)
