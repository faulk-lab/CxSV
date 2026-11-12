#!/usr/bin/env python3
"""D5 positive-control event table: CN segments, candidates and CN loci of the PC inside each
planted event window, per swept nested_min_z. Windows: the event's non-flank segments
(events yaml; p/q flanks excluded), or addons/<id>/truth/truth_events.bed for L01-L04.
Supersedes the PC section of d5_collect.py (its windows matched nothing).
  d5_pc_events.py DIR OUT.tsv   (DIR: z<Z>/{segments,cn_loci,cn_candidates}.tsv, events yaml, <L>.truth_events.bed)"""
import csv, sys
from collections import defaultdict
import yaml
D, OUT = sys.argv[1:3]
ev = yaml.safe_load(open(f"{D}/events_E01-E13+L01-L04.yaml"))["events"]
win = {}
for e in ev:
    byc = defaultdict(list)
    for k, v in (e.get("segments") or {}).items():
        if k in ("p", "q") or k.startswith(("p_", "q_")):
            continue
        byc[v[0]] += [int(v[1]), int(v[2])]
    if not byc:
        for line in open(f"{D}/{e['id']}.truth_events.bed"):
            f = line.rstrip("\n").split("\t"); byc[f[0]] += [int(f[1]), int(f[2])]
    win[e["id"]] = {c: (min(v), max(v)) for c, v in byc.items()}
inwin = lambda w, ch, a, b: ch in w and int(b) > w[ch][0] and int(a) < w[ch][1]
rd = lambda p: list(csv.DictReader(open(p), delimiter="\t"))
rows = []
for z in ["5.0", "5.25", "5.5", "5.75", "6.0"]:
    segs, loci, cands = (rd(f"{D}/z{z}/{f}.tsv") for f in ("segments", "cn_loci", "cn_candidates"))
    for eid, w in win.items():
        ss = [s for s in segs if inwin(w, s["chrom"], s["start"], s["end"])]
        ll = [x for x in loci if inwin(w, x["locus_id"].split(":")[0], *x["locus_id"].split(":")[1].split("-"))]
        cc = [x for x in cands if any(inwin(w, p.split(":")[0], *p.split(":")[1].split("-")) for p in x["parts"].split(";") if ":" in p)]
        rows.append(dict(z=z, event=eid, n_segments=len(ss),
                         segments=";".join(f"{s['chrom']}:{s['start']}-{s['end']} {s['state']}/{s['ploidy']}{' nested' if s['nested']=='True' else ''} z={s['z']}" for s in ss),
                         candidates=";".join(f"{x['decision']}:{x['cls']}" for x in cc),
                         cn_loci=";".join(f"{x['locus_id']}={x['cn_class']}({x.get('written','')})" for x in ll)))
w = csv.DictWriter(open(OUT, "w"), fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
