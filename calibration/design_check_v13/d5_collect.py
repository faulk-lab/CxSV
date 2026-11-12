#!/usr/bin/env python3
"""panel-v1.3 design check D5: collect the nested_min_z sweep.

Segmentation level (all 1,619 panel libraries + PC): nested segments (cn_segment.py nested-level
scan, column nested=True) and all abnormal segments per genome, by z and PoN stratum.
Linkage level (113-genome subset + PC; cn_linkage.py against the panel table of the same z):
CN candidates by decision, CN loci written, and loci whose segments include a nested segment.
Positive control: segments and CN loci in each planted event window (events_E01-E13+L01-L04.yaml).
  d5_collect.py D5 SUBSET_TASKS EVENTS_YAML OUT_DIR
"""
import csv, json, os, sys
from collections import Counter, defaultdict
import yaml

D5, TASKS, EVENTS, OUT = sys.argv[1:5]
ZS = ["5.0", "5.25", "5.5", "5.75", "6.0"]
REL = "/mnt/sargasso1/CxSV/build/panel-v1.2-r2/results/panel/release"
os.makedirs(OUT, exist_ok=True)
rd = lambda p: list(csv.DictReader([l for l in open(p) if not l.startswith("#")], delimiter="\t"))
mem = rd(f"{REL}/cn/pon/pon_members.tsv")
subset = {l.split()[0] for l in open(TASKS)}
keys = [m["key"] for m in mem] + ["extra/pos_control"]
strat = {m["key"]: m["stratum"] for m in mem}; strat["extra/pos_control"] = "ONT_LRSC_R10"

# --- segmentation level ---
per = []
for z in ZS:
    for k in keys:
        segs = rd(f"{D5}/seg/z{z}/{k.replace('/', '__')}.tsv")
        nested = [s for s in segs if s["nested"] == "True"]
        per.append(dict(z=z, key=k, stratum=strat[k], in_subset=k in subset, segments=len(segs),
                        nested_segments=len(nested),
                        nested_autosomal=sum(s["chrom"] not in ("chrX", "chrY") for s in nested)))
w = csv.DictWriter(open(f"{OUT}/d5_seg_per_genome.tsv", "w"), fieldnames=list(per[0]), delimiter="\t"); w.writeheader(); w.writerows(per)
agg = defaultdict(lambda: [0, 0, 0, 0])
for r in per:
    if r["key"] == "extra/pos_control":
        continue
    for grp in ("all", r["stratum"]) + (("subset",) if r["in_subset"] else ()):
        a = agg[(r["z"], grp)]; a[0] += 1; a[1] += r["nested_segments"]; a[2] += r["segments"]; a[3] += r["nested_segments"] > 0
rows = [dict(z=z, group=g, genomes=a[0], nested_per_genome=round(a[1] / a[0], 4), segments_per_genome=round(a[2] / a[0], 3),
             genomes_with_nested=a[3]) for (z, g), a in sorted(agg.items())]
w = csv.DictWriter(open(f"{OUT}/d5_seg_summary.tsv", "w"), fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)

# --- linkage level ---
lrows = []
for z in ZS:
    for k in sorted(subset) + ["extra/pos_control"]:
        d = f"{D5}/link/z{z}/{k.replace('/', '__')}"
        if not os.path.exists(f"{d}/cn_loci.tsv"):
            lrows.append(dict(z=z, key=k, ok=False)); continue
        c = rd(f"{d}/cn_candidates.tsv"); loci = rd(f"{d}/cn_loci.tsv")
        dec = Counter(x["decision"] for x in c)
        written = [x for x in loci if x.get("written", "yes") == "yes"] if loci and "written" in loci[0] else loci
        segs = {(s["chrom"], s["start"], s["end"]) for s in rd(f"{D5}/seg/z{z}/{k.replace('/', '__')}.tsv") if s["nested"] == "True"}
        def has_nested(x):
            for part in (x.get("segments") or "").split(";"):
                if ":" in part:
                    ch, rest = part.split(":", 1); a, b = rest.split(":")[0].split("-")
                    if (ch, a, b) in segs:
                        return True
            return False
        lrows.append(dict(z=z, key=k, ok=True, candidates=len(c), locus=dec.get("locus", 0), cn_only=dec.get("cn_only", 0),
                          public_cn=dec.get("public_cn", 0), simple_cnv=dec.get("simple_cnv", 0), cn_loci=len(loci),
                          cn_loci_nested=sum(has_nested(x) for x in loci),
                          cn_only_nested=sum(x["nested"] == "True" for x in c if x["decision"] == "cn_only"),
                          candidates_nested=sum(x["nested"] == "True" for x in c),
                          loci_classes=";".join(sorted(f"{x['locus_id']}={x['cn_class']}" for x in loci))))
fn = sorted({f for r in lrows for f in r}, key=lambda f: list(lrows[0]).index(f) if f in lrows[0] else 99)
w = csv.DictWriter(open(f"{OUT}/d5_link_per_genome.tsv", "w"), fieldnames=fn, delimiter="\t"); w.writeheader(); w.writerows(lrows)
lsum = []
for z in ZS:
    g = [r for r in lrows if r["z"] == z and r["key"] != "extra/pos_control" and r["ok"]]
    lsum.append(dict(z=z, genomes=len(g), **{f"{f}_per_genome": round(sum(r[f] for r in g) / len(g), 4)
                     for f in ("candidates", "candidates_nested", "locus", "cn_loci", "cn_loci_nested", "cn_only", "cn_only_nested")}))
w = csv.DictWriter(open(f"{OUT}/d5_link_summary.tsv", "w"), fieldnames=list(lsum[0]), delimiter="\t"); w.writeheader(); w.writerows(lsum)

# --- positive control events ---
ev = yaml.safe_load(open(EVENTS))["events"]
win = {}
for e in ev:
    iv = [v for v in (e.get("segments") or {}).values()]
    byc = defaultdict(list)
    for ch, a, b in iv:
        byc[ch] += [int(a), int(b)]
    if not iv:   # literature add-ons L01-L04: windows from addons/<id>/truth/truth_events.bed
        bed = os.path.join(os.path.dirname(EVENTS), e["id"], "truth", "truth_events.bed")
        for line in open(bed):
            f = line.split("\t")
            if len(f) >= 3 and not line.startswith(("#", "track")):
                byc[f[0]] += [int(f[1]), int(f[2])]
    win[e["id"]] = {ch: (min(v), max(v)) for ch, v in byc.items()}
prow = []
for z in ZS:
    segs = rd(f"{D5}/seg/z{z}/extra__pos_control.tsv")
    loci = rd(f"{D5}/link/z{z}/extra__pos_control/cn_loci.tsv") if os.path.exists(f"{D5}/link/z{z}/extra__pos_control/cn_loci.tsv") else []
    cands = rd(f"{D5}/link/z{z}/extra__pos_control/cn_candidates.tsv") if loci or os.path.exists(f"{D5}/link/z{z}/extra__pos_control/cn_candidates.tsv") else []
    for eid, w_ in win.items():
        def inwin(ch, a, b, pad=50000):
            return ch in w_ and int(b) > w_[ch][0] - pad and int(a) < w_[ch][1] + pad
        ss = [s for s in segs if inwin(s["chrom"], s["start"], s["end"])]
        ll = [x for x in loci if inwin(*([x["locus_id"].split(":")[0]] + x["locus_id"].split(":")[1].split("-")))]
        cc = [x for x in cands if any(inwin(p.split(":")[0], *p.split(":")[1].split("-")) for p in x["parts"].split(";") if ":" in p)]
        prow.append(dict(z=z, event=eid, segments=";".join(f"{s['chrom']}:{s['start']}-{s['end']} {s['state']}/{s['ploidy']}{' nested' if s['nested']=='True' else ''} z{s['z']}" for s in ss),
                         candidates=";".join(f"{x['decision']}:{x['cls']}" for x in cc),
                         cn_loci=";".join(f"{x['locus_id']}={x['cn_class']}" for x in ll)))
w = csv.DictWriter(open(f"{OUT}/d5_pc_events.tsv", "w"), fieldnames=list(prow[0]), delimiter="\t"); w.writeheader(); w.writerows(prow)
print(open(f"{OUT}/d5_seg_summary.tsv").read()); print(open(f"{OUT}/d5_link_summary.tsv").read())
for r in prow:
    if r["segments"] or r["cn_loci"]:
        print(r)
