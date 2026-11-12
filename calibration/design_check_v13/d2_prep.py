#!/usr/bin/env python3
"""Build design.tsv / bounds.tsv for d2_check.py (panel-v1.3 design check)."""
import csv, json, os, sys
OUT = sys.argv[1]
B12 = "/mnt/sargasso1/CxSV/build/panel-v1.2-r2/results/panel"
BC = "/mnt/sargasso1/CxSV/calibration/panel_v12_r2/bound_check"
DC2 = "/mnt/sargasso1/CxSV/calibration/panel_v12_r2/design_check_v2"
RUN = "/mnt/sargasso1/CxSV/patient_test_runs/monkeymeth_20261110"
PC12 = "/home/minknow/scratch/item15/v12r2bundle/pc_full_L_tag/run/results/patient/samples/pos_control/pos_control/cn"
TEST = ["BOC_1882", "BOC_2057", "BOC_2865", "BOC_3537", "BOC_3559", "BOC_4433", "BOC_5002", "BOC_5883", "BOC_7308"]
rd = lambda p: list(csv.DictReader(open(p), delimiter="\t"))
pg = rd(f"{DC2}/per_genome.tsv")
sub = [r for r in pg if r["in_subset"] == "True"]
s100 = {f"{r['cohort']}/{r['sample_id']}" for r in rd("/home/minknow/scratch/v13/design_check/repo/calibration/cn_linkage/subset100.tsv")}
extra = [r["key"] for r in sub if r["key"] not in s100]
json.dump(dict(subset100=sorted(s100), l01_matchers_not_in_subset100=extra, n=len(sub)), open(f"{OUT}/design_set.json", "w"), indent=1)
ind = {r["key"]: r["individual_id"] for r in rd(f"{B12}/release/cn/pon/pon_members.tsv")}
rows = []
for r in sub:
    k = r["key"]; c, s = k.split("/")
    rows.append(dict(set="subset", key=k, individual_id=ind[k], segments=f"{B12}/samples/{c}/{s}/cn/segments.tsv",
                     stats=f"{B12}/samples/{c}/{s}/cn/stats.json", bound="1.12"))
    for b in ("1.18", "1.25"):
        rows.append(dict(set="subset", key=k, individual_id=ind[k], segments=f"{BC}/design_f{b}/segments/{c}__{s}.tsv",
                         stats=f"{BC}/design_f{b}/segments/{c}__{s}.json", bound=b))
for t in TEST:
    S = f"{RUN}/{t}/results/patient/samples/monkeymeth_20261110/{t}/cn"
    rows.append(dict(set="test_genome", key=t, individual_id=t, segments=f"{S}/segments.tsv", stats=f"{S}/stats.json", bound="1.12"))
rows.append(dict(set="positive_control", key="pos_control", individual_id="pos_control", segments=f"{PC12}/segments.tsv", stats=f"{PC12}/stats.json", bound="1.12"))
for b in ("1.18", "1.25"):
    rows.append(dict(set="positive_control", key="pos_control", individual_id="pos_control",
                     segments=f"{BC}/design_f{b}/segments/extra__pos_control.tsv", stats=f"{BC}/design_f{b}/segments/extra__pos_control.json", bound=b))
for r in rows:
    assert os.path.exists(r["segments"]) and os.path.exists(r["stats"]), r
w = csv.DictWriter(open(f"{OUT}/design.tsv", "w"), fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
# members per bound
brows = [dict(bound="1.12", panel_segments=f"{B12}/release/cn/panel_cn_segments.tsv", members=f"{B12}/release/cn/pon/pon_members.tsv")]
for b in ("1.18", "1.25"):
    D = f"{BC}/design_f{b}"
    mem = []
    for m in rd(f"{D}/members.tsv"):
        c, s = m["key"].split("/")
        st = json.load(open(f"{D}/segments/{c}__{s}.json"))
        mem.append(dict(key=m["key"], individual_id=m["individual_id"], sex_call=st["sex_call"], x_ratio=st["x_ratio"], y_ratio=st["y_ratio"]))
    p = f"{OUT}/members_f{b}.tsv"
    w = csv.DictWriter(open(p, "w"), fieldnames=list(mem[0]), delimiter="\t"); w.writeheader(); w.writerows(mem)
    brows.append(dict(bound=b, panel_segments=f"{D}/panel_cn_segments.tsv", members=p))
w = csv.DictWriter(open(f"{OUT}/bounds.tsv", "w"), fieldnames=list(brows[0]), delimiter="\t"); w.writeheader(); w.writerows(brows)
print(len(sub), "subset genomes;", len(extra), "L01 matchers outside subset100;", len(rows), "design rows")
