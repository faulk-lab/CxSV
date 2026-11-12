#!/usr/bin/env python3
"""panel_cn_segments.tsv of one swept z, written exactly as rule panel_cn_segments (cn.smk)."""
import csv, json, os, sys
D5, Z = sys.argv[1:3]
REL = "/mnt/sargasso1/CxSV/build/panel-v1.2-r2/results/panel/release"
mem = list(csv.DictReader(open(f"{REL}/cn/pon/pon_members.tsv"), delimiter="\t"))
out = f"{D5}/panel_tables/panel_cn_segments.z{Z}.tsv"; os.makedirs(os.path.dirname(out), exist_ok=True)
head = None; missing = []
with open(out + ".tmp", "w") as fo:
    for m in mem:
        tag = m["key"].replace("/", "__")
        p = f"{D5}/seg/z{Z}/{tag}"
        if not os.path.exists(p + ".json"):
            missing.append(m["key"]); continue
        st = json.load(open(p + ".json"))
        with open(p + ".tsv") as fh:
            h = fh.readline().rstrip("\n")
            if head is None:
                head = h; fo.write("key\tindividual_id\tstratum\tsex\t" + h + "\n")
            for line in fh:
                fo.write(f"{m['key']}\t{m['individual_id']}\t{st['stratum']}\t{st.get('sex_call', st['sex'])}\t{line}")
if missing:
    sys.exit(f"z{Z}: {len(missing)} libraries missing, e.g. {missing[:3]}")
os.replace(out + ".tmp", out); print(out)
