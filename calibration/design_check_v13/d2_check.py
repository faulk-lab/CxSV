#!/usr/bin/env python3
"""panel-v1.3 design check D2: chrX/chrY copy-number public counting, v1.2 rule vs D2 rule.

v1.2 rule = origin/main code (cxsv_cnlink.cn_public + cn_linkage.others): a chrX/chrY segment
is public when >= 1 other sex-matched panel library carries it (same direction, reciprocal
overlap >= 0.5). D2 rule = branch code (PR #43): >= 2 distinct sex-matched individuals, and the
libraries of a sex-discordant individual (cxsv_cnlink.sexchrom_discordant) count for no one.
Both rules exclude the sample's own individual (cn_linkage.others).

  d2_check.py --old OLD_REPO --new NEW_REPO --design design.tsv --out OUT
design.tsv: set, key, individual_id, segments, stats, bound   (one row per genome and bound)
bounds.tsv: bound, panel_segments, members   (members: key individual_id sex_call x_ratio y_ratio)
"""
import argparse, csv, importlib.util, json, os, sys
from collections import defaultdict

L01 = ("chrX", 101_700_000, 102_550_000)


def load(repo, tag):
    """Load cxsv_cnlink + cn_linkage of `repo` as separate module objects."""
    for m in ("cxsv_cnlink",):
        sys.modules.pop(m, None)
    sys.path.insert(0, f"{repo}/workflow/lib")
    spec = importlib.util.spec_from_file_location(f"cnlink_{tag}", f"{repo}/workflow/lib/cxsv_cnlink.py")
    cl = importlib.util.module_from_spec(spec); spec.loader.exec_module(cl)
    sys.modules["cxsv_cnlink"] = cl
    spec = importlib.util.spec_from_file_location(f"cnlinkage_{tag}", f"{repo}/workflow/scripts/cn_linkage.py")
    lk = importlib.util.module_from_spec(spec); spec.loader.exec_module(lk)
    sys.path.pop(0)
    sys.modules.pop("cxsv_cnlink", None)
    return cl, lk


def rd(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader([l for l in fh if not l.startswith("#")], delimiter="\t"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--old", required=True); ap.add_argument("--new", required=True)
    ap.add_argument("--design", required=True); ap.add_argument("--bounds", required=True)
    ap.add_argument("--calling-json", required=True, help="branch cn_params (calling.yaml as JSON)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    cl_old, lk_old = load(a.old, "old")
    cl_new, lk_new = load(a.new, "new")
    cnp = json.load(open(a.calling_json))["cn_params"]
    P_old = {**cl_old.CN_DEFAULTS, **{k: v for k, v in cnp["linkage"].items() if k in cl_old.CN_DEFAULTS}}
    P_old["public_min_genomes"] = 1
    P_new = {**cl_new.CN_DEFAULTS, **cnp["linkage"]}
    bounds = {r["bound"]: r for r in rd(a.bounds)}
    design = rd(a.design)
    seg_rows, summ, l01_rows, disc_rows = [], [], [], []
    cache = {}
    for b, br in bounds.items():
        members = rd(br["members"])
        excl = cl_new.sexchrom_discordant(members, P_new)
        for k, why in sorted(excl.items()):
            disc_rows.append(dict(bound=b, key=k, reason=why))
        cache[b] = (br["panel_segments"], excl)
    # others per (bound, individual) built lazily; panel table read once per bound
    panel_rows = {}
    for b, (ps, excl) in cache.items():
        panel_rows[b] = rd(ps)
    def others_for(b, individual, rule):
        excl = cache[b][1]
        by, sex = {}, {}
        for r in panel_rows[b]:
            if r["individual_id"] == individual:
                continue
            k = (r["key"], r["individual_id"])
            sex[k] = (r.get("sex") or "").strip()
            by.setdefault(k, []).append(dict(chrom=r["chrom"], start=int(r["start"]), end=int(r["end"]),
                                             state=int(r["state"]), ploidy=int(r["ploidy"])))
        if rule == "old":
            return [(k, ind, s, sex[(k, ind)]) for (k, ind), s in by.items()]
        return [(k, ind, s, sex[(k, ind)], k not in excl) for (k, ind), s in by.items()]
    for g in design:
        b = g["bound"]
        st = json.load(open(g["stats"]))
        psex = st.get("sex"); scall = st.get("sex_call", psex)
        segs = [s for s in lk_new.read_segments(g["segments"]) if s["chrom"] in ("chrX", "chrY") and s["state"] != s["ploidy"]]
        o_old = others_for(b, g["individual_id"], "old"); o_new = others_for(b, g["individual_id"], "new")
        n_pub_old = {"chrX": 0, "chrY": 0}; n_pub_new = {"chrX": 0, "chrY": 0}; n = {"chrX": 0, "chrY": 0}
        for s in segs:
            c = {"segs": [s]}
            _, co = cl_old.cn_public(c, o_old, P_old, sex=psex)
            _, cn_ = cl_new.cn_public(c, o_new, P_new, sex=psex)
            po, pn = co[0] >= 1, cn_[0] >= cl_new.min_genomes(s["chrom"], P_new)
            n[s["chrom"]] += 1; n_pub_old[s["chrom"]] += po; n_pub_new[s["chrom"]] += pn
            seg_rows.append(dict(set=g["set"], bound=b, key=g["key"], sex_call=scall, ploidy_sex=psex,
                                 chrom=s["chrom"], start=s["start"], end=s["end"], state=s["state"], ploidy=s["ploidy"],
                                 direction="gain" if s["state"] > s["ploidy"] else "loss", nested=s["nested"],
                                 count_v12=co[0], public_v12=po, count_d2=cn_[0], public_d2=pn))
        summ.append(dict(set=g["set"], bound=b, key=g["key"], sex_call=scall, ploidy_sex=psex,
                         x_ratio=st.get("x_ratio"), y_ratio=st.get("y_ratio"),
                         chrX_segments=n["chrX"], chrX_public_v12=n_pub_old["chrX"], chrX_public_d2=n_pub_new["chrX"],
                         chrY_segments=n["chrY"], chrY_public_v12=n_pub_old["chrY"], chrY_public_d2=n_pub_new["chrY"]))
        if g["set"] == "positive_control":
            l01 = [s for s in segs if s["chrom"] == L01[0] and s["end"] > L01[1] and s["start"] < L01[2]]
            c = {"segs": l01}
            pu_o, co = cl_old.cn_public(c, o_old, P_old, sex=psex) if l01 else (False, [])
            pu_n, cn_ = cl_new.cn_public(c, o_new, P_new, sex=psex) if l01 else (False, [])
            # which individuals match, per segment (D2 rule, before / after the discordance exclusion)
            def who(s, oo, excl_on):
                out = []
                for o in oo:
                    if psex and len(o) > 3 and o[3] != psex:
                        continue
                    if excl_on and len(o) > 4 and not o[4]:
                        continue
                    sign = s["state"] > s["ploidy"]
                    if any(t["chrom"] == s["chrom"] and (t["state"] > t["ploidy"]) == sign and
                           cl_new._recip(s["start"], s["end"], t["start"], t["end"]) >= P_new["public_overlap"] for t in o[2]):
                        out.append(o[0] + ("" if len(o) < 5 or o[4] else "[sex-discordant]"))
                return out
            l01_rows.append(dict(bound=b, pc_sex_call=scall, pc_x_ratio=st.get("x_ratio"),
                                 l01_segments=";".join(f"{s['chrom']}:{s['start']}-{s['end']} {s['state']}/{s['ploidy']}{' nested' if s['nested'] else ''}" for s in l01),
                                 counts_v12=",".join(map(str, co)), public_v12=pu_o, counts_d2=",".join(map(str, cn_)), public_d2=pu_n,
                                 private_d2=(not pu_n) and bool(l01),
                                 matching_sex_matched=";".join(sorted({x for s in l01 for x in who(s, o_new, False)})),
                                 matching_after_d2_exclusion=";".join(sorted({x for s in l01 for x in who(s, o_new, True)}))))
    def wr(name, rows):
        if not rows:
            return
        with open(f"{a.out}/{name}", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t"); w.writeheader(); w.writerows(rows)
    wr("d2_segments.tsv", seg_rows); wr("d2_per_genome.tsv", summ); wr("d2_l01.tsv", l01_rows); wr("d2_discordant.tsv", disc_rows)
    # summary by set, bound, sex call
    agg = defaultdict(lambda: defaultdict(float))
    for r in summ:
        k = (r["set"], r["bound"], r["sex_call"])
        agg[k]["genomes"] += 1
        for f in ("chrX_segments", "chrX_public_v12", "chrX_public_d2", "chrY_segments", "chrY_public_v12", "chrY_public_d2"):
            agg[k][f] += r[f]
    rows = []
    for (s, b, x), v in sorted(agg.items()):
        g = v["genomes"]
        rows.append(dict(set=s, bound=b, sex_call=x, genomes=int(g),
                         **{f: int(v[f]) for f in ("chrX_segments", "chrX_public_v12", "chrX_public_d2", "chrY_segments", "chrY_public_v12", "chrY_public_d2")},
                         chrX_public_v12_per_genome=round(v["chrX_public_v12"] / g, 3), chrX_public_d2_per_genome=round(v["chrX_public_d2"] / g, 3)))
    wr("d2_summary_by_sex.tsv", rows)
    print(json.dumps(l01_rows, indent=1)); print(len(disc_rows), "discordant library rows")


if __name__ == "__main__":
    main()
