#!/usr/bin/env python3
"""Design check of the sex-stratified chrX/chrY copy-number reference (roadmap item 15,
panel-v1.2-r2; PI decision 2026-11-09). Calibration only, not part of a run.

With the code of --code (its workflow/scripts/cn_pon.py and cn_segment.py, its
config/calling.yaml), on the panel-v1.2 depth tracks (cn_depth_v1):

  1. PoN per stratum from the panel members (as rule panel_cn_pon: same members, same
     library-QC exclusions) -> OUT/pon/
  2. every panel library segmented against its stratum's new PoN (as rule cn_segment)
     -> OUT/segments/<cohort>__<sample>.{tsv,json}; the panel table as rule
     panel_cn_segments (sex column = sex call) -> OUT/panel_cn_segments.tsv
  3. the positive control's depth (--pc-depth) segmented the same way; its segments in the
     L01 window and their panel counts under the new (sex-matched) and the old (no sex
     matching) counting rule, with the matching genomes -> OUT/l01_public.tsv
  4. per genome, old (--old-samples: <root>/<cohort>/<sample>/cn/segments.tsv of the
     panel-v1.2 candidate) vs new: chrX / chrY gains and losses, the share of chrX gains that
     are state 2 at ploidy 1, and whether the autosomal segments are identical
     -> OUT/per_genome.tsv; summaries by design set, platform stratum and sex
     -> OUT/summary_by_sex.tsv, OUT/summary.json

Design sets: 'subset' = calibration/cn_linkage/subset100.tsv (the Phase 1 subset) plus
--extra-keys (the 14 genomes that matched L01 in panel-v1.2); 'panel' = every library.

  sex_design_check.py --code REPO --members pon_input.tsv --depth-root cn_depth_v1 \
      --old-samples panel-v1.2/results/panel/samples --pc-depth PC/cn/depth.regions.bed.gz \
      --extra-keys K1,K2,... --out OUT [--threads 16]
"""
import argparse
import csv
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import yaml

L01 = ("chrX", 101_700_000, 102_550_000)


def read_tsv(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def run(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    return r.returncode, r.stderr[-2000:]


def seg_rows(p):
    return [dict(r, start=int(r["start"]), end=int(r["end"]), state=int(r["state"]), ploidy=int(r["ploidy"]))
            for r in read_tsv(p)] if os.path.exists(p) else None


def counts(segs, chroms):
    g = sum(1 for s in segs if s["chrom"] in chroms and s["state"] > s["ploidy"])
    l = sum(1 for s in segs if s["chrom"] in chroms and s["state"] < s["ploidy"])
    g21 = sum(1 for s in segs if s["chrom"] in chroms and s["state"] == 2 and s["ploidy"] == 1)
    return g, l, g21


def auto_key(segs):
    return sorted((s["chrom"], s["start"], s["end"], s["state"], s["ploidy"]) for s in segs
                  if s["chrom"] not in ("chrX", "chrY"))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--code", required=True)
    ap.add_argument("--members", required=True, help="panel-v1.2 panel/cn/pon_input.tsv")
    ap.add_argument("--depth-root", required=True)
    ap.add_argument("--old-samples", required=True)
    ap.add_argument("--pc-depth", required=True)
    ap.add_argument("--pc-stratum", default="ONT_LRSC_R10")
    ap.add_argument("--extra-keys", default="")
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=16)
    a = ap.parse_args()
    os.makedirs(f"{a.out}/segments", exist_ok=True)
    scripts = f"{a.code}/workflow/scripts"
    sys.path.insert(0, f"{a.code}/workflow/lib")
    import cxsv_cnlink as cl  # noqa: E402

    # 0. cn_params.json as rule cn_params_json
    cfg = yaml.safe_load(open(f"{a.code}/config/calling.yaml"))
    params = f"{a.out}/cn_params.json"
    with open(params, "w") as fh:
        json.dump({"cn_params": cfg["cn_params"], "cxsv_params": cfg["cxsv_params"],
                   "junction_params": cfg["junction_params"]}, fh, sort_keys=True)
    members = read_tsv(a.members)
    for m in members:
        m["depth"] = f"{a.depth_root}/{m['key']}/depth.regions.bed.gz"
    mem = f"{a.out}/members.tsv"
    with open(mem, "w") as fh:
        fh.write("key\tindividual_id\tstratum\tdepth\n")
        for m in members:
            fh.write(f"{m['key']}\t{m['individual_id']}\t{m['stratum']}\t{m['depth']}\n")

    # 1. PoN
    rc, err = run([sys.executable, f"{scripts}/cn_pon.py", "--members", mem, "--params", params,
                   "--out-dir", f"{a.out}/pon", "--threads", str(a.threads)])
    open(f"{a.out}/cn_pon.log", "w").write(err)
    if rc:
        sys.exit(f"cn_pon.py failed: {err}")

    # 2. + 3. segmentation
    def seg(job):
        key, depth, stratum = job
        base = f"{a.out}/segments/{key.replace('/', '__')}"
        return key, run([sys.executable, f"{scripts}/cn_segment.py", "--depth", depth, "--pon",
                         f"{a.out}/pon/{stratum}.npz", "--params", params, "--stratum", stratum,
                         "--default-stratum", stratum, "--out", base + ".tsv", "--stats", base + ".json"])
    jobs = [(m["key"], m["depth"], m["stratum"]) for m in members] + [("extra/pos_control", a.pc_depth, a.pc_stratum)]
    failed = []
    with ThreadPoolExecutor(a.threads) as ex:
        for key, (rc, err) in ex.map(seg, jobs):
            if rc:
                failed.append((key, err))
    if failed:
        open(f"{a.out}/segment_failures.txt", "w").write("\n".join(f"{k}\t{e}" for k, e in failed))
    by_key = {m["key"]: m for m in members}

    def new(key):
        base = f"{a.out}/segments/{key.replace('/', '__')}"
        return seg_rows(base + ".tsv"), json.load(open(base + ".json"))

    # panel table (rule panel_cn_segments)
    others_new, others_old = [], []
    with open(f"{a.out}/panel_cn_segments.tsv", "w") as fo:
        fo.write("key\tindividual_id\tstratum\tsex\tchrom\tstart\tend\tstate\tploidy\n")
        for m in members:
            segs, st = new(m["key"])
            for s in segs:
                fo.write(f"{m['key']}\t{m['individual_id']}\t{st['stratum']}\t{st['sex_call']}\t"
                         f"{s['chrom']}\t{s['start']}\t{s['end']}\t{s['state']}\t{s['ploidy']}\n")
            others_new.append((m["key"], m["individual_id"], segs, st["sex_call"]))
            old = seg_rows(f"{a.old_samples}/{m['key']}/cn/segments.tsv") or []
            others_old.append((m["key"], m["individual_id"], old))

    # 3. L01
    pc, pst = new("extra/pos_control")
    win = [s for s in pc if s["chrom"] == L01[0] and s["end"] > L01[1] and s["start"] < L01[2]]
    P = {**cl.CN_DEFAULTS, **(cfg["cn_params"].get("linkage") or {})}
    rows = []
    for s in win:
        c = {"segs": [s]}
        for rule, oth, sex in (("new_sex_matched", others_new, pst["sex"]), ("new_no_sex_matching", [o[:3] for o in others_new], None),
                               ("old_panel_v1.2", others_old, None)):
            n = cl.cn_public(c, oth, P, sex=sex)[1][0]
            who = []
            sign = s["state"] > s["ploidy"]
            for o in oth:
                if sex and len(o) > 3 and o[3] != sex:
                    continue
                for t in o[2]:
                    if t["chrom"] == s["chrom"] and (t["state"] > t["ploidy"]) == sign and \
                            cl._recip(s["start"], s["end"], t["start"], t["end"]) >= P["public_overlap"]:
                        who.append(f"{o[0]}({t['start']}-{t['end']} {t['state']}/{t['ploidy']}"
                                   f"{' ' + o[3] if len(o) > 3 else ''})")
                        break
            rows.append(dict(pc_segment=f"{s['chrom']}:{s['start']}-{s['end']}", state=s["state"], ploidy=s["ploidy"],
                             nested=s.get("nested"), rule=rule, n_genomes=n, genomes=";".join(who)))
    cand_all = {"segs": win}
    pub_new = cl.cn_public(cand_all, others_new, P, sex=pst["sex"]) if win else (None, [])
    with open(f"{a.out}/l01_public.tsv", "w") as fh:
        cols = ["pc_segment", "state", "ploidy", "nested", "rule", "n_genomes", "genomes"]
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    # 4. per genome old vs new
    subset = {f"{r['cohort']}/{r['sample_id']}" for r in
              read_tsv(f"{a.code}/calibration/cn_linkage/subset100.tsv")}
    extra = {k for k in a.extra_keys.split(",") if k}
    per = []
    for (key, ind, newsegs, sexcall), (_, _, oldsegs) in zip(others_new, others_old):
        st = new(key)[1]
        old_st = f"{a.old_samples}/{key}/cn/stats.json"
        old_sex = json.load(open(old_st)).get("sex") if os.path.exists(old_st) else ""
        gx_o, lx_o, g21_o = counts(oldsegs, ("chrX",))
        gx_n, lx_n, g21_n = counts(newsegs, ("chrX",))
        gy_o, ly_o, _ = counts(oldsegs, ("chrY",))
        gy_n, ly_n, _ = counts(newsegs, ("chrY",))
        ao, an = auto_key(oldsegs), auto_key(newsegs)
        per.append(dict(key=key, individual_id=ind, stratum=by_key[key]["stratum"], sex_call=sexcall,
                        ploidy_sex=st["sex"], old_sex=old_sex, x_ratio=st.get("x_ratio"), y_ratio=st.get("y_ratio"),
                        in_subset=key in subset or key.split("/")[1] in extra,
                        chrX_gain_old=gx_o, chrX_gain_new=gx_n, chrX_loss_old=lx_o, chrX_loss_new=lx_n,
                        chrX_gain_s2p1_old=g21_o, chrX_gain_s2p1_new=g21_n,
                        chrY_old=gy_o + ly_o, chrY_new=gy_n + ly_n,
                        auto_old=len(ao), auto_new=len(an), auto_identical=ao == an))
    cols = list(per[0])
    with open(f"{a.out}/per_genome.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in per:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    summ = []
    for setname, sel in (("subset", [r for r in per if r["in_subset"]]), ("panel", per)):
        groups = {}
        for r in sel:
            for st in (r["stratum"], "all"):
                groups.setdefault((st, r["sex_call"]), []).append(r)
        for (st, sx), rs in sorted(groups.items()):
            n = len(rs)
            sm = lambda k: sum(r[k] for r in rs)  # noqa: E731
            summ.append(dict(set=setname, stratum=st, sex_call=sx, genomes=n,
                             chrX_gain_old=round(sm("chrX_gain_old") / n, 3), chrX_gain_new=round(sm("chrX_gain_new") / n, 3),
                             chrX_loss_old=round(sm("chrX_loss_old") / n, 3), chrX_loss_new=round(sm("chrX_loss_new") / n, 3),
                             s2p1_share_old=round(sm("chrX_gain_s2p1_old") / max(1, sm("chrX_gain_old")), 3),
                             s2p1_share_new=round(sm("chrX_gain_s2p1_new") / max(1, sm("chrX_gain_new")), 3),
                             chrY_old=round(sm("chrY_old") / n, 3), chrY_new=round(sm("chrY_new") / n, 3),
                             auto_old=round(sm("auto_old") / n, 2), auto_new=round(sm("auto_new") / n, 2),
                             auto_identical=sum(r["auto_identical"] for r in rs)))
    cols = list(summ[0])
    with open(f"{a.out}/summary_by_sex.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in summ:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    json.dump(dict(failures=len(failed), pc_stats=pst, pc_l01_segments=[f"{s['chrom']}:{s['start']}-{s['end']} {s['state']}/{s['ploidy']}" for s in win],
                   l01_public_new=pub_new[0], l01_counts_new=pub_new[1],
                   subset_genomes=sum(r["in_subset"] for r in per), panel_genomes=len(per)),
              open(f"{a.out}/summary.json", "w"), indent=1)
    print(open(f"{a.out}/summary.json").read())


if __name__ == "__main__":
    main()
