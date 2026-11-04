#!/usr/bin/env python3
"""Roadmap item 13, steps 1-2: characterise the oversized clusters of panel genomes and
simulate candidate clustering rules on every genome.

For each panel sample this script

  1. re-implements `bedtools cluster -d D` on the sample's own cxsv/merged_sv.bed and
     checks that it reproduces the released cxsv/clusters.bed partition exactly
     (`v0_partition_identical` in summary.json);
  2. describes every oversized cluster (span > max_cluster_size) of the released
     partition: records, types, the largest record interval, and the records whose
     interval alone exceeds max_cluster_size ("giant" records);
  3. writes one row per giant record with its merged-VCF evidence (SUPP_VEC, per-caller
     genotype, breakend-notation ALT or not, caller ids) and the original caller
     evidence (Sniffles2 SUPPORT / COVERAGE / AF, cuteSV RE, SVIM SUPPORT);
  4. re-clusters the sample under candidate rules (VARIANTS), runs the unchanged
     classifier (workflow/lib/cxsv_classes.analyse_locus) on every cluster, and reports
     per rule: clusters, oversized clusters, blind span, CxSV loci, and CxSV loci that
     are identical to / new relative to / lost from the current rule.

Read-only on the panel build. Output: one directory per sample under --out.

Usage
  characterise.py --manifest sample_manifest.tsv --samples-root DIR --code REPO \
      --fai hs1.fa.fai --out OUTDIR [--shard i/k] [--threads N] [--samples S1,S2]
"""
import argparse
import csv
import json
import multiprocessing as mp
import os
import sys
import time
from collections import defaultdict

G = {}

# Discovery parameters as fingerprinted in config/calling.yaml (panel-v1.0).
CLUSTER_D = 50000
MAX_CLUSTER = 500000
BND_WINDOW = 1
PARAMS = {"collapse_distance": 300, "ins_len_ratio": 0.5, "composite_fraction": 0.25,
          "min_breakpoints": 3, "chromoanag_min_window": 10000, "chromoanag_max_window": 50000,
          "chromoanag_min_svs": 3, "extreme_min_breakends": 6, "min_interleaved": 3,
          "interleave_alt_min_intra": 3, "interleave_alt_min_inter": 4}
SPAN = ("DEL", "DUP", "INV")

# Candidate rules. ends_over: span records whose interval is longer than this link only at
# their two breakends (None = never); bnd_ends: records written in breakend (bracket) ALT
# notation link only at their breakends whatever their length; drop_single_over: records
# longer than this supported by ONE caller are left out of clustering and of the locus;
# drop_all_over: records longer than this are left out whatever their support (the
# rule adopted: cxsv_params.long_span_rule exclude).
VARIANTS = {
    "V0_current":     dict(ends_over=None,   bnd_ends=False, drop_single_over=None),
    "S500_ends":      dict(ends_over=500000, bnd_ends=False, drop_single_over=None),
    "S500_ends_x1":   dict(ends_over=500000, bnd_ends=False, drop_single_over=500000),
    "B_S500_ends":    dict(ends_over=500000, bnd_ends=True,  drop_single_over=None),
    "B_S500_ends_x1": dict(ends_over=500000, bnd_ends=True,  drop_single_over=500000),
    "S50_ends":       dict(ends_over=50000,  bnd_ends=False, drop_single_over=None),
    "S500_exclude":   dict(ends_over=None,   bnd_ends=False, drop_single_over=None, drop_all_over=500000),
}


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_bed(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 6:
                rows.append((f[0], int(f[1]), int(f[2]), f[3], f[4], f[5]))
    return rows


def bt_cluster(rows, d):
    """bedtools cluster -d d on rows (any order): [(row, cluster_index)] in `bedtools
    sort` order (chrom lexicographic, start, end). A row joins the open cluster when
    its start is at most d bp past the cluster's largest end so far."""
    srt = sorted(rows, key=lambda r: (r[0], r[1], r[2]))
    out, cid, chrom, max_end = [], 0, None, None
    for r in srt:
        if r[0] != chrom or r[1] - max_end > d:
            cid += 1
            chrom, max_end = r[0], r[2]
        else:
            max_end = max(max_end, r[2])
        out.append((r, cid))
    return out


def merged_span(iv):
    tot, cur = 0, None
    for c, s, e in sorted(iv):
        if cur and cur[0] == c and s <= cur[2]:
            cur[2] = max(cur[2], e)
        else:
            if cur:
                tot += cur[2] - cur[1]
            cur = [c, s, e]
    if cur:
        tot += cur[2] - cur[1]
    return tot


def vcf_info(path, detail):
    """{id: (supp_vec, breakend-notation ALT)} for every merged-VCF record; records whose
    id is a key of `detail` also get their per-caller evidence written into detail[id]."""
    out = {}
    with open(path) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 8)
            vid, alt, info = f[2], f[4], f[7]
            i = info.find("SUPP_VEC=")
            sv = info[i + 9:i + 12] if i >= 0 else "..."
            bnd_alt = ("[" in alt or "]" in alt)
            out[vid] = (sv, bnd_alt)
            if vid in detail:
                full = line.rstrip("\n").split("\t")
                fmt = full[8].split(":")
                d = {"supp_vec": sv, "bnd_alt": bnd_alt, "qual": full[5]}
                for name, s in zip(("cuteSV", "SVIM", "Sniffles2"), full[9:12]):
                    c = dict(zip(fmt, s.split(":")))
                    d[f"{name}_gt"] = c.get("GT", ".")
                    d[f"{name}_id"] = c.get("ID", "NaN")
                    d[f"{name}_dr"] = c.get("DR", ".")
                detail[vid] = d
    return out


def caller_evidence(sdir, want):
    """want: {caller: {caller_record_id: merged_record_id}} -> {merged_id: {fields}}."""
    files = {"cuteSV": "cutesv_filtered.vcf", "SVIM": "svim_filtered.vcf", "Sniffles2": "sniffles_filtered.vcf"}
    out = defaultdict(dict)
    for caller, ids in want.items():
        p = f"{sdir}/callers/{files[caller]}"
        if not ids or not os.path.exists(p):
            continue
        with open(p) as fh:
            for line in fh:
                if line[0] == "#":
                    continue
                f = line.split("\t", 3)
                if f[2] not in ids:
                    continue
                full = line.rstrip("\n").split("\t")
                info = dict(kv.split("=", 1) if "=" in kv else (kv, "1") for kv in full[7].split(";"))
                fmt = full[8].split(":") if len(full) > 9 else []
                smp = dict(zip(fmt, full[9].split(":"))) if len(full) > 9 else {}
                o = out[ids[f[2]]]
                o[f"{caller}_svtype_orig"] = info.get("SVTYPE", ".")
                o[f"{caller}_support"] = info.get("SUPPORT", info.get("RE", "."))
                o[f"{caller}_precise"] = "PRECISE" if "PRECISE" in info else ("IMPRECISE" if "IMPRECISE" in info else ".")
                if caller == "Sniffles2":
                    o["Sniffles2_coverage"] = info.get("COVERAGE", ".")
                    o["Sniffles2_af"] = info.get("AF", ".")
                    o["Sniffles2_dv_dr"] = f"{smp.get('DV', '.')},{smp.get('DR', '.')}"
                elif caller == "cuteSV":
                    o["cuteSV_dv_dr"] = f"{smp.get('DV', '.')},{smp.get('DR', '.')}"
    return out


def transform(rows, ev, v):
    """Rows of merged_sv.bed under rule v; ev: {id: (supp_vec, bnd_alt)}.
    Returns (rows, ids left out)."""
    out, dropped = [], set()
    for r in rows:
        chrom, s, e, t, vid, inter = r
        L = e - s
        if t in SPAN and inter == "intra":
            sv, bnd_alt = ev.get(vid, ("...", False))
            if v["drop_single_over"] is not None and L > v["drop_single_over"] and sv.count("1") == 1:
                dropped.add(vid)
                continue
            if v.get("drop_all_over") is not None and L > v["drop_all_over"]:
                dropped.add(vid)
                continue
            if (v["ends_over"] is not None and L > v["ends_over"]) or (v["bnd_ends"] and bnd_alt):
                out.append((chrom, max(0, s - BND_WINDOW), s + BND_WINDOW, t, vid, inter))
                out.append((chrom, max(0, e - BND_WINDOW), e + BND_WINDOW, t, vid + "_end", inter))
                continue
        out.append(r)
    return out, dropped


def canon(vid):
    for suf in ("_mate", "_dest", "_end"):
        if vid.endswith(suf):
            return vid[: -len(suf)]
    return vid


def classify(clustered, rows_len, by_id, order, dropped):
    """[(cluster id, cluster dict, analyse_locus result or None when oversized)]."""
    cc = G["cc"]
    cl = {}
    for r, cid in clustered:
        c = cl.setdefault(cid, {"chrom": r[0], "start": r[1], "end": r[2], "ids": set(), "maxlen": 0})
        c["start"] = min(c["start"], r[1])
        c["end"] = max(c["end"], r[2])
        k = canon(r[4])
        c["ids"].add(k)
        c["maxlen"] = max(c["maxlen"], rows_len.get(k, 0))
    res = []
    for cid, c in cl.items():
        if c["end"] - c["start"] > MAX_CLUSTER:
            res.append((cid, c, None))
            continue
        js = [by_id[i] for i in sorted(c["ids"]) if i in by_id and i not in dropped]
        res.append((cid, c, cc.analyse_locus(js, c["chrom"], order, PARAMS)))
    return res


def one(task):
    cohort, sample, meta = task
    a, cj = G["a"], G["cj"]
    sdir = f"{a.samples_root}/{cohort}/{sample}"
    odir = f"{a.out}/{cohort}/{sample}"
    if os.path.exists(f"{odir}/done.ok") and not a.force:
        return (cohort, sample, "skip")
    os.makedirs(odir, exist_ok=True)
    t0 = time.time()
    rows = read_bed(f"{sdir}/cxsv/merged_sv.bed")
    rows_len = {}
    for r in rows:
        if r[3] in SPAN and r[5] == "intra":
            rows_len[r[4]] = r[2] - r[1]
    rel = defaultdict(set)
    with open(f"{sdir}/cxsv/clusters.bed") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rel[f[6]].add((f[0], int(f[1]), int(f[2]), f[4]))
    v0 = bt_cluster(rows, CLUSTER_D)
    mine = defaultdict(set)
    for r, cid in v0:
        mine[cid].add((r[0], r[1], r[2], r[4]))
    identical = set(map(frozenset, rel.values())) == set(map(frozenset, mine.values()))

    detail = {vid: None for vid, L in rows_len.items() if L > MAX_CLUSTER}
    ev = vcf_info(f"{sdir}/merged_sv.vcf", detail)
    want = defaultdict(dict)
    for vid, d in detail.items():
        if d:
            for c in ("cuteSV", "SVIM", "Sniffles2"):
                if d.get(f"{c}_id", "NaN") not in ("NaN", "NA", "."):
                    want[c][d[f"{c}_id"]] = vid
    cev = caller_evidence(sdir, want)

    order = cj.contig_order(a.fai)
    by_id = {}
    for j in cj.parse_vcf(f"{sdir}/merged_sv.vcf", order):
        by_id.setdefault(j.record_ids[0], j)

    cl0 = defaultdict(list)
    for r, cid in v0:
        cl0[cid].append(r)
    over_rows, over_iv = [], {}
    for cid, rs in cl0.items():
        s, e = min(r[1] for r in rs), max(r[2] for r in rs)
        if e - s <= MAX_CLUSTER:
            continue
        lens = [r[2] - r[1] for r in rs if r[3] in SPAN and r[5] == "intra"]
        types = defaultdict(int)
        for r in rs:
            types[r[3]] += 1
        gids = [r[4] for r in rs if rows_len.get(r[4], 0) > MAX_CLUSTER]
        over_iv[cid] = (rs[0][0], s, e)
        over_rows.append({"cohort": cohort, "sample": sample, "platform": meta["platform"],
                          "coverage_x": meta["coverage_x"], "cluster": cid, "chrom": rs[0][0], "start": s,
                          "end": e, "span": e - s, "n_rows": len(rs), "n_records": len({canon(r[4]) for r in rs}),
                          "types": ",".join(f"{k}:{types[k]}" for k in sorted(types)),
                          "max_record_len": max(lens) if lens else 0, "n_giant": len(gids),
                          "n_giant_single_caller": sum(1 for g in gids if ev.get(g, ("...",))[0].count("1") == 1),
                          "n_giant_bnd_notation": sum(1 for g in gids if ev.get(g, ("...", False))[1]),
                          "n_long_50k_500k": sum(1 for x in lens if CLUSTER_D < x <= MAX_CLUSTER)})

    summary, loci_by_v = {}, {}
    for vname, v in VARIANTS.items():
        if vname == "V0_current":
            clustered, dropped = v0, set()
        else:
            trows, dropped = transform(rows, ev, v)
            clustered = bt_cluster(trows, CLUSTER_D)
        res = classify(clustered, rows_len, by_id, order, dropped)
        over = [(c["chrom"], c["start"], c["end"]) for _, c, x in res if x is None]
        loci, n_cl = {}, 0
        for cid, c, x in res:
            if x is not None and x["is_cxsv"]:
                n_cl += 1
                key = frozenset(i for i in c["ids"] if i not in dropped)
                loci[key] = (c["chrom"], c["start"], c["end"], x["class"], x["n_junctions"], c["maxlen"],
                             ",".join(x["criteria"]))
        loci_by_v[vname] = loci
        resid = defaultdict(int)
        for ch, s, e in over:
            for cid, (c2, s2, e2) in over_iv.items():
                if ch == c2 and s < e2 and e > s2:
                    resid[cid] += min(e, e2) - max(s, s2)
        for orow in over_rows:
            orow[f"resid_{vname}"] = resid.get(orow["cluster"], 0)
        summary[vname] = {"n_clusters": len(res), "n_oversized": len(over), "blind_span": merged_span(over),
                          "largest_oversized": max([e - s for _, s, e in over], default=0),
                          "n_cxsv": len(loci), "n_cxsv_clusters": n_cl, "n_dropped_records": len(dropped)}
    base = loci_by_v["V0_current"]
    for vname, loci in loci_by_v.items():
        same = sum(1 for k in loci if k in base)
        summary[vname].update({"cxsv_same": same, "cxsv_new": len(loci) - same,
                               "cxsv_lost": sum(1 for k in base if k not in loci),
                               "cxsv_new_with_giant_member": sum(1 for k, v in loci.items()
                                                                 if k not in base and v[5] > MAX_CLUSTER)})
    for orow in over_rows:
        for vname, loci in loci_by_v.items():
            orow[f"cxsv_inside_{vname}"] = sum(1 for v in loci.values()
                                               if v[0] == orow["chrom"] and v[1] < orow["end"] and v[2] > orow["start"])

    with open(f"{odir}/summary.json", "w") as fh:
        json.dump({"cohort": cohort, "sample": sample, **meta, "v0_partition_identical": identical,
                   "n_rows": len(rows), "n_released_clusters": len(rel), "variants": summary,
                   "seconds": round(time.time() - t0, 1)}, fh)
    cols = list(over_rows[0].keys()) if over_rows else ["cohort"]
    with open(f"{odir}/oversized_clusters.tsv", "w") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(over_rows)
    gcols = ["cohort", "sample", "platform", "coverage_x", "id", "chrom", "start", "end", "len", "svtype",
             "supp_vec", "n_callers", "bnd_alt", "qual"]
    extra = sorted(({k for d in detail.values() if d for k in d} | {k for d in cev.values() for k in d})
                   - set(gcols))
    with open(f"{odir}/giant_records.tsv", "w") as fh:
        fh.write("\t".join(gcols + extra) + "\n")
        for r in rows:
            d = detail.get(r[4])
            if d:
                allv = {**d, **cev.get(r[4], {})}
                b = [cohort, sample, meta["platform"], meta["coverage_x"], r[4], r[0], r[1], r[2], r[2] - r[1],
                     r[3], d["supp_vec"], d["supp_vec"].count("1"), int(d["bnd_alt"]), d["qual"]]
                fh.write("\t".join(map(str, b + [allv.get(k, ".") for k in extra])) + "\n")
    with open(f"{odir}/loci_by_variant.tsv", "w") as fh:
        fh.write("variant\tchrom\tstart\tend\tclass\tn_junctions\tmax_member_len\tcriteria\tin_current\n")
        for vname, loci in loci_by_v.items():
            for k, v in loci.items():
                if vname != "V0_current" and k in base:
                    continue          # unchanged loci are written once, under V0
                fh.write("\t".join(map(str, (vname,) + v + (int(k in base),))) + "\n")
        for k, v in base.items():     # current loci each rule loses
            for vname, loci in loci_by_v.items():
                if k not in loci:
                    fh.write("\t".join(map(str, (f"LOST:{vname}",) + v + (1,))) + "\n")
    open(f"{odir}/done.ok", "w").close()
    return (cohort, sample, "ok")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--samples-root", required=True)
    ap.add_argument("--code", required=True, help="repository checkout (workflow/lib)")
    ap.add_argument("--fai", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--shard", default="0/1")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--samples", default="")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    sys.path.insert(0, f"{a.code}/workflow/lib")
    import cxsv_classes as cc
    import cxsv_junctions as cj
    G.update(a=a, cc=cc, cj=cj)     # inherited by the forked workers
    i, k = map(int, a.shard.split("/"))
    want = set(a.samples.split(",")) if a.samples else None
    tasks = []
    for n, r in enumerate(read_tsv(a.manifest)):
        if (want is not None and r["sample_id"] not in want) or (want is None and n % k != i):
            continue
        tasks.append((r["cohort"], r["sample_id"], {"platform": r["platform"], "coverage_x": r["coverage_x"],
                                                     "individual_id": r["individual_id"]}))
    ctx = mp.get_context("fork")
    with ctx.Pool(a.threads) as pool:
        for n, res in enumerate(pool.imap_unordered(one, tasks)):
            if n % 50 == 0:
                print(n, res, time.strftime("%H:%M:%S"), flush=True)
    print("done", len(tasks), flush=True)


if __name__ == "__main__":
    main()
