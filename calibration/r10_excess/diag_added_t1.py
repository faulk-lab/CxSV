#!/usr/bin/env python3
"""Diagnose Tier-1 calls added by panel-v1.1 (roadmap item 13 follow-up, R10 excess).

For every panel genome: v1.1 LOFO Tier-1 loci (panel_rates.sh scores) minus the
genome's v1.0 LOFO Tier-1 loci (item-10 calibration, same Locus_ID) = added
Tier-1 loci. For each added locus: class, privacy state and f
(Panel_Junctions_Seen / Panel_Junctions), size, repeat / segdup / satellite
context, whether it lies in a v1.0 FILTERED_oversized_cluster span (formerly
hidden), and per junction: callers (SUPP_VEC), variant reads per caller (merged
VCF FORMAT DR, second value) relative to the genome's coverage, panel carriers.
Also per-genome call-set statistics (all genomes) to see why some genomes were
already noisy in v1.0.

    diag_added_t1.py --v10 BUILD10 --v11 BUILD11 --cal CAL_ROOT --old-scores OLD.tsv.gz
        --satellite BED --out-prefix PREFIX [--threads N]
"""
import argparse
import bisect
import csv
import gzip
import os
import sys
from collections import defaultdict
from multiprocessing import Pool

T1 = "Tier1_high_priority"
CALLERS = ("cuteSV", "SVIM", "Sniffles2")     # SURVIVOR file order (SUPP_VEC)


def read_bed(path):
    d = defaultdict(list)
    for line in open(path):
        f = line.split("\t")
        d[f[0]].append((int(f[1]), int(f[2])))
    return {c: sorted(v) for c, v in d.items()}


def overlap_bp(iv, chrom, s, e):
    v = iv.get(chrom)
    if not v:
        return 0
    i = bisect.bisect_left(v, (s, -1)) - 1
    tot = 0
    for a, b in v[max(i, 0):]:
        if a >= e:
            break
        tot += max(0, min(b, e) - max(a, s))
    return tot


def merged_records(path):
    """id -> (svtype, supp_vec, variant reads per caller) for every merged record; plus genome stats."""
    rec, st = {}, defaultdict(int)
    with open(path) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            info = dict(x.split("=", 1) for x in f[7].split(";") if "=" in x)
            sv, sv_vec = info.get("SVTYPE", "NA"), info.get("SUPP_VEC", "000")
            fmt = f[8].split(":")
            k = fmt.index("DR") if "DR" in fmt else None
            alt = []
            for col in f[9:12]:
                v = col.split(":")
                try:
                    alt.append(int(v[k].split(",")[1]) if k is not None else None)
                except (ValueError, IndexError):
                    alt.append(None)
            n = sv_vec.count("1")
            st["records"] += 1
            st[f"records_{sv}"] += 1
            st[f"callers_{n}"] += 1
            if sv in ("TRA", "BND") or (sv == "INV" and n == 1):
                st[f"bndlike_callers_{n}"] += 1
            if n == 1:
                c1 = CALLERS[sv_vec.index('1')]
                st[f"single_{c1}"] += 1
                if sv == "TRA":
                    try:
                        q = float(f[5])
                    except ValueError:
                        q = -1
                    qb = "q<10" if q < 10 else ("q10-19" if q < 20 else "q>=20")
                    st[f"single_{c1}_TRA_{qb}"] += 1
            ids = [col.split(":")[fmt.index("ID")] for col in f[9:12]] if "ID" in fmt else [f[2]]
            rec[f[2]] = (sv, sv_vec, alt, ids, f[5])
    return rec, st


def caller_support(p, need):
    """record id -> (QUAL, supporting reads) from the three filtered caller VCFs (only ids in need)."""
    out = {}
    for fn, keys in (("sniffles_filtered.vcf", ("SUPPORT",)), ("cutesv_filtered.vcf", ("RE",)),
                     ("svim_filtered.vcf", ("SUPPORT",))):
        path = f"{p}/callers/{fn}"
        if not os.path.exists(path):
            continue
        with open(path) as fh:
            for line in fh:
                if line[0] == "#":
                    continue
                f = line.split("\t", 8)
                if f[2] not in need:
                    continue
                info = dict(x.split("=", 1) for x in f[7].split(";") if "=" in x)
                v = next((info[k] for k in keys if k in info), None)
                try:
                    q = float(f[5])
                except ValueError:
                    q = None
                out[f[2]] = (q, int(v) if v and v.isdigit() else None)
    return out


def one(task):
    a, cohort, sample, cov, group, old_t1 = task
    p10 = f"{a.v10}/results/panel/samples/{cohort}/{sample}"
    p11 = f"{a.v11}/results/panel/samples/{cohort}/{sample}"
    rec, st = merged_records(f"{p11}/merged_sv.vcf")
    st.update(cohort=cohort, sample_id=sample, coverage_x=cov, group=group)
    sc = list(csv.DictReader(open(f"{a.cal}/persample/{cohort}/{sample}/scores_LOFO.tsv"), delimiter="\t"))
    st["T1_v1.1"] = sum(r["Tier"] == T1 for r in sc)
    st["T1_v1.0"] = len(old_t1)
    added = [r for r in sc if r["Tier"] == T1 and r["Locus_ID"] not in old_t1]
    st["T1_added"] = len(added)
    if a.all_t1:
        added = [r for r in sc if r["Tier"] in (T1, "Tier2")]
    st["T1_lost"] = len(old_t1 - {r["Locus_ID"] for r in sc if r["Tier"] == T1})
    if not added:
        return st, []
    status = {f"{r['Chr']}:{r['Start']}-{r['End']}": r for r in
              csv.DictReader(open(f"{a.cal}/leaveout/{cohort}/{sample}/status_LOFO.tsv"), delimiter="\t")}
    jst = defaultdict(list)
    for r in csv.DictReader(open(f"{a.cal}/leaveout/{cohort}/{sample}/jstatus_LOFO.tsv"), delimiter="\t"):
        for lid in r["locus_id"].split(","):
            jst[lid].append(r)
    jun = defaultdict(list)
    for r in csv.DictReader(open(f"{p11}/junctions/junctions.tsv"), delimiter="\t"):
        for lid in r["locus_id"].split(","):
            jun[lid].append(r)
    over = defaultdict(list)
    for r in csv.DictReader(open(f"{p10}/cxsv/cxsv_summary.tsv"), delimiter="\t"):
        if "FILTERED_oversized_cluster" in r["CxSV_Criteria"]:
            over[r["Chr"]].append((int(r["Start"]), int(r["End"])))
    v10_ids = set()
    with open(f"{p10}/cxsv/cxsv_master_table.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            v10_ids.add(f"{r['Chr']}:{r['Start']}-{r['End']}")
    need = set()
    for r in added:
        for j in jun.get(r["Locus_ID"], []):
            for rid in j["record_ids"].split(","):
                x = rec.get(rid)
                if x:
                    need.update(i for i in x[3] if i not in ("NaN", "NA", "."))
    sup = caller_support(p11, need)
    out = []
    for r in added:
        L = r["Locus_ID"]
        s0 = status.get(L, {})
        chrom, s, e = r["Chr"], int(r["Start"]), int(r["End"])
        js = jun.get(L, [])
        ncall, alts, svt, svimq, svim_only, jd = [], [], defaultdict(int), [], 0, []
        for j in js:
            jq = []
            ncall.append(j["supp_vec"].count("1"))
            svt[j["svtype"]] += 1
            svim_only += j["supp_vec"] == "010"
            best = 0
            for rid in j["record_ids"].split(","):
                x = rec.get(rid)
                if x:
                    for i in x[3]:
                        if i in sup:
                            best = max(best, sup[i][1] or 0)
                            if i.startswith("svim."):
                                svimq.append(sup[i][0])
                                jq.append(sup[i][0])
            alts.append(best)
            jd.append(f"{j['svtype']}|{j['supp_vec']}|{best}|{max(jq) if jq else 'NA'}")
        carriers = [int(x["panel_carriers"]) if x["panel_carriers"] not in ("", "NA", ".") else 0 for x in jst.get(L, [])]
        pj, ps = s0.get("Panel_Junctions", "NA"), s0.get("Panel_Junctions_Seen", "NA")
        f = (int(ps) / int(pj)) if pj not in ("NA", "0", "") else None
        sat = overlap_bp(a.sat, chrom, s, e)
        hid = overlap_bp({chrom: sorted(over[chrom])}, chrom, s, e) if over.get(chrom) else 0
        covf = float(cov) if cov not in ("", "NA") else None
        out.append({
            "cohort": cohort, "sample_id": sample, "group": group, "coverage_x": cov, "Locus_ID": L, "Size": r["Size"],
            "Tier": r["Tier"], "added_T1": "yes" if (r["Tier"] == T1 and L not in old_t1) else "no",
            "junction_detail": ";".join(jd),
            "Cmplx_Class": r["Cmplx_Class"], "SV_Types": r["SV_Types"], "Panel_Status": r["Panel_Status"],
            "Panel_Junctions": pj, "Panel_Junctions_Seen": ps, "f_seen": "NA" if f is None else f"{f:.3f}",
            "Total_Score": r["Total_Score"], "Top_Features": r["Top_Features"],
            "in_v1.0_loci": "yes" if L in v10_ids else "no",
            "v1.0_oversized_overlap_frac": f"{hid / max(1, e - s):.3f}",
            "N_Junctions": len(js), "junction_svtypes": ",".join(f"{k}:{v}" for k, v in sorted(svt.items())),
            "junctions_1caller": sum(n == 1 for n in ncall), "junctions_2plus_callers": sum(n >= 2 for n in ncall),
            "max_callers": max(ncall) if ncall else 0,
            "median_alt_reads": sorted(alts)[len(alts) // 2] if alts else "NA",
            "median_alt_over_cov": f"{sorted(alts)[len(alts) // 2] / covf:.3f}" if alts and covf else "NA",
            "min_alt_reads": min(alts) if alts else "NA",
            "junctions_svim_only": svim_only,
            "median_svim_qual": sorted(svimq)[len(svimq) // 2] if svimq else "NA",
            "junctions_panel_carriers_0": sum(c == 0 for c in carriers),
            "max_panel_carriers": max(carriers) if carriers else "NA",
            "satellite_frac": f"{sat / max(1, e - s):.3f}", "SegDup_Overlap": s0.get("SegDup_Overlap", "NA"),
            "Genomic_Context": s0.get("Genomic_Context", "NA"), "Repeat_Annotation": s0.get("Repeat_Annotation", "NA"),
            "read_support_value": r.get("read_support_value", "NA"), "caller_agreement_value": r.get("caller_agreement_value", "NA"),
            "artefact_context_value": r.get("artefact_context_value", "NA"),
        })
    return st, out


def main():
    ap = argparse.ArgumentParser()
    for k in ("--v10", "--v11", "--cal", "--old-scores", "--satellite", "--out-prefix", "--per-genome-old"):
        ap.add_argument(k)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--all-t1", action="store_true", help="describe every Tier-1 and Tier-2 locus, not only added Tier 1")
    a = ap.parse_args()
    a.sat = read_bed(a.satellite)
    old = defaultdict(set)
    with gzip.open(a.old_scores, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["mode"] == "LOFO" and r["Tier"] == T1:
                old[(r["cohort"], r["sample_id"])].add(r["Locus_ID"])
    grp = {(r["cohort"], r["sample_id"]): r["group"] for r in csv.DictReader(open(a.per_genome_old), delimiter="\t")}
    man = list(csv.DictReader(open(f"{a.v11}/results/panel/panel/sample_manifest.tsv"), delimiter="\t"))
    tasks = [(a, m["cohort"], m["sample_id"], m["coverage_x"], grp[(m["cohort"], m["sample_id"])],
              old[(m["cohort"], m["sample_id"])]) for m in man]
    sts, loci = [], []
    with Pool(a.threads) as p:
        for i, (st, out) in enumerate(p.imap_unordered(one, tasks, chunksize=4)):
            sts.append(st)
            loci += out
            if i % 200 == 0:
                print(i, file=sys.stderr)
    keys = ["cohort", "sample_id", "group", "coverage_x", "T1_v1.0", "T1_v1.1", "T1_added", "T1_lost"]
    keys += sorted({k for s in sts for k in s} - set(keys))
    with open(a.out_prefix + "_genomes.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t", lineterminator="\n", restval=0)
        w.writeheader()
        w.writerows(sts)
    with open(a.out_prefix + "_added_loci.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(loci[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(loci)
    print(f"genomes {len(sts)}, added Tier-1 loci {len(loci)}")


if __name__ == "__main__":
    main()
