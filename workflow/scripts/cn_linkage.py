#!/usr/bin/env python3
"""Copy-number-linked loci of one sample: the second discovery track (rule
sample_cn_linkage, workflow/rules/cn.smk; docs/CN_LINKAGE.md). It adds loci and
removes none.

  cn_linkage.py --vcf merged_sv.vcf --clusters cxsv/junction_clusters.bed \\
      --summary cxsv/junction_summary.tsv --cxsv-vcf cxsv/junction_only.vcf \\
      --single cxsv/single_caller_filtered.tsv --segments cn/segments.tsv \\
      --panel-segments panel_cn_segments.tsv --individual ID --fai ref.fa.fai \\
      --calling calling_params.json --scripts workflow/scripts --work cn/work \\
      --out-clusters cxsv/clusters.bed --out-summary cxsv/cxsv_summary.tsv \\
      --out-vcf cxsv/cxsv_only.vcf --out-candidates cn/cn_candidates.tsv \\
      --out-loci cn/cn_loci.tsv --out-cn-only cn/cn_only.tsv

1. Candidates (cxsv_cnlink.candidates): abnormal segments grouped within max_gap,
   every merged_sv.vcf record attached by breakend position (distance-free; the records
   discovery sets aside as long spans or svim_lowq included), structure class, the
   complexity rule and the segment-level panel frequency (segments of every OTHER
   individual in --panel-segments; a candidate whose every segment is seen there is
   public_cn and gets no locus).
2. A candidate with decision "locus" becomes one locus per contig: its records go
   through vcf_to_bed.py (--span-link-max 0, every record kept) and classify_cxsv.py
   (--max-cluster-size 0) exactly as discovery runs them; the row is then IsCxSV=YES,
   CxSV_Criteria "CN_linked;<criteria>", Cmplx_Class = the structure class. A locus
   whose id (Chr:Start-End) equals a junction-track locus is not added.
3. clusters.bed / cxsv_summary.tsv = junction track + CN loci; cxsv_only.vcf =
   junction track + the records of the CN loci. Every candidate is listed in
   cn_candidates.tsv; the CN loci with their segments in cn_loci.tsv (read by
   locus_segments.py for depth at copy-number boundaries); candidates with no attached
   junction in cn_only.tsv (the patient report lists them, untiered).
"""
import argparse
import csv
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import cxsv_cnlink as cl  # noqa: E402

ID_BASE = 900000            # cluster ids of the CN loci (junction-track ids are small integers)
SEG_KEYS = ("start_mask_before", "start_mask_after", "end_mask_before", "end_mask_after")


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader([l for l in fh if not l.startswith("#")], delimiter="\t"))


def read_segments(path):
    out = []
    for r in read_tsv(path):
        s = dict(chrom=r["chrom"], start=int(r["start"]), end=int(r["end"]), state=int(r["state"]),
                 ploidy=int(r["ploidy"]), cn=float(r.get("cn") or 0), mean_ratio=float(r.get("mean_ratio") or 0),
                 n_bins=int(float(r.get("n_bins") or 0)), z=float(r.get("z") or 0), nested=r.get("nested") == "True")
        for k in SEG_KEYS:
            s[k] = int(float(r.get(k) or 0))
        out.append(s)
    return out


def others(path, individual):
    by = {}
    if not path or path == "none" or not os.path.exists(path):
        return []
    for r in read_tsv(path):
        if r["individual_id"] == individual:
            continue
        by.setdefault((r["key"], r["individual_id"]), []).append(
            dict(chrom=r["chrom"], start=int(r["start"]), end=int(r["end"]), state=int(r["state"]), ploidy=int(r["ploidy"])))
    return [(k, ind, s) for (k, ind), s in by.items()]


def sh(cmd, log):
    with open(log, "a") as fh:
        fh.write("$ " + " ".join(cmd) + "\n")
        fh.flush()
        subprocess.run(cmd, check=True, stdout=fh, stderr=fh)


def seg_str(segs):
    return ";".join(f"{s['chrom']}:{s['start']}-{s['end']}:{s['state']}/{s['ploidy']}" for s in segs)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    for k in ("vcf", "clusters", "summary", "cxsv-vcf", "single", "segments", "fai", "calling", "scripts", "work",
              "out-clusters", "out-summary", "out-vcf", "out-candidates", "out-loci", "out-cn-only"):
        ap.add_argument(f"--{k}", required=True)
    ap.add_argument("--panel-segments", default="none")
    ap.add_argument("--individual", default="")
    a = ap.parse_args()
    cfg = json.load(open(a.calling))
    c, j = cfg["cxsv_params"], cfg["junction_params"]
    P = {**cl.CN_DEFAULTS, **(cfg.get("cn_params", {}).get("linkage") or {}), "span_link_max": c["span_link_max"]}
    os.makedirs(a.work, exist_ok=True)
    log = f"{a.work}/cn_linkage.log"
    open(log, "w").close()
    segs = read_segments(a.segments)
    recs = cl.read_records(a.vcf)
    single = {r["ID"] for r in read_tsv(a.single)} if os.path.exists(a.single) else set()
    cands = cl.candidates(segs, recs, others(a.panel_segments, a.individual), P, single_filtered=single)
    summ0 = read_tsv(a.summary)
    header = open(a.summary).readline().rstrip("\n").split("\t")
    existing = {f"{r['Chr']}:{r['Start']}-{r['End']}" for r in summ0 if r.get("IsCxSV") == "YES"}
    # 2. one locus per contig of every "locus" candidate
    cn_rows, meta, nid = [], {}, ID_BASE
    for k in [x for x in cands if x["decision"] == "locus"]:
        mini, bed = f"{a.work}/cand{k['id']}.vcf", f"{a.work}/cand{k['id']}.bed"
        ids = {r["id"] for r in k["records"]}
        with open(a.vcf) as fi, open(mini, "w") as fo:
            for line in fi:
                if line.startswith("#") or line.split("\t", 3)[2] in ids:
                    fo.write(line)
        sh([sys.executable, f"{a.scripts}/vcf_to_bed.py", "--vcf", mini, "--out", bed, "--window", str(c["bnd_window"]),
            "--span-link-max", "0", "--long-span-rule", "exclude", "--long-out", os.devnull,
            "--single-caller-bnd-rule", "keep", "--svim-min-qual", str(c.get("svim_min_qual_discovery", 10)),
            "--filtered-out", os.devnull], log)
        by_chr = {}
        for line in open(bed):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 6:
                by_chr.setdefault(f[0], []).append(f[:6])
        for chrom, rows in by_chr.items():
            nid += 1
            for f in sorted(rows, key=lambda x: int(x[1])):
                cn_rows.append(f + [str(nid)])
            meta[str(nid)] = (k, chrom)
    written, keep = [], set()
    if cn_rows:
        cn_bed, cn_sum = f"{a.work}/cn_clusters.bed", f"{a.work}/cn_summary.tsv"
        with open(cn_bed, "w") as fh:
            fh.writelines("\t".join(f) + "\n" for f in cn_rows)
        sh([sys.executable, f"{a.scripts}/classify_cxsv.py", "--bed", cn_bed, "--vcf", a.vcf, "--fai", a.fai,
            "--output", cn_sum, "--collapse-distance", str(j["collapse_distance"]),
            "--ins-len-ratio", str(j["ins_len_ratio"]), "--composite-fraction", str(c["composite_fraction"]),
            "--min-bp", str(c["min_breakpoints"]), "--max-window", str(c["cluster_distance"]),
            "--chromoanag-min", str(c["chromoanag_min_window"]), "--chromoanag-max", str(c["chromoanag_max_window"]),
            "--chromoanag-min-sv", str(c["chromoanag_min_svs"]), "--extreme-min-breakends", str(c["extreme_min_breakends"]),
            "--min-interleaved", str(c["min_interleaved"]), "--interleave-alt-min-intra", str(c["interleave_alt_min_intra"]),
            "--interleave-alt-min-inter", str(c["interleave_alt_min_inter"]), "--max-cluster-size", "0"], log)
        for r in read_tsv(cn_sum):
            k, chrom = meta[r["ClusterID"]]
            lid = f"{r['Chr']}:{r['Start']}-{r['End']}"
            row = dict(locus_id=lid, candidate=k["id"], cn_class=k["cls"], classifier_class=r.get("Cmplx_Class", "NA"),
                       n_records=len(k["records"]), n_long=k["n_long"], n_single_filtered=k["n_single_filtered"],
                       segments=seg_str([s for s in k["segs"] if s["chrom"] == chrom]))
            if lid in existing:
                row["written"] = "no: same id as a junction-track locus"
            else:
                r["IsCxSV"] = "YES"
                r["CxSV_Criteria"] = "CN_linked;" + (r.get("CxSV_Criteria") or "")
                r["Cmplx_Class"] = k["cls"].replace("/", "_")
                keep.add(r["ClusterID"])
                existing.add(lid)
                row["written"] = "yes"
                summ0.append(r)
            written.append(row)
    # 3. combined discovery outputs
    with open(a.out_clusters, "w") as fo:
        fo.write(open(a.clusters).read())
        fo.writelines("\t".join(f) + "\n" for f in cn_rows if f[6] in keep)
    with open(a.out_summary, "w", newline="") as fo:
        w = csv.DictWriter(fo, header, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(summ0)
    lines = open(a.cxsv_vcf).readlines()
    head, body = [l for l in lines if l.startswith("#")], [l for l in lines if not l.startswith("#")]
    have = {l.split("\t", 3)[2] for l in body}
    for cid in sorted(keep):
        for r in meta[cid][0]["records"]:
            if r["id"] not in have:
                have.add(r["id"])
                body.append(r["line"])
    order = {l.split("\t")[0]: i for i, l in enumerate(open(a.fai))}
    body.sort(key=lambda l: (order.get(l.split("\t", 1)[0], 1 << 30), int(l.split("\t", 2)[1])))
    with open(a.out_vcf, "w") as fo:
        fo.writelines(head + body)
    cols = ["candidate", "decision", "cls", "nested", "inverted", "parts", "n_segments", "segments", "n_records",
            "record_types", "n_long", "n_single_filtered", "public_counts", "span"]
    with open(a.out_candidates, "w") as fh, open(a.out_cn_only, "w") as fo:
        fh.write("\t".join(cols) + "\n")
        fo.write("\t".join(cols) + "\n")
        for k in cands:
            row = [k["id"], k["decision"], k["cls"], k["nested"], k["inverted"],
                   ";".join(f"{x}:{b}-{e}" for x, b, e in k["parts"]), len(k["segs"]), seg_str(k["segs"]),
                   len(k["records"]), ",".join(sorted(r["svtype"] for r in k["records"])), k["n_long"],
                   k["n_single_filtered"], ",".join(map(str, k["public_counts"])), k["span"]]
            line = "\t".join(map(str, row)) + "\n"
            fh.write(line)
            if k["decision"] == "cn_only":
                fo.write(line)
    lcols = ["locus_id", "candidate", "cn_class", "classifier_class", "n_records", "n_long", "n_single_filtered",
             "segments", "written"]
    with open(a.out_loci, "w") as fh:
        fh.write("\t".join(lcols) + "\n")
        for r in written:
            fh.write("\t".join(str(r[c]) for c in lcols) + "\n")
    n = {d: sum(k["decision"] == d for k in cands) for d in ("locus", "cn_only", "simple_cnv", "public_cn")}
    print(f"{len(segs)} segments, {len(cands)} candidates {n}; CN loci written {sum(r['written'] == 'yes' for r in written)}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
