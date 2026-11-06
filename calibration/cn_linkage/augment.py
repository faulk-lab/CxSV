"""Write a genome's discovery outputs with the copy-number-linked loci added
(roadmap item 15, Phase 1; second discovery track: it adds loci, removes nothing).

Starting from the genome's own discovery outputs (cxsv/clusters.bed, cxsv_summary.tsv,
cxsv_only.vcf of the panel build or the patient run), for every candidate with
decision "locus" (cnlink.candidates):

  1 the attached records' BED rows are made by workflow/scripts/vcf_to_bed.py itself
    (on a VCF of just those records, --span-link-max 0, every record kept), one new
    cluster per contig (a translocation's two ends are two loci, as in discovery);
  2 workflow/scripts/classify_cxsv.py classifies the new clusters alone
    (--max-cluster-size 0: a copy-number-linked locus may exceed 500 kb); each row is
    then set IsCxSV=YES, CxSV_Criteria "CN_linked;<criteria>", Cmplx_Class = the
    structure-model class;
  3 clusters.bed and cxsv_summary.tsv = original + new; extract_junctions.py,
    annotate_cxsv.py on them (same scripts and parameters as discovery);
    cxsv_only.vcf = original + the attached records (readmitted long / svim_lowq
    records included).
A new locus whose id (Chr:Start-End) equals an existing locus id is skipped (reported).
"""
import csv
import os
import subprocess

import yaml


def sh(cmd, log):
    with open(log, "a") as fh:
        fh.write("$ " + cmd + "\n")
        fh.flush()
        subprocess.run(cmd, shell=True, check=True, stdout=fh, stderr=fh, executable="/bin/bash")


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader([l for l in fh if not l.startswith("#")], delimiter="\t"))


def load_calling(code):
    return yaml.safe_load(open(f"{code}/config/calling.yaml"))


def write_mini_vcf(src_vcf, recs, out):
    ids = {r["id"] for r in recs}
    with open(src_vcf) as fi, open(out, "w") as fo:
        for line in fi:
            if line.startswith("#"):
                fo.write(line)
            elif line.split("\t", 3)[2] in ids:
                fo.write(line)


def augment(src_dir, out_dir, merged_vcf, cands, code, fai, log, id_base=900000):
    """Returns per-candidate rows describing the loci written."""
    cfg = load_calling(code)
    c, j = cfg["cxsv_params"], cfg["junction_params"]
    res = os.environ.get("CXSV_RESOURCES", "/mnt/sargasso1/CxSV/resources")
    an = {k: (v if str(v).startswith("/") else f"{res}/{os.path.basename(str(v))}")
          for k, v in cfg["annotation"].items() if k.endswith("_bed")}
    S = f"{code}/workflow/scripts"
    os.makedirs(f"{out_dir}/cxsv", exist_ok=True)
    os.makedirs(f"{out_dir}/junctions", exist_ok=True)
    os.makedirs(f"{out_dir}/cn", exist_ok=True)
    if not os.path.exists(f"{out_dir}/merged_sv.vcf"):      # read by scorers of whole runs (score_pos_control.py)
        try:
            os.link(merged_vcf, f"{out_dir}/merged_sv.vcf")
        except OSError:
            import shutil
            shutil.copy(merged_vcf, f"{out_dir}/merged_sv.vcf")
    summ0 = read_tsv(f"{src_dir}/cxsv/cxsv_summary.tsv")
    existing = {f"{r['Chr']}:{r['Start']}-{r['End']}" for r in summ0 if r["IsCxSV"] == "YES"}
    loci = [k for k in cands if k["decision"] == "locus"]
    report = []
    cn_rows = []
    cid_meta = {}
    nid = id_base
    for k in loci:
        mini = f"{out_dir}/cn/cand{k['id']}.vcf"
        write_mini_vcf(merged_vcf, k["records"], mini)
        bed = f"{out_dir}/cn/cand{k['id']}.bed"
        sh(f"python {S}/vcf_to_bed.py --vcf {mini} --out {bed} --window {c['bnd_window']} --span-link-max 0 "
           f"--long-span-rule exclude --long-out /dev/null --single-caller-bnd-rule keep "
           f"--svim-min-qual {c.get('svim_min_qual_discovery', 10)} --filtered-out /dev/null", log)
        by_chr = {}
        for line in open(bed):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 6:
                by_chr.setdefault(f[0], []).append(f[:6])
        for chrom, rows in by_chr.items():
            nid += 1
            for f in sorted(rows, key=lambda x: int(x[1])):
                cn_rows.append(f + [str(nid)])
            cid_meta[str(nid)] = k
    if not cn_rows:
        for name in ("clusters.bed", "cxsv_summary.tsv", "cxsv_only.vcf", "cxsv_master_table.tsv", "cxsv_loci.bed"):
            if os.path.exists(f"{src_dir}/cxsv/{name}"):
                sh(f"cp {src_dir}/cxsv/{name} {out_dir}/cxsv/{name}", log)
        sh(f"cp {src_dir}/junctions/junctions.tsv {out_dir}/junctions/junctions.tsv", log)
        return report
    cn_bed = f"{out_dir}/cn/cn_clusters.bed"
    with open(cn_bed, "w") as fh:
        for f in cn_rows:
            fh.write("\t".join(f) + "\n")
    cn_sum = f"{out_dir}/cn/cn_summary.tsv"
    sh(f"python {S}/classify_cxsv.py --bed {cn_bed} --vcf {merged_vcf} --fai {fai} --output {cn_sum} "
       f"--collapse-distance {j['collapse_distance']} --ins-len-ratio {j['ins_len_ratio']} "
       f"--composite-fraction {c['composite_fraction']} --min-bp {c['min_breakpoints']} --max-window {c['cluster_distance']} "
       f"--chromoanag-min {c['chromoanag_min_window']} --chromoanag-max {c['chromoanag_max_window']} "
       f"--chromoanag-min-sv {c['chromoanag_min_svs']} --extreme-min-breakends {c['extreme_min_breakends']} "
       f"--min-interleaved {c['min_interleaved']} --interleave-alt-min-intra {c['interleave_alt_min_intra']} "
       f"--interleave-alt-min-inter {c['interleave_alt_min_inter']} --max-cluster-size 0", log)
    rows = read_tsv(cn_sum)
    header = open(f"{src_dir}/cxsv/cxsv_summary.tsv").readline().rstrip("\n").split("\t")
    keep_ids = set()
    for r in rows:
        k = cid_meta[r["ClusterID"]]
        lid = f"{r['Chr']}:{r['Start']}-{r['End']}"
        rep = dict(candidate=k["id"], locus_id=lid, cls=k["cls"], classifier_class=r.get("Cmplx_Class", "NA"),
                   classifier_iscxsv=r["IsCxSV"], n_records=len(k["records"]), n_long=k["n_long"],
                   n_single_filtered=k["n_single_filtered"])
        if lid in existing:
            rep["written"] = "no: same id as an existing locus"
        else:
            r["IsCxSV"] = "YES"
            r["CxSV_Criteria"] = "CN_linked;" + (r.get("CxSV_Criteria") or "")
            r["Cmplx_Class"] = k["cls"].replace("/", "_")
            keep_ids.add(r["ClusterID"])
            existing.add(lid)
            rep["written"] = "yes"
        report.append(rep)
    with open(f"{out_dir}/cxsv/clusters.bed", "w") as fo:
        fo.write(open(f"{src_dir}/cxsv/clusters.bed").read())
        for f in cn_rows:
            if f[6] in keep_ids:
                fo.write("\t".join(f) + "\n")
    with open(f"{out_dir}/cxsv/cxsv_summary.tsv", "w", newline="") as fo:
        w = csv.DictWriter(fo, header, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(summ0)
        w.writerows([r for r in rows if r["ClusterID"] in keep_ids])
    # cxsv_only.vcf: original + attached records of the written loci
    have = set()
    lines = open(f"{src_dir}/cxsv/cxsv_only.vcf").readlines()
    head = [l for l in lines if l.startswith("#")]
    body = [l for l in lines if not l.startswith("#")]
    for l in body:
        have.add(l.split("\t", 3)[2])
    add = []
    for cid in keep_ids:
        for r in cid_meta[cid]["records"]:
            if r["id"] not in have:
                have.add(r["id"])
                add.append(r["line"])
    order = {l.split("\t")[0]: i for i, l in enumerate(open(fai))}
    body = sorted(body + add, key=lambda l: (order.get(l.split("\t", 1)[0], 999), int(l.split("\t", 2)[1])))
    with open(f"{out_dir}/cxsv/cxsv_only.vcf", "w") as fo:
        fo.writelines(head + body)
    sh(f"python {S}/extract_junctions.py --vcf {merged_vcf} --clusters {out_dir}/cxsv/clusters.bed "
       f"--summary {out_dir}/cxsv/cxsv_summary.tsv --fai {fai} --collapse-distance {j['collapse_distance']} "
       f"--ins-len-ratio {j['ins_len_ratio']} --out {out_dir}/junctions/junctions.tsv", log)
    sh(f"python {S}/annotate_cxsv.py --summary {out_dir}/cxsv/cxsv_summary.tsv --vcf {out_dir}/cxsv/cxsv_only.vcf "
       f"--genes {an['genes_bed']} --exons {an['exons_bed']} --repeats {an['repeats_bed']} "
       f"--segdups {an['segdup_bed']} --out-table {out_dir}/cxsv/cxsv_master_table.tsv --out-bed {out_dir}/cxsv/cxsv_loci.bed", log)
    return report
