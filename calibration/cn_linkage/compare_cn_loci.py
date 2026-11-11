#!/usr/bin/env python3
"""Re-run cn_linkage.py (this branch) on finished patient runs and compare CN loci with the
run's own cn/cn_loci.tsv (panel-v1.2 code). Two variants: 'locus' = item g only (v1.2 public
rule: public_min_genomes_sexchrom 1, no discordance exclusion); 'v13' = item g + D2 defaults.
Usage: compare_cn_loci.py REPO RELEASE FAI OUTDIR SAMPLE_DIR [SAMPLE_DIR ...]"""
import csv, json, os, subprocess, sys
repo, rel, fai, out = sys.argv[1:5]
dirs = sys.argv[5:]
def rt(p):
    with open(p) as fh:
        return list(csv.DictReader([l for l in fh if not l.startswith("#")], delimiter="\t"))
rows = []
for S in dirs:
    name = os.path.basename(S.rstrip("/"))
    stats = json.load(open(f"{S}/cn/stats.json"))
    base = json.load(open(os.path.join(S, "..", "..", "..", "cn_params.json")))
    old = rt(f"{S}/cn/cn_loci.tsv")
    for var in ("locus", "v13"):
        W = f"{out}/{name}/{var}"; os.makedirs(W, exist_ok=True)
        p = json.loads(json.dumps(base))
        if var == "locus":
            p["cn_params"]["linkage"].update(public_min_genomes_sexchrom=1, sexchrom_discordance_x=99, sexchrom_discordance_y=99)
            mem = "none"
        else:
            p["cn_params"]["linkage"].update(public_min_genomes_sexchrom=2, sexchrom_discordance_x=0.10, sexchrom_discordance_y=0.10)
            mem = f"{rel}/cn/pon/pon_members.tsv"
        json.dump(p, open(f"{W}/params.json", "w"))
        cmd = [sys.executable, f"{repo}/workflow/scripts/cn_linkage.py", "--vcf", f"{S}/merged_sv.vcf",
               "--clusters", f"{S}/cxsv/junction_clusters.bed", "--summary", f"{S}/cxsv/junction_summary.tsv",
               "--cxsv-vcf", f"{S}/cxsv/junction_only.vcf", "--single", f"{S}/cxsv/single_caller_filtered.tsv",
               "--segments", f"{S}/cn/segments.tsv", "--panel-segments", f"{rel}/cn/panel_cn_segments.tsv",
               "--pon-members", mem, "--individual", name, "--stats", f"{S}/cn/stats.json", "--fai", fai,
               "--calling", f"{W}/params.json", "--scripts", f"{repo}/workflow/scripts", "--work", f"{W}/work",
               "--out-clusters", f"{W}/clusters.bed", "--out-summary", f"{W}/cxsv_summary.tsv", "--out-vcf", f"{W}/cxsv_only.vcf",
               "--out-candidates", f"{W}/cn_candidates.tsv", "--out-loci", f"{W}/cn_loci.tsv", "--out-cn-only", f"{W}/cn_only.tsv"]
        r = subprocess.run(cmd, capture_output=True, text=True)
        open(f"{W}/stderr.log", "w").write(r.stderr)
        if r.returncode:
            rows.append(dict(sample=name, variant=var, error=r.stderr[-400:])); continue
        new = rt(f"{W}/cn_loci.tsv")
        def zero(L): return [x["locus_id"] for x in L if int(x["locus_id"].rsplit("-",1)[1]) - int(x["locus_id"].rsplit(":",1)[1].split("-")[0]) <= 2]
        def noseg(L): return [x["locus_id"] for x in L if not x["segments"]]
        rows.append(dict(sample=name, variant=var, sex_call=stats.get("sex_call"),
                         old_loci=len(old), old_written=sum(x["written"]=="yes" for x in old),
                         old_2bp=len(zero(old)), old_no_segment=len(noseg(old)),
                         new_loci=len(new), new_written=sum(x["written"]=="yes" for x in new),
                         new_2bp=len(zero(new)), new_no_segment=len(noseg(new)),
                         old_candidates_locus=sum(1 for x in rt(f"{S}/cn/cn_candidates.tsv") if x["decision"]=="locus"),
                         new_candidates_locus=sum(1 for x in rt(f"{W}/cn_candidates.tsv") if x["decision"]=="locus"),
                         new_ids=";".join(x["locus_id"]+"="+x["cn_class"] for x in new),
                         old_ids=";".join(x["locus_id"]+"="+x["cn_class"] for x in old)))
keys = sorted({k for r in rows for k in r}, key=lambda k: ["sample","variant","error"].index(k) if k in ("sample","variant","error") else 9)
with open(f"{out}/compare_cn_loci.tsv", "w") as fh:
    w = csv.DictWriter(fh, keys, delimiter="\t", extrasaction="ignore"); w.writeheader(); w.writerows(rows)
for r in rows:
    print({k: v for k, v in r.items() if k not in ("new_ids", "old_ids")})
