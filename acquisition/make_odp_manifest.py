
import csv
D="/mnt/sargasso1/CxSV/datasets"
ena={r["run_accession"]:r for r in csv.DictReader(open(f"{D}/manifests/manifest_PRJEB89727.tsv"),delimiter="\t")}
probe=list(csv.DictReader(open(f"{D}/manifests/odp_probe.tsv"),delimiter="\t"))
rows=[]; fallback=[]
for p in probe:
    e=ena[p["run_accession"]]
    rec=dict(run_accession=p["run_accession"], individual_id=e["individual_id"],
             filename=p["run_accession"]+".sra", url=p["url"], md5="",
             size_bytes=p["size_bytes"], format="SRA",
             est_coverage_x=e["est_coverage_x"], coverage_source="ena_base_count",
             passes_coverage_floor=e["passes_coverage_floor"],
             ena_fallback_url=e["url"], ena_md5=e["md5"], ena_bytes=e["size_bytes"])
    (rows if p["status"]=="ok" else fallback).append(rec)
for name,rs in [("manifest_PRJEB89727_odp.tsv",rows),("manifest_PRJEB89727_ena_fallback.tsv",fallback)]:
    with open(f"{D}/manifests/{name}","w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0].keys()),delimiter="\t"); w.writeheader(); w.writerows(rs)
    print(name, len(rs), "%.2f TB"%(sum(int(r["size_bytes"] or r["ena_bytes"]) for r in rs)/1e12))
print("passing floor in odp set:", sum(1 for r in rows if r["passes_coverage_floor"]=="yes"))
