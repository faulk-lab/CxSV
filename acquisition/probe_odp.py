
import csv, concurrent.futures as cf, urllib.request, sys
man="/mnt/sargasso1/CxSV/datasets/manifests/manifest_PRJEB89727.tsv"
rows=list(csv.DictReader(open(man),delimiter="\t"))
print("runs:",len(rows),flush=True)
def check(r):
    acc=r["run_accession"]
    u=f"https://sra-pub-run-odp.s3.amazonaws.com/sra/{acc}/{acc}"
    req=urllib.request.Request(u, headers={"Range":"bytes=0-0"})
    try:
        with urllib.request.urlopen(req,timeout=60) as resp:
            cr=resp.headers.get("Content-Range","")
            size=int(cr.split("/")[-1]) if "/" in cr else 0
            return acc, u, size, "ok"
    except Exception as e:
        return acc, u, 0, f"{type(e).__name__}:{getattr(e,'code','')}"
out=[]
with cf.ThreadPoolExecutor(32) as ex:
    for i,res in enumerate(ex.map(check,rows),1):
        out.append(res)
        if i%200==0: print(i,res,flush=True)
w=csv.writer(open("/mnt/sargasso1/CxSV/datasets/manifests/odp_probe.tsv","w",newline=""),delimiter="\t")
w.writerow(["run_accession","url","size_bytes","status"]); w.writerows(out)
ok=[r for r in out if r[3]=="ok"]
print("available:",len(ok),"of",len(out))
print("total TB: %.2f"%(sum(r[2] for r in ok)/1e12))
bad=[r for r in out if r[3]!="ok"]
print("missing:",len(bad), bad[:5])
