#!/usr/bin/env python3
"""Which external SV catalogue, if any, contains each patient CxSV locus (roadmap item 4).

Every junction of every patient CxSV locus (junctions/junctions.tsv, locus_id)
is matched against each catalogue's junction-site table with the item-3
matcher (workflow/lib/cxsv_matching.py, SiteIndex.query): svtype and contigs
agree and both ends lie within --tolerance bp of the site (INS: lengths within
--ins-len-ratio). A catalogue *contains* a locus when at least
--min-junction-fraction of its junctions match a site there.

Before use, each table's sha256 must equal its manifest's output_sha256
(acquisition/catalogues manifests); a mismatch stops the run.

Catalogues overlap each other and the panel (some individuals are in several),
so their frequencies are reported side by side and never summed.

Outputs
-------
--out-loci       one row per locus: Locus_ID, then per catalogue <C>_Junctions_Matched
                 (k/n), <C>_Max_AF, <C>_Max_Group_AF, <C>_Contains; then
                 Catalogues_Containing and Catalogue_Max_AF (over containing catalogues)
--out-junctions  one row per (locus, junction): coordinates, supp_vec, and per
                 catalogue <C>_sites, <C>_max_af, <C>_best_site
--out-check      the verification record (catalogue, version, table, sha256, status)

Usage: match_catalogues.py --junctions J --loci cxsv_panel_status.tsv
           [--catalogue NAME:VERSION:TABLE:MANIFEST ...] --out-loci ... --out-junctions ... --out-check ...
"""
import argparse
import csv
import hashlib
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import cxsv_junctions as cj  # noqa: E402
import cxsv_matching as cm   # noqa: E402


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def manifest_field(path, field):
    with open(path, newline="") as fh:
        for r in csv.reader(fh, delimiter="\t"):
            if len(r) >= 2 and r[0] == field:
                return r[1]
    return None


def fmt(x):
    return "." if x is None else f"{x:.6g}"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--junctions", required=True)
    ap.add_argument("--loci", required=True, help="cxsv_panel_status.tsv (or any table with Chr, Start, End)")
    ap.add_argument("--catalogue", action="append", default=[], help="NAME:VERSION:TABLE:MANIFEST")
    ap.add_argument("--tolerance", type=int, default=500)
    ap.add_argument("--ins-len-ratio", type=float, default=0.5)
    ap.add_argument("--min-junction-fraction", type=float, default=0.75)
    ap.add_argument("--out-loci", required=True)
    ap.add_argument("--out-junctions", required=True)
    ap.add_argument("--out-check", required=True)
    a = ap.parse_args()

    cats = []
    for spec in a.catalogue:
        name, version, table, manifest = spec.split(":", 3)
        cats.append({"name": name, "version": version, "table": table, "manifest": manifest})

    # verify every table against its manifest before reading it
    checks, bad = [], []
    for c in cats:
        want = manifest_field(c["manifest"], "output_sha256")
        mver = manifest_field(c["manifest"], "version")
        got = sha256(c["table"])
        ok = want == got and (mver is None or mver == c["version"])
        checks.append([c["name"], c["version"], c["table"], got, want or ".", "VERIFIED" if ok else "MISMATCH"])
        if not ok:
            bad.append(f"{c['name']} {c['version']}: {c['table']} sha256 {got} != manifest {want} "
                       f"(manifest version {mver})")
    with open(a.out_check, "w") as fh:
        fh.write("catalogue\tversion\ttable\tsha256\tmanifest_sha256\tstatus\n")
        for r in checks:
            fh.write("\t".join(r) + "\n")
    if bad:
        sys.exit("external catalogue(s) do not match their manifests; refusing them:\n  " + "\n  ".join(bad))

    junctions = cj.read_table(a.junctions)
    by_locus = defaultdict(list)
    for j in junctions:
        for L in j["loci"]:
            by_locus[L].append(j)
    with open(a.loci, newline="") as fh:
        loci = [f"{r['Chr']}:{r['Start']}-{r['End']}" for r in csv.DictReader(fh, delimiter="\t")]
    query = defaultdict(set)
    for L in loci:
        for j in by_locus.get(L, []):
            query[j["chrom1"]].add(j["pos1"])
    query = {k: sorted(v) for k, v in query.items()}

    per_cat = {}
    for c in cats:
        sites = cm.read_catalogue_sites(c["table"], query, a.tolerance)
        per_cat[c["name"]] = cm.SiteIndex(sites) if sites else None
        print(f"{c['name']} {c['version']}: {len(sites)} sites near patient CxSV junctions", file=sys.stderr)

    lcols = ["Locus_ID"]
    for c in cats:
        n = c["name"]
        lcols += [f"{n}_Junctions_Matched", f"{n}_Max_AF", f"{n}_Max_Group_AF", f"{n}_Contains"]
    lcols += ["Catalogues_Containing", "Catalogue_Max_AF"]
    jcols = ["Locus_ID", "chrom1", "pos1", "chrom2", "pos2", "svtype", "svlen", "supp_vec"]
    for c in cats:
        jcols += [f"{c['name']}_sites", f"{c['name']}_max_af", f"{c['name']}_best_site"]

    with open(a.out_loci, "w", newline="") as fl, open(a.out_junctions, "w", newline="") as fj:
        wl = csv.writer(fl, delimiter="\t", lineterminator="\n")
        wj = csv.writer(fj, delimiter="\t", lineterminator="\n")
        wl.writerow(lcols)
        wj.writerow(jcols)
        for L in loci:
            js = by_locus.get(L, [])
            row, containing, contain_af = [L], [], []
            hits_by_cat = {}
            for c in cats:
                idx = per_cat[c["name"]]
                per = [idx.query(j, a.tolerance, a.ins_len_ratio) if idx else [] for j in js]
                hits_by_cat[c["name"]] = per
                k = sum(1 for h in per if h)
                afs = [s["af"] for h in per for s in h if s["af"] is not None]
                gafs = [v for h in per for s in h for v in s["group_af"].values()]
                contains = bool(js) and k / len(js) >= a.min_junction_fraction - 1e-9
                max_af = max(afs) if afs else None
                row += [f"{k}/{len(js)}", fmt(max_af), fmt(max(gafs) if gafs else None), "yes" if contains else "no"]
                if contains:
                    containing.append(c["name"])
                    if max_af is not None:
                        contain_af.append(max_af)
            row += [",".join(containing) or "none", fmt(max(contain_af) if contain_af else None)]
            wl.writerow(row)
            for i, j in enumerate(js):
                r = [L, j["chrom1"], j["pos1"], j["chrom2"], j["pos2"], j["svtype"], j["svlen"], j["supp_vec"]]
                for c in cats:
                    h = hits_by_cat[c["name"]][i]
                    afs = [s["af"] for s in h if s["af"] is not None]
                    best = max(h, key=lambda s: (s["af"] or 0.0)) if h else None
                    r += [len(h), fmt(max(afs) if afs else None),
                          f"{best['svtype']}:{best['chrom1']}:{best['pos1']}-{best['chrom2']}:{best['pos2']}"
                          if best else "."]
                wj.writerow(r)
    print(f"{len(loci)} loci matched against {len(cats)} catalogue(s) "
          f"(tolerance {a.tolerance} bp, min junction fraction {a.min_junction_fraction})", file=sys.stderr)


if __name__ == "__main__":
    main()
