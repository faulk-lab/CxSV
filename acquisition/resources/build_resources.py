#!/usr/bin/env python3
"""Fetch and build the annotation + scoring resources (roadmap item 5).

Every file used by the scorer, and the rebuilt repeat annotation, comes from a
free source listed in sources.tsv, is downloaded once into
<root>/_downloads/<YYYYMMDD>/ and is recorded with its URL, retrieval time and
checksums in <root>/resources_manifest.tsv (one row per built file).

Built under --root (default /mnt/sargasso1/CxSV/resources):

  repeats.final.bed            per-element RepeatMasker records (chrom start end class)
                               from the T2T consortium's CHM13v2 RepeatMasker
                               4.1.2p1.2022Apr14 annotation, replacing UCSC "joined"
                               spans; the previous file is kept as repeats.final.bed.bak
  scoring/clingen_dosage.tsv   gene_symbol gene_id hi_score ts_score last_evaluated
  scoring/gnomad_constraint.tsv  gene_symbol transcript loeuf pli transcript_choice
                               (MANE Select, else Ensembl canonical; lowest LOEUF if several)
  scoring/hpo_genes_to_phenotype.txt, scoring/hp.obo   HPO release (optional HPO match)
  scoring/par.bed              CHM13v2 pseudo-autosomal regions (chrX/chrY)
  resources_manifest.tsv       file bytes sha256 md5 source version url source_sha256
                               retrieved_utc licence

OMIM is not fetched: it needs a licence key (scoring.yaml omim, off by default).

Usage:
  python acquisition/resources/build_resources.py [--root DIR] [--downloads DIR]
         [--hpo-tag v2026-09-01] [--download]
Standard library only.
"""
import argparse
import csv
import datetime as dt
import hashlib
import json
import os
import shutil
import sys
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))


def digests(path):
    md5, sha = hashlib.md5(), hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            md5.update(block)
            sha.update(block)
    return md5.hexdigest(), sha.hexdigest()


def utc_mtime(path):
    return dt.datetime.fromtimestamp(os.path.getmtime(path), dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def atomic_write(path, write):
    tmp = path + ".tmp"
    with open(tmp, "w", newline="") as fh:
        write(fh)
    os.replace(tmp, path)


def latest_hpo_tag():
    with urllib.request.urlopen("https://api.github.com/repos/obophenotype/human-phenotype-ontology/releases/latest",
                                timeout=60) as r:
        return json.load(r)["tag_name"]


def fetch(url, dest):
    tmp = dest + ".part"
    with urllib.request.urlopen(url, timeout=600) as r, open(tmp, "wb") as fh:
        shutil.copyfileobj(r, fh, 1 << 20)
    os.replace(tmp, dest)


# ─── builders ────────────────────────────────────────────────────────────────
def build_repeats(src, out):
    rows = []
    with open(src) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 7 or line.startswith("#"):
                continue
            rows.append((f[0], int(f[1]), int(f[2]), f[6]))
    rows.sort()
    atomic_write(out, lambda fh: fh.writelines(f"{c}\t{s}\t{e}\t{k}\n" for c, s, e, k in rows))
    return {"records": len(rows)}


def build_clingen(src, out):
    version, header, rows = None, None, []
    with open(src, newline="") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#Gene Symbol"):
                header = line[1:].split("\t")
                continue
            if line.startswith("#"):
                txt = line[1:].strip()
                for fmt in ("%d %b,%Y", "%d %B,%Y", "%d %b, %Y"):
                    try:
                        version = dt.datetime.strptime(txt, fmt).strftime("%Y-%m-%d")
                    except ValueError:
                        pass
                continue
            if header and line:
                rows.append(dict(zip(header, line.split("\t"))))
    if not header:
        sys.exit(f"{src}: no '#Gene Symbol' header line")

    def write(fh):
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["gene_symbol", "gene_id", "hi_score", "ts_score", "last_evaluated"])
        for r in sorted(rows, key=lambda r: r["Gene Symbol"]):
            w.writerow([r["Gene Symbol"], r.get("Gene ID", ""), r.get("Haploinsufficiency Score", "").strip(),
                        r.get("Triplosensitivity Score", "").strip(), r.get("Date Last Evaluated", "")])
    atomic_write(out, write)
    return {"genes": len(rows), "version": version or "unknown"}


def build_constraint(src, out):
    best = {}
    with open(src, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            choice = "mane_select" if r.get("mane_select") == "true" else (
                "canonical" if r.get("canonical") == "true" else None)
            loeuf = r.get("lof.oe_ci.upper", "NA")
            if choice is None or loeuf in ("NA", ""):
                continue
            rank = 0 if choice == "mane_select" else 1
            key = (rank, float(loeuf))
            g = r["gene"]
            if g not in best or key < best[g][0]:
                best[g] = (key, r["transcript"], loeuf, r.get("lof.pLI", "NA"), choice)

    def write(fh):
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["gene_symbol", "transcript", "loeuf", "pli", "transcript_choice"])
        for g in sorted(best):
            _, tx, loeuf, pli, choice = best[g]
            w.writerow([g, tx, loeuf, pli, choice])
    atomic_write(out, write)
    return {"genes": len(best)}


def obo_version(path):
    with open(path) as fh:
        for line in fh:
            if line.startswith("data-version:"):
                return line.split(":", 1)[1].strip()
            if line.startswith("[Term]"):
                break
    return "unknown"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--root", default="/mnt/sargasso1/CxSV/resources")
    ap.add_argument("--downloads", default=None, help="download dir (default <root>/_downloads/<today>)")
    ap.add_argument("--hpo-tag", default=None, help="HPO release tag (default: latest release)")
    ap.add_argument("--download", action="store_true", help="fetch files missing from --downloads")
    a = ap.parse_args()
    dl = a.downloads or os.path.join(a.root, "_downloads", dt.date.today().strftime("%Y%m%d"))
    os.makedirs(dl, exist_ok=True)
    os.makedirs(os.path.join(a.root, "scoring"), exist_ok=True)
    with open(os.path.join(HERE, "sources.tsv"), newline="") as fh:
        sources = {r["name"]: r for r in csv.DictReader(fh, delimiter="\t")}
    hpo_tag = a.hpo_tag
    for s in sources.values():
        if "{HPO_TAG}" in s["url"]:
            hpo_tag = hpo_tag or latest_hpo_tag()
            s["url"] = s["url"].replace("{HPO_TAG}", hpo_tag)
        path = os.path.join(dl, s["file"])
        if not os.path.exists(path):
            if not a.download:
                sys.exit(f"missing {path} (run with --download)")
            print(f"downloading {s['url']}", file=sys.stderr)
            fetch(s["url"], path)
        s["path"] = path
        s["md5"], s["sha256"] = digests(path)
        s["retrieved_utc"] = utc_mtime(path)
    # source list next to the downloads
    def write_sources(fh):
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["name", "version", "file", "url", "retrieved_utc", "bytes", "md5", "sha256", "licence"])
        for s in sources.values():
            w.writerow([s["name"], s["version"], s["file"], s["url"], s["retrieved_utc"],
                        os.path.getsize(s["path"]), s["md5"], s["sha256"], s["licence"]])
    atomic_write(os.path.join(dl, "SOURCES.tsv"), write_sources)

    R = a.root
    built = []
    rep = os.path.join(R, "repeats.final.bed")
    if os.path.exists(rep) and not os.path.exists(rep + ".bak"):
        shutil.copy2(rep, rep + ".bak")
        print(f"kept the previous repeat annotation as {rep}.bak", file=sys.stderr)
    info = build_repeats(sources["RepeatMasker_CHM13v2"]["path"], rep)
    built.append(("repeats.final.bed", "RepeatMasker_CHM13v2", "4.1.2p1.2022Apr14",
                  f"per-element records, column 4 = repeat class; {info['records']} records"))
    info = build_clingen(sources["ClinGen_dosage"]["path"], os.path.join(R, "scoring", "clingen_dosage.tsv"))
    built.append(("scoring/clingen_dosage.tsv", "ClinGen_dosage", info["version"], f"{info['genes']} genes"))
    info = build_constraint(sources["gnomAD_constraint"]["path"], os.path.join(R, "scoring", "gnomad_constraint.tsv"))
    built.append(("scoring/gnomad_constraint.tsv", "gnomAD_constraint", "v4.1", f"{info['genes']} genes"))
    for name, out in (("HPO_genes_to_phenotype", "scoring/hpo_genes_to_phenotype.txt"),
                      ("HPO_ontology", "scoring/hp.obo"), ("PAR_CHM13v2", "scoring/par.bed")):
        tmp = os.path.join(R, out) + ".tmp"
        shutil.copyfile(sources[name]["path"], tmp)
        os.replace(tmp, os.path.join(R, out))
        ver = hpo_tag if name.startswith("HPO") else sources[name]["version"]
        note = f"data-version {obo_version(os.path.join(R, out))}" if name == "HPO_ontology" else "copied as is"
        built.append((out, name, ver, note))

    def write_manifest(fh):
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["file", "bytes", "sha256", "md5", "source", "source_version", "source_url",
                    "source_file", "source_sha256", "retrieved_utc", "licence", "note"])
        for out, name, ver, note in built:
            s = sources[name]
            md5, sha = digests(os.path.join(R, out))
            w.writerow([out, os.path.getsize(os.path.join(R, out)), sha, md5, name, ver, s["url"], s["file"],
                        s["sha256"], s["retrieved_utc"], s["licence"], note])
        for out in ("genes.final.bed", "exons.filtered.bed", "segdup.final.bed", "repeats.final.bed.bak"):
            p = os.path.join(R, out)
            if os.path.exists(p):
                md5, sha = digests(p)
                w.writerow([out, os.path.getsize(p), sha, md5, "pre-existing", ".", ".", ".", ".", utc_mtime(p), ".",
                            "built before roadmap item 5 (UCSC hs1 tracks); provenance not recorded"
                            + ("; superseded UCSC 'joined' RepeatMasker spans" if out.endswith(".bak") else "")])
    atomic_write(os.path.join(R, "resources_manifest.tsv"), write_manifest)
    print(f"wrote {len(built)} resources and {os.path.join(R, 'resources_manifest.tsv')}", file=sys.stderr)


if __name__ == "__main__":
    main()
