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
  scoring/clingen_region_dosage.tsv   ClinGen curated regions (roadmap item 9), lifted
                               GRCh38 -> hs1 with UCSC hg38ToHs1.over.chain.gz (the chain
                               the gnomAD-SV catalogue liftover uses): chrom start end
                               (0-based half-open) isca_id region_name cytoband hi_score
                               ts_score last_evaluated hg38_location lift_status span_ratio
                               interval_pieces interval_mapped_fraction
  scoring/clingen_region_dosage.liftover.tsv   every source region with its lift status
                               (lifted or why not), the chain's checksums in the header
  resources_manifest.tsv       file bytes sha256 md5 source version url source_sha256
                               retrieved_utc licence

OMIM is not fetched: it needs a licence key (scoring.yaml omim, off by default).

Downloads are reused: a source file missing from --downloads is taken from the
newest earlier <root>/_downloads/<date>/ that has it (unless --refresh), else it
is fetched into --downloads. SOURCES.tsv in --downloads lists the file actually
used for every source. A source with expected_md5 in sources.tsv is verified.

Region liftover (build_clingen_regions): each region's two end bases are lifted
independently with acquisition/catalogues/liftover.lift_points (the catalogue
code; default -minMatch, so an end base must sit in an aligned chain block),
re-oriented, and kept when both ends land on the source chromosome and the hs1
span is 0.5-2x the GRCh38 span. The whole interval is also lifted with liftOver
-multiple -minMatch=0.1: interval_pieces > 1 means the region is split across
chain blocks in hs1 (reported, not filtered), and a region whose end base does
not lift is kept when >= 95% of its bases lift onto the source chromosome
(lifted_interval). Regions without
coordinates ('tbd') are not lifted.

Usage:
  python acquisition/resources/build_resources.py [--root DIR] [--downloads DIR]
         [--hpo-tag v2026-09-01] [--download] [--refresh] [--liftover-bin liftOver]
Standard library only; the region liftover needs UCSC liftOver (envs/cxsv-catalogues.yaml).
"""
import argparse
import csv
import datetime as dt
import hashlib
import glob
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import urllib.request
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
CATALOGUES = os.path.join(HERE, "..", "catalogues")


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


# ─── ClinGen region list (roadmap item 9) ────────────────────────────────────
_LOC = re.compile(r"^(chr[0-9XY]+):([0-9,]+)-([0-9,]+)$")
REGION_COLUMNS = ["chrom", "start", "end", "isca_id", "region_name", "cytoband", "hi_score", "ts_score",
                  "last_evaluated", "hg38_location", "lift_status", "span_ratio", "interval_pieces",
                  "interval_mapped_fraction"]


def header_date(txt):
    """'28 Apr,2026' (ClinGen file header) -> '2026-04-28', else None."""
    for fmt in ("%d %b,%Y", "%d %B,%Y", "%d %b, %Y"):
        try:
            return dt.datetime.strptime(txt.strip(), fmt).strftime("%Y-%m-%d")
        except ValueError:
            pass
    return None


def parse_clingen_regions(lines):
    """ClinGen_region_curation_list_GRCh38.tsv -> (version, regions). Each region:
    isca_id, region_name, cytoband, hg38_location, chrom/start/end (0-based
    half-open GRCh38; None when the location is 'tbd' or unparseable), hi_score,
    ts_score, last_evaluated. The 'Genomic Location' is 1-based inclusive."""
    version, header, regions = None, None, []
    for line in lines:
        line = line.rstrip("\n")
        if line.startswith("#ISCA ID"):
            header = line[1:].split("\t")
            continue
        if line.startswith("#"):
            version = version or header_date(line[1:])
            continue
        if not line.strip():
            continue
        if not header:
            raise ValueError("ClinGen region list: data before the '#ISCA ID' header line")
        r = dict(zip(header, line.split("\t")))
        loc = (r.get("Genomic Location") or "").strip()
        m = _LOC.match(loc)
        chrom = start = end = None
        if m:
            chrom, start, end = m.group(1), int(m.group(2).replace(",", "")) - 1, int(m.group(3).replace(",", ""))
            if end <= start:
                chrom = start = end = None
        regions.append({"isca_id": r["ISCA ID"], "region_name": r.get("ISCA Region Name", ""),
                        "cytoband": r.get("cytoBand", ""), "hg38_location": loc, "chrom": chrom, "start": start,
                        "end": end, "hi_score": (r.get("Haploinsufficiency Score") or "").strip(),
                        "ts_score": (r.get("Triplosensitivity Score") or "").strip(),
                        "last_evaluated": r.get("Date Last Evaluated", "")})
    if header is None:
        raise ValueError("ClinGen region list: no '#ISCA ID' header line")
    return version, regions


def lift_regions(regions, lift_points, lift_intervals=None, min_ratio=0.5, max_ratio=2.0,
                 min_interval_fraction=0.95):
    """Lift regions GRCh38 -> hs1 by their two end bases.

    lift_points({(chrom, pos1)}) -> {(chrom, pos1): (chrom', pos1')} (1-based
    bases; acquisition/catalogues/liftover.lift_points). lift_intervals (QC,
    optional) {isca_id: (chrom, start, end)} -> {isca_id: [(chrom', s', e'), ...]}.
    Returns (all regions with lift_status and hs1 coordinates where lifted,
    Counter of statuses). Statuses: lifted, lifted_span_changed (10-50%),
    lifted_interval (an end base did not lift or landed on another contig, but
    >= min_interval_fraction of the region's bases lift onto the source
    chromosome; extent of those pieces), unlifted_no_coordinates,
    unlifted_end_unmapped, unlifted_split_contig, unlifted_span_changed (hs1 span
    outside min_ratio-max_ratio x the GRCh38 span)."""
    import junction_sites  # acquisition/catalogues (contig names)
    placed = [r for r in regions if r["chrom"]]
    pts = {(r["chrom"], r["start"] + 1) for r in placed} | {(r["chrom"], r["end"]) for r in placed}
    mapped = lift_points(pts)
    pieces = lift_intervals({r["isca_id"]: (r["chrom"], r["start"], r["end"]) for r in placed}) \
        if lift_intervals else {}
    out, stats = [], Counter()
    for r in regions:
        n = dict(r, src_chrom=r["chrom"], src_start=r["start"], src_end=r["end"], span_ratio="",
                 interval_pieces="", interval_mapped_fraction="")
        n["chrom"] = n["start"] = n["end"] = None
        if not r["chrom"]:
            n["lift_status"] = "unlifted_no_coordinates"
        else:
            a, b = mapped.get((r["chrom"], r["start"] + 1)), mapped.get((r["chrom"], r["end"]))
            ca = junction_sites.normalise_contig(a[0]) if a else None
            cb = junction_sites.normalise_contig(b[0]) if b else None
            if lift_intervals:
                pv = pieces.get(r["isca_id"], [])
                n["interval_pieces"] = len(pv)
                n["interval_mapped_fraction"] = f"{sum(e - s for _, s, e in pv) / (r['end'] - r['start']):.3f}"
            src_c = junction_sites.normalise_contig(r["chrom"])
            if ca is None or cb is None or ca != cb or ca != src_c:
                n["lift_status"] = "unlifted_end_unmapped" if ca is None or cb is None else "unlifted_split_contig"
                # fallback: the whole interval lifts (>= min_interval_fraction of its bases) onto the source chromosome
                pv = [(c, s, e) for c, s, e in pieces.get(r["isca_id"], []) if junction_sites.normalise_contig(c) == src_c]
                size = r["end"] - r["start"]
                if pv and sum(e - s for _, s, e in pv) >= min_interval_fraction * size:
                    lo, hi = min(s for _, s, _ in pv), max(e for _, _, e in pv)
                    ratio = (hi - lo) / size
                    if min_ratio <= ratio <= max_ratio:
                        n["chrom"], n["start"], n["end"] = src_c, lo, hi
                        n["span_ratio"] = f"{ratio:.3f}"
                        n["lift_status"] = "lifted_interval"
            else:
                lo, hi = sorted((a[1], b[1]))
                ratio = (hi - lo + 1) / (r["end"] - r["start"])
                n["span_ratio"] = f"{ratio:.3f}"
                if not (min_ratio <= ratio <= max_ratio):
                    n["lift_status"] = "unlifted_span_changed"
                else:
                    n["chrom"], n["start"], n["end"] = ca, lo - 1, hi
                    n["lift_status"] = "lifted" if abs(ratio - 1) <= 0.1 else "lifted_span_changed"
        stats[n["lift_status"]] += 1
        if lift_intervals and r["chrom"] and n["chrom"] and int(n["interval_pieces"] or 0) > 1:
            stats["lifted_split_in_hs1 (QC: >1 chain piece)"] += 1
        out.append(n)
    return out, stats


def liftover_intervals(liftover_bin, chain, min_match=0.1):
    """lift_intervals for lift_regions: UCSC liftOver -multiple on whole regions."""
    def run(ivs):
        with tempfile.TemporaryDirectory() as tmp:
            bed, out, unm = (os.path.join(tmp, x) for x in ("in.bed", "out.bed", "unmapped.bed"))
            with open(bed, "w") as fh:
                for k, (c, s, e) in sorted(ivs.items()):
                    fh.write(f"{c}\t{s}\t{e}\t{k}\n")
            subprocess.run([liftover_bin, "-multiple", f"-minMatch={min_match}", bed, chain, out, unm],
                           check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            res = {}
            with open(out) as fh:
                for line in fh:
                    c, s, e, k = line.split("\t")[:4]
                    res.setdefault(k.strip(), []).append((c, int(s), int(e)))
        return res
    return run


def write_region_tables(lifted, out, report, chain_info):
    def score_key(r):
        return (r["chrom"], r["start"], r["end"], r["isca_id"])

    def write(fh):
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(REGION_COLUMNS)
        for r in sorted((r for r in lifted if r["chrom"]), key=score_key):
            w.writerow([r[c] for c in REGION_COLUMNS])
    atomic_write(out, write)

    def write_report(fh):
        for k, v in chain_info.items():
            fh.write(f"## {k}: {v}\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["isca_id", "region_name", "hg38_location", "hi_score", "ts_score", "lift_status", "hs1_chrom",
                    "hs1_start", "hs1_end", "span_ratio", "interval_pieces", "interval_mapped_fraction"])
        for r in sorted(lifted, key=lambda r: r["isca_id"]):
            w.writerow([r["isca_id"], r["region_name"], r["hg38_location"], r["hi_score"], r["ts_score"],
                        r["lift_status"], r["chrom"] or "", "" if r["start"] is None else r["start"],
                        "" if r["end"] is None else r["end"], r["span_ratio"], r["interval_pieces"],
                        r["interval_mapped_fraction"]])
    atomic_write(report, write_report)


def build_clingen_regions(src, chain, liftover_bin, out, report, chain_info):
    sys.path.insert(0, CATALOGUES)
    import liftover  # acquisition/catalogues/liftover.py, the catalogue liftover
    with open(src) as fh:
        version, regions = parse_clingen_regions(fh)
    lifted, stats = lift_regions(regions, lambda pts: liftover.lift_points(pts, chain, liftover_bin),
                                 liftover_intervals(liftover_bin, chain))
    write_region_tables(lifted, out, report, chain_info)
    scored = Counter()
    for r in lifted:
        if r["hi_score"] in ("1", "2", "3") or r["ts_score"] in ("1", "2", "3"):
            scored["kept" if r["chrom"] else "lost"] += 1
    return {"version": version or "unknown", "regions": len(regions), "stats": stats, "scored": scored}


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
    ap.add_argument("--refresh", action="store_true",
                    help="never reuse an earlier _downloads/<date>/ copy (fetch into --downloads)")
    ap.add_argument("--liftover-bin", default="liftOver", help="UCSC liftOver (region list liftover)")
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
        if not os.path.exists(path) and not a.refresh:
            earlier = sorted(glob.glob(os.path.join(a.root, "_downloads", "*", s["file"])))
            path = earlier[-1] if earlier else path
        if not os.path.exists(path):
            if not a.download:
                sys.exit(f"missing {path} (run with --download)")
            print(f"downloading {s['url']}", file=sys.stderr)
            fetch(s["url"], path)
        s["path"] = path
        s["md5"], s["sha256"] = digests(path)
        s["retrieved_utc"] = utc_mtime(path)
        if (s.get("expected_md5") or "").strip() and s["md5"] != s["expected_md5"].strip():
            sys.exit(f"{path}: md5 {s['md5']} != expected {s['expected_md5']} (sources.tsv)")
    # source list next to the downloads
    def write_sources(fh):
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["name", "version", "file", "url", "retrieved_utc", "bytes", "md5", "sha256", "licence",
                    "download_dir"])
        for s in sources.values():
            w.writerow([s["name"], s["version"], s["file"], s["url"], s["retrieved_utc"],
                        os.path.getsize(s["path"]), s["md5"], s["sha256"], s["licence"],
                        os.path.dirname(s["path"])])
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
    ch = sources["UCSC_hg38ToHs1_chain"]
    chain_info = {"chain": ch["file"], "chain_url": ch["url"], "chain_path": ch["path"], "chain_md5": ch["md5"],
                  "chain_sha256": ch["sha256"], "liftover": "end bases via acquisition/catalogues/liftover.py "
                  "(liftOver default -minMatch); kept if both ends on the source chromosome and hs1 span 0.5-2x; "
                  "QC: whole interval with -multiple -minMatch=0.1",
                  "source": sources["ClinGen_regions"]["file"], "source_sha256": sources["ClinGen_regions"]["sha256"]}
    info = build_clingen_regions(sources["ClinGen_regions"]["path"], ch["path"], a.liftover_bin,
                                 os.path.join(R, "scoring", "clingen_region_dosage.tsv"),
                                 os.path.join(R, "scoring", "clingen_region_dosage.liftover.tsv"), chain_info)
    st = info["stats"]
    note = (f"{info['regions']} regions; " + ", ".join(f"{k} {v}" for k, v in sorted(st.items()))
            + f"; HI/TS 1-3 regions kept {info['scored']['kept']}, lost {info['scored']['lost']}"
            + f"; lifted with {ch['file']} md5 {ch['md5']} sha256 {ch['sha256']}")
    print("ClinGen regions: " + note, file=sys.stderr)
    built.append(("scoring/clingen_region_dosage.tsv", "ClinGen_regions", info["version"], note))
    built.append(("scoring/clingen_region_dosage.liftover.tsv", "ClinGen_regions", info["version"],
                  "per-region liftover status (all source regions)"))
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
