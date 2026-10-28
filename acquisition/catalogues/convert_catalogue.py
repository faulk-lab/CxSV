#!/usr/bin/env python3
"""Convert one downloaded SV catalogue to the shared junction-site schema.

    python3 convert_catalogue.py --catalogue HGSVC3 \
        --root /mnt/sargasso1/CxSV/datasets/external_catalogues \
        --fai /mnt/sargasso1/CxSV/refs/hs1.fa.fai [--threads 8]

Reads <root>/<catalogue>/<version>/ (as filled by download_catalogues.py) and
writes, next to the source files:

    <catalogue>.<version>.junction_sites.tsv.gz      bgzipped, tabix-indexed
    <catalogue>.<version>.junction_sites.tsv.gz.tbi  (-s1 -b2 -e2, header skipped)
    <catalogue>.<version>.manifest.tsv               field / value
    <catalogue>.<version>.liftover_stats.tsv         (GRCh38 sources only)

Needs bcftools, bgzip and tabix on PATH (conda env envs/cxsv-catalogues.yaml),
and UCSC liftOver for gnomAD-SV. The per-record logic is in adapters.py; the
schema rules are in junction_sites.py; see README.md for every choice made.

Rules applied to every catalogue:
* SVs shorter than 50 bp (net length) are dropped (``dropped_below_min_len``).
* Sites with AC = 0 (no carrier left) are dropped.
* Alleles split from one multi-allelic record that give the same junction
  site (same ends, type and svlen) are merged, summing AC (HPRC graph VCF).
* Exact duplicate rows (same ends, type, svlen, ac, an) are written once.
"""
import argparse
import csv
import datetime as dt
import hashlib
import os
import subprocess
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import adapters  # noqa: E402
import junction_sites as js  # noqa: E402
import liftover  # noqa: E402

# Catalogue metadata. Licence and cohort facts are from the official sources
# cited in docs/EXTERNAL_CATALOGUES.md.
CATALOGUES = {
    "CoLoRSdb": {
        "version": "v1.2.0", "native_build": "CHM13 (T2T-CHM13v2.0, maskedY + rCRS)",
        "licence": "CC-BY-4.0", "platform": "PacBio HiFi (pbsv + Jasmine)",
        "af_groups": [],
    },
    "HGSVC3": {
        "version": "HGSVC2024v1.0", "native_build": "T2T-CHM13v2.0 (hs1)",
        "licence": "IGSR/HGSVC3 data reuse statement (Fort Lauderdale; open, no access control)",
        "platform": "PacBio HiFi + ONT ultra-long phased Verkko assemblies (PAV)",
        "af_groups": [],
    },
    "HPRC": {
        "version": "release2-mc-v2.1", "native_build": "CHM13v2.0 (Minigraph-Cactus reference path)",
        "licence": "Public domain (HPRC Data Use: not patented nor subject to copyright)",
        "platform": "PacBio HiFi + ONT ultra-long phased assemblies (Minigraph-Cactus graph, vcfwave)",
        "af_groups": [],
    },
    "gnomAD-SV": {
        "version": "v4.1", "native_build": "GRCh38 (lifted to hs1)",
        "licence": "CC0-1.0", "platform": "Illumina short-read WGS (GATK-SV)",
        "af_groups": [p.upper() for p in adapters.GNOMAD_POPS],
    },
}


def run_query(cmd_pipeline):
    """Yield tab-split lines from a shell pipeline (bcftools ...)."""
    p = subprocess.Popen(["bash", "-o", "pipefail", "-c", cmd_pipeline],
                         stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    for line in p.stdout:
        yield line.rstrip("\n").split("\t")
    if p.wait() != 0:
        raise RuntimeError(f"pipeline failed ({p.returncode}): {cmd_pipeline}")


def nz(v):
    return None if v in (".", "") else v


# ----------------------------------------------------------- per catalogue
def convert_colorsdb(d, threads, stats):
    vcf = os.path.join(d, "CoLoRSdb.CHM13.v1.2.0.pbsv.jasmine.vcf.gz")
    q = (f"bcftools norm --threads {threads} -m -any {vcf} -Ou | bcftools query -f "
         "'%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER\\t%INFO/SVTYPE\\t%INFO/END\\t"
         "%INFO/SVLEN\\t%INFO/AC\\t%INFO/AN\\t%INFO/AF\\t%INFO/NS\\n'")
    sites, max_ns = [], 0
    for f in run_query(q):
        stats["n_records_read"] += 1
        stats[f"filter_{f[4]}"] += 1
        stats[f"source_svtype_{f[5]}"] += 1
        if nz(f[11]):
            max_ns = max(max_ns, int(f[11]))
        s, why = adapters.vcf_sv_sites(f[0], f[1], f[2], f[3], nz(f[5]), nz(f[6]),
                                       nz(f[7]), ac=nz(f[8]), an=nz(f[9]), af=nz(f[10]))
        if why:
            stats[f"dropped_{why}"] += 1
        sites += s
    stats["max_NS_individuals"] = max_ns
    return sites


def convert_hgsvc3(d, threads, stats):
    sites = []
    for name in ("variants_T2T-CHM13_sv_insdel_sym_HGSVC2024v1.0.vcf.gz",
                 "variants_T2T-CHM13_sv_inv_sym_HGSVC2024v1.0.vcf.gz"):
        vcf = os.path.join(d, name)
        q = (f"bcftools +fill-tags {vcf} -Ou -- -t AC,AN | bcftools query -f "
             "'%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER\\t%INFO/SVTYPE\\t%INFO/SVLEN\\t"
             "%INFO/AC\\t%INFO/AN\\n'")
        for f in run_query(q):
            stats["n_records_read"] += 1
            stats[f"filter_{f[4]}"] += 1
            stats[f"source_svtype_{f[5]}"] += 1
            s, why = adapters.vcf_sv_sites(f[0], f[1], f[2], f[3], nz(f[5]), None,
                                           nz(f[6]), ac=nz(f[7]), an=nz(f[8]))
            if why:
                stats[f"dropped_{why}"] += 1
            sites += s
    bed = os.path.join(d, "complex_variants_T2T-CHM13_HGSVC2024v1.0.bed.gz")
    import gzip
    with gzip.open(bed, "rt") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            stats["n_records_read"] += 1
            stats["source_svtype_CPX"] += 1
            s, why, tiled = adapters.hgsvc3_cpx_sites(
                r["#CHROM"], r["POS"], r["END"], r["SVLEN"], r["STRUCT_REF"],
                nz(r["AC"]), nz(r["AN"]), nz(r["AF"]))
            if why:
                stats[f"dropped_{why}"] += 1
            stats["cpx_tiled_decomposed" if tiled else "cpx_outer_span_only"] += 1
            sites += s
    return sites


def convert_hprc(d, threads, stats):
    vcf = os.path.join(d, "hprc-v2.1-mc-chm13.wave.vcf.gz")
    # Drop the GRCh38 haplotype (a reference path, not an individual); bcftools
    # view recomputes INFO/AC and INFO/AN for the remaining samples.
    q = (f"bcftools view --threads {threads} -s ^GRCh38 {vcf} -Ou | "
         f"bcftools norm --threads {threads} -m -any -Ou | "
         "bcftools view -G -i 'abs(strlen(REF)-strlen(ALT))>=50 || INFO/INV=1' -Ou | "
         "bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/AC\\t%INFO/AN\\t%INFO/INV\\n'")
    sites, group, group_key = [], [], None
    for f in run_query(q):
        stats["n_records_read"] += 1
        inv = f[6] == "1"
        if inv:
            stats["source_INV_flag"] += 1
        # alleles split from one multi-allelic record share CHROM, POS and REF
        k = (f[0], f[1], f[2])
        if k != group_key:
            sites += merge_record_alleles(group, stats)
            group, group_key = [], k
        s, why = adapters.vcf_sv_sites(f[0], f[1], f[2], f[3], None, None, None,
                                       inv_flag=inv, ac=nz(f[4]), an=nz(f[5]))
        if why:
            stats[f"dropped_{why}"] += 1
        group += s
    sites += merge_record_alleles(group, stats)
    return sites


def merge_record_alleles(group, stats):
    """Merge alleles of ONE source record that give the same junction site.

    Different ALT sequences of equal length at the same anchor (common at
    tandem repeats in a graph VCF) are one junction site; a haplotype carries
    at most one ALT of a record, so their AC values add up (AN is shared).
    """
    merged = {}
    for s in group:
        k = (s["chrom1"], s["pos1"], s["chrom2"], s["pos2"], s["svtype"], s["svlen"])
        if k in merged:
            m = merged[k]
            m["ac"] = (m["ac"] or 0) + (s["ac"] or 0)
            stats["alleles_merged_same_record"] += 1
        else:
            merged[k] = dict(s)
    return list(merged.values())


def convert_gnomad(d, threads, stats, liftover_bin, workdir):
    vcf = os.path.join(d, "gnomad.v4.1.sv.sites.vcf.gz")
    pops = adapters.GNOMAD_POPS
    fmt = ("%CHROM\\t%POS\\t%ALT\\t%FILTER\\t%INFO/SVTYPE\\t%INFO/END\\t%INFO/SVLEN\\t"
           "%INFO/CHR2\\t%INFO/END2\\t%INFO/CPX_INTERVALS\\t%INFO/AC\\t%INFO/AN\\t%INFO/AF"
           + "".join(f"\\t%INFO/AF_{p}" for p in pops) + "\\n")
    q = f"bcftools query -f '{fmt}' {vcf}"
    sites = []
    for f in run_query(q):
        stats["n_records_read"] += 1
        stats[f"filter_{f[3]}"] += 1
        stats[f"source_svtype_{f[4]}"] += 1
        pop_af = {p.upper(): float(v) for p, v in zip(pops, f[13:]) if nz(v)}
        s, why = adapters.gnomad_sites(f[0], f[1], f[2], f[3], f[4], nz(f[5]), nz(f[6]),
                                       nz(f[7]), nz(f[8]), nz(f[9]), nz(f[10]),
                                       nz(f[11]), nz(f[12]), pop_af)
        if why:
            stats[f"dropped_{why}"] += 1
        sites += s
    chain = os.path.join(d, "hg38ToHs1.over.chain.gz")
    lifted, lstats = liftover.lift_sites(sites, chain, liftover_bin, workdir)
    return lifted, lstats


# ----------------------------------------------------------------- output
def file_hashes(path):
    md5, sha = hashlib.md5(), hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            md5.update(chunk)
            sha.update(chunk)
    return md5.hexdigest(), sha.hexdigest()


def dedupe(sites, order):
    seen, out = set(), []
    for s in sites:
        s = js.orient(s, order)
        k = (s["chrom1"], s["pos1"], s["chrom2"], s["pos2"], s["svtype"],
             s["svlen"], s["ac"], s["an"])
        if k in seen:
            continue
        seen.add(k)
        out.append(s)
    return out


def read_sources(d):
    p = os.path.join(d, "SOURCES.tsv")
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--catalogue", required=True, choices=sorted(CATALOGUES))
    ap.add_argument("--root", required=True)
    ap.add_argument("--fai", required=True, help="hs1.fa.fai (defines contig order)")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--liftover-bin", default="liftOver")
    ap.add_argument("--code-version", default="unknown",
                    help="git commit of this code, recorded in the manifest")
    args = ap.parse_args()

    meta = CATALOGUES[args.catalogue]
    name, ver = args.catalogue, meta["version"]
    d = os.path.join(args.root, name, ver)
    order = js.ContigOrder.from_fai(args.fai)
    stats, lstats = Counter(), None

    if name == "CoLoRSdb":
        sites = convert_colorsdb(d, args.threads, stats)
    elif name == "HGSVC3":
        sites = convert_hgsvc3(d, args.threads, stats)
    elif name == "HPRC":
        sites = convert_hprc(d, args.threads, stats)
    else:
        sites, lstats = convert_gnomad(d, args.threads, stats, args.liftover_bin, d)

    n_before = len(sites)
    sites = dedupe(sites, order)
    stats["duplicate_rows_collapsed"] = n_before - len(sites)

    prefix = os.path.join(d, f"{name}.{ver}")
    out_gz = prefix + ".junction_sites.tsv.gz"
    tmp_tsv = prefix + ".junction_sites.tsv.tmp"
    with open(tmp_tsv, "w", newline="") as fh:
        counts = js.write_sites(fh, sites, name, ver, meta["af_groups"], order)
    subprocess.run(f"bgzip -@ {args.threads} -c {tmp_tsv} > {out_gz}.tmp", shell=True, check=True)
    os.remove(tmp_tsv)
    os.replace(out_gz + ".tmp", out_gz)
    subprocess.run(["tabix", "-f", "-S", "1", "-s", "1", "-b", "2", "-e", "2", out_gz], check=True)
    o_md5, o_sha = file_hashes(out_gz)

    rows = [("schema_version", js.SCHEMA_VERSION), ("source", name), ("version", ver),
            ("native_build", meta["native_build"]), ("output_build", "hs1 (T2T-CHM13v2.0, UCSC names)"),
            ("licence", meta["licence"]), ("platform", meta["platform"]),
            ("af_groups", ",".join(meta["af_groups"]) or "none"),
            ("generated_utc", dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")),
            ("code_version", args.code_version),
            ("output_file", os.path.basename(out_gz)), ("output_md5", o_md5),
            ("output_sha256", o_sha)]
    for src in read_sources(d):
        rows += [(f"input:{src['file']}:md5", src["computed_md5"]),
                 (f"input:{src['file']}:sha256", src["computed_sha256"]),
                 (f"input:{src['file']}:verified", src["verified"])]
    rows += [(k, v) for k, v in sorted(stats.items())]
    for t in js.SVTYPES:
        rows.append((f"n_sites_{t}", counts.get(t, 0)))
    rows.append(("n_sites_total", sum(counts.values())))

    if lstats is not None:
        lp = prefix + ".liftover_stats.tsv"
        cols = ["svtype", "n_in", "n_out", "lost_unmapped_end", "lost_split_contig",
                "loss_pct", "span_changed_gt50pct"]
        tot = Counter()
        with open(lp + ".tmp", "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for t in list(js.SVTYPES) + ["ALL"]:
                c = tot if t == "ALL" else lstats.get(t, Counter())
                if t != "ALL":
                    tot.update(c)
                if c["n_in"] == 0:
                    continue
                loss = 100.0 * (c["n_in"] - c["n_out"]) / c["n_in"]
                fh.write("\t".join(str(x) for x in [t, c["n_in"], c["n_out"],
                         c["lost_unmapped_end"], c["lost_split_contig"], f"{loss:.3f}",
                         c["span_changed_gt50pct"]]) + "\n")
        os.replace(lp + ".tmp", lp)
        rows += [("liftover_chain", "hg38ToHs1.over.chain.gz (UCSC)"),
                 ("liftover_tool", "UCSC liftOver, default -minMatch=0.95, 1-bp point per breakpoint"),
                 ("liftover_sites_in", tot["n_in"]), ("liftover_sites_out", tot["n_out"]),
                 ("liftover_loss_pct", f"{100.0 * (tot['n_in'] - tot['n_out']) / max(tot['n_in'], 1):.3f}")]

    mp = prefix + ".manifest.tsv"
    with open(mp + ".tmp", "w") as fh:
        fh.write("field\tvalue\n")
        for k, v in rows:
            fh.write(f"{k}\t{v}\n")
    os.replace(mp + ".tmp", mp)
    print(f"{name} {ver}: {sum(counts.values())} sites {dict(counts)} -> {out_gz}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
