#!/usr/bin/env python3
"""Cluster every panel sample's junctions into panel junction sites with frequencies.

Inputs
------
--manifest     TSV: cohort, sample_id, individual_id, platform, coverage_x, junctions
               (one row per panel sample; junctions = its junctions/junctions.tsv)
--individuals  optional registry individuals.tsv (individual_id, sex, superpopulation, ...);
               without it, sex is "unknown" (chrX counted diploid) and no
               superpopulation groups are written
--fai          reference .fai: contig order for the (chrom1,pos1) <= (chrom2,pos2) convention
--version      written to every row's `version` column (panel release id, or
               the calling fingerprint's short hash)

Outputs
-------
--sites        panel_junction_sites.tsv in the shared junction-site schema
               (docs/JUNCTIONS.md): chrom1 pos1 chrom2 pos2 svtype svlen ac an af
               source version, then af_<platform>... af_<superpopulation>...
--carriers     panel_junction_carriers.tsv: the same sites with their coordinate
               envelope and carrier individuals (what the patient matcher reads)
--panel-individuals  panel_individuals.tsv: individual_id, sex, superpopulation,
               platforms, cohorts, samples (denominators and carrier breakdowns)

Clustering is type-aware and leader-based within --tolerance bp at each end
(workflow/lib/cxsv_matching.py). Frequencies: ac = sum of carrier dosages
(het 1, hom 2, capped by ploidy), an = sum of ploidy over all panel
individuals (chrX/chrY hemizygous in males; PAR not modelled), af = ac/an.
Conventions shared with the catalogue track: svlen unsigned, BND svlen '.',
non-BND junctions shorter than --min-svlen (sv_params.min_sv_len, 50 bp) dropped.
"""
import argparse
import csv
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import cxsv_junctions as cj  # noqa: E402
import cxsv_matching as cm   # noqa: E402


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_individuals(manifest, registry):
    reg = {}
    if registry:
        for r in read_tsv(registry):
            reg[r["individual_id"]] = r
    inds = {}
    for s in manifest:
        i = inds.setdefault(s["individual_id"], {"platforms": set(), "cohorts": set(), "samples": set()})
        i["platforms"].add(s["platform"])
        i["cohorts"].add(s["cohort"])
        i["samples"].add(s["sample_id"])
    for iid, i in inds.items():
        r = reg.get(iid, {})
        sex = (r.get("sex") or "unknown").lower()
        i["sex"] = sex if sex in ("male", "female") else "unknown"
        i["superpopulation"] = r.get("superpopulation") or "unknown"
    return inds


def resolve_version(version, fingerprint):
    """--version wins; else "fp-" + the first 12 hex digits of the fingerprint's sha256."""
    if version:
        return version
    if fingerprint:
        import json
        with open(fingerprint) as fh:
            return "fp-" + json.load(fh)["sha256"][:12]
    raise SystemExit("give --version or --fingerprint")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--individuals", default=None)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--tolerance", type=int, required=True)
    ap.add_argument("--ins-len-ratio", type=float, default=0.5)
    ap.add_argument("--min-svlen", type=int, default=50,
                    help="drop non-BND junctions shorter than this (sv_params.min_sv_len)")
    ap.add_argument("--source", default="cxsv_panel")
    ap.add_argument("--version", default=None)
    ap.add_argument("--fingerprint", default=None,
                    help="calling_fingerprint.json; version = fp-<12 hex> when --version is not given")
    ap.add_argument("--sites", required=True)
    ap.add_argument("--carriers", required=True)
    ap.add_argument("--panel-individuals", required=True)
    a = ap.parse_args()
    a.version = resolve_version(a.version, a.fingerprint)

    manifest = read_tsv(a.manifest)
    order = cj.contig_order(a.fai)
    inds = load_individuals(manifest, a.individuals if a.individuals not in (None, "", "None") else None)

    obs, n_short = [], 0
    for s in manifest:
        for r in cj.read_table(s["junctions"]):
            if cm.too_short(r["svtype"], r["svlen"], a.min_svlen):
                n_short += 1
                continue
            obs.append({"svtype": r["svtype"], "chrom1": r["chrom1"], "pos1": r["pos1"],
                        "chrom2": r["chrom2"], "pos2": r["pos2"], "svlen": r["svlen"],
                        "dosage": r["dosage"], "individual": s["individual_id"],
                        "sample": f"{s['cohort']}/{s['sample_id']}"})
    sites = cm.cluster_sites(obs, a.tolerance, a.ins_len_ratio, order)

    groups = sorted({g for i in inds.values() for g in cm.groups_of(i)},
                    key=lambda g: (0 if g.startswith("platform:") else 1, g))
    gcols = [cm.group_column(g) for g in groups]
    with open(a.sites, "w") as fs, open(a.carriers, "w") as fc:
        fs.write("\t".join(list(cm.SITE_COLUMNS) + gcols) + "\n")
        fc.write("\t".join(cm.CARRIER_COLUMNS) + "\n")
        for s in sites:
            k = cm.site_counts(s, inds, groups)
            base = [s["chrom1"], s["pos1"], s["chrom2"], s["pos2"], s["svtype"],
                    cm.svlen_field(s["svtype"], s["svlen"])]
            row = base + [k["ac"], k["an"], cm.fmt_af(k["ac"], k["an"]), a.source, a.version]
            row += [cm.fmt_af(k["groups"][g]["ac"], k["groups"][g]["an"]) for g in groups]
            fs.write("\t".join(map(str, row)) + "\n")
            crow = base + [s["min_pos1"], s["max_pos1"], s["min_pos2"], s["max_pos2"], s["min_svlen"],
                           s["max_svlen"], s["n_observations"], len(k["carriers"]), ",".join(k["carriers"])]
            fc.write("\t".join(map(str, crow)) + "\n")

    with open(a.panel_individuals, "w") as fh:
        fh.write("individual_id\tsex\tsuperpopulation\tplatforms\tcohorts\tsamples\n")
        for iid in sorted(inds):
            i = inds[iid]
            fh.write(f"{iid}\t{i['sex']}\t{i['superpopulation']}\t{','.join(sorted(i['platforms']))}\t"
                     f"{','.join(sorted(i['cohorts']))}\t{','.join(sorted(i['samples']))}\n")

    by_n = defaultdict(int)
    for s in sites:
        by_n[min(len({m['individual'] for m in s['members']}), 5)] += 1
    print(f"panel sites: {n_short} junctions < {a.min_svlen} bp dropped; "
          f"{len(obs)} junction observations from {len(manifest)} samples / "
          f"{len(inds)} individuals -> {len(sites)} sites (tolerance {a.tolerance} bp); carriers "
          + ", ".join(f"{'>=' if n == 5 else ''}{n}: {c}" for n, c in sorted(by_n.items())), file=sys.stderr)


if __name__ == "__main__":
    main()
