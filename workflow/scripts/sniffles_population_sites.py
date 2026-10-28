#!/usr/bin/env python3
"""Sniffles2 population-mode genotypes as junction sites, compared with the clustered panel sites.

Sniffles2 multi-sample calling (`sniffles --input <tsv of .snf files> --vcf ...`)
re-genotypes every candidate SV in every panel sample from the evidence stored
in the per-sample .snf files. This script turns its VCF into the shared
junction-site schema (source=sniffles2_population) and compares each site's
carriers and allele frequency with the panel's clustered junction sites
(build_panel_sites.py), which come from the three-caller merged calls.

Genotypes: a sample carries a site when its GT has a non-reference allele;
dosage = number of non-reference alleles (capped by ploidy). an sums ploidy over
panel individuals with a called (non-missing) genotype, so missing genotypes do
not count as reference. One individual sequenced twice counts once (its
highest dosage).

Inputs
------
--vcf                 Sniffles2 population VCF
--samples             TSV: snf, sniffles_id, individual_id (rule panel_sniffles_population)
--panel-individuals   panel_individuals.tsv (sex, superpopulation, platforms)
--carriers            panel_junction_carriers.tsv (clustered sites with envelopes)
--fai                 contig order

Outputs
-------
--sites      sniffles_population_sites.tsv (shared schema + af_<group>)
--compare    panel_frequency_comparison.tsv: one row per Sniffles2 site with the
             best-matching clustered site (same svtype and contigs, both ends within
             --tolerance bp of its envelope), carriers and af on both sides
--summary    panel_frequency_comparison_summary.tsv: counts and agreement statistics
"""
import argparse
import csv
import os
import re
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import cxsv_junctions as cj  # noqa: E402
import cxsv_matching as cm   # noqa: E402


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def gt_dosage(gt):
    alleles = [a for a in re.split(r"[/|]", gt) if a not in (".", "")]
    if not alleles:
        return None
    return sum(a != "0" for a in alleles)


def pearson(x, y):
    n = len(x)
    if n < 2:
        return float("nan")
    mx, my = sum(x) / n, sum(y) / n
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    return sxy / (sxx * syy) ** 0.5 if sxx > 0 and syy > 0 else float("nan")


def ranks(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return r


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
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--panel-individuals", required=True)
    ap.add_argument("--carriers", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--tolerance", type=int, default=500)
    ap.add_argument("--ins-len-ratio", type=float, default=0.5)
    ap.add_argument("--min-svlen", type=int, default=50)
    ap.add_argument("--version", default=None)
    ap.add_argument("--fingerprint", default=None,
                    help="calling_fingerprint.json; version = fp-<12 hex> when --version is not given")
    ap.add_argument("--sites", required=True)
    ap.add_argument("--compare", required=True)
    ap.add_argument("--summary", required=True)
    a = ap.parse_args()
    a.version = resolve_version(a.version, a.fingerprint)

    order = cj.contig_order(a.fai)
    smap = {r["sniffles_id"]: r["individual_id"] for r in read_tsv(a.samples)}
    inds = {r["individual_id"]: {"sex": r["sex"], "superpopulation": r["superpopulation"],
                                  "platforms": set(r["platforms"].split(","))}
            for r in read_tsv(a.panel_individuals)}
    groups = sorted({g for i in inds.values() for g in cm.groups_of(i)},
                    key=lambda g: (0 if g.startswith("platform:") else 1, g))
    index = cm.SiteIndex(cm.read_carriers(a.carriers))

    sn_sites = []
    with open(a.vcf) as fh:
        samples = []
        for line in fh:
            if line.startswith("##"):
                if line.startswith("##contig=<ID="):
                    order.setdefault(line[13:].split(",", 1)[0].rstrip(">\n"), len(order))
                continue
            if line.startswith("#"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            f = line.rstrip("\n").split("\t")
            if f[6] not in ("PASS", "."):
                continue
            j = cj.record_to_junction(f[:8] + ["GT"], order)
            if j is None or cm.too_short(j.svtype, j.svlen, a.min_svlen):
                continue
            keys = f[8].split(":")
            gi = keys.index("GT") if "GT" in keys else None
            dose, called = {}, set()
            for sid, col in zip(samples, f[9:]):
                ind = smap.get(sid, sid)
                if gi is None:
                    continue
                d = gt_dosage(col.split(":")[gi])
                if d is None:
                    continue
                called.add(ind)
                if d > 0:
                    dose[ind] = max(dose.get(ind, 0), d)
            if not dose:
                continue
            tot_ac = tot_an = 0
            per = {g: [0, 0] for g in groups}
            for iid in called:
                ind = inds.get(iid, {"sex": "unknown", "superpopulation": "unknown", "platforms": set()})
                p = cm.site_ploidy(j.chrom1, j.chrom2, ind["sex"])
                ac = min(dose.get(iid, 0), p) if p else dose.get(iid, 0)
                p = p or ac
                tot_ac += ac
                tot_an += p
                for g in cm.groups_of(ind):
                    if g in per:
                        per[g][0] += ac
                        per[g][1] += p
            sn_sites.append({"j": j, "id": f[2], "carriers": sorted(dose), "ac": tot_ac, "an": tot_an,
                             "per": per})
    sn_sites.sort(key=lambda s: (order.get(s["j"].chrom1, 1 << 30), s["j"].pos1,
                                 order.get(s["j"].chrom2, 1 << 30), s["j"].pos2, s["j"].svtype))

    with open(a.sites, "w") as fs:
        fs.write("\t".join(list(cm.SITE_COLUMNS) + [cm.group_column(g) for g in groups]) + "\n")
        for s in sn_sites:
            j = s["j"]
            # Sniffles2 population VCFs carry no INFO/AF, so af = ac/an (shared convention)
            row = [j.chrom1, j.pos1, j.chrom2, j.pos2, j.svtype, cm.svlen_field(j.svtype, j.svlen),
                   s["ac"], s["an"], cm.fmt_af(s["ac"], s["an"]), "sniffles2_population", a.version]
            row += [cm.fmt_af(*s["per"][g]) for g in groups]
            fs.write("\t".join(map(str, row)) + "\n")

    n_panel = len(inds)
    matched_clustered = set()
    comp = []
    for s in sn_sites:
        j = s["j"]
        jd = {"svtype": j.svtype, "chrom1": j.chrom1, "pos1": j.pos1, "chrom2": j.chrom2,
              "pos2": j.pos2, "svlen": j.svlen}
        hits = index.query(jd, a.tolerance, a.ins_len_ratio)
        best = None
        if hits:
            best = min(hits, key=lambda h: (abs(h["pos1"] - j.pos1) + abs(h["pos2"] - j.pos2), h["pos1"]))
            matched_clustered.add((best["svtype"], best["chrom1"], best["pos1"], best["chrom2"], best["pos2"],
                                   best["svlen"]))
        sc = set(s["carriers"])
        cc = best["carrier_set"] if best else set()
        comp.append({
            "svtype": j.svtype, "chrom1": j.chrom1, "pos1": j.pos1, "chrom2": j.chrom2, "pos2": j.pos2,
            "svlen": cm.svlen_field(j.svtype, j.svlen), "sniffles_id": s["id"], "sniffles_carriers": len(sc),
            "sniffles_af": cm.fmt_af(s["ac"], s["an"]),
            "clustered_site": (f"{best['chrom1']}:{best['pos1']}-{best['chrom2']}:{best['pos2']}"
                               if best else "."),
            "clustered_carriers": len(cc) if best else 0,
            "clustered_carrier_freq": f"{len(cc) / n_panel:.6f}" if (best and n_panel) else "0.000000",
            "carriers_both": len(sc & cc), "carriers_sniffles_only": len(sc - cc),
            "carriers_clustered_only": len(cc - sc),
        })
    with open(a.compare, "w", newline="") as fh:
        w = csv.DictWriter(fh, list(comp[0].keys()) if comp else ["none"], delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(comp)

    n_clustered = sum(1 for _ in open(a.carriers)) - 1
    m = [r for r in comp if r["clustered_site"] != "."]
    eq = sum(r["sniffles_carriers"] == r["clustered_carriers"] for r in m)
    more = sum(r["sniffles_carriers"] > r["clustered_carriers"] for r in m)
    less = sum(r["sniffles_carriers"] < r["clustered_carriers"] for r in m)
    xs = [r["sniffles_carriers"] for r in m]
    ys = [r["clustered_carriers"] for r in m]
    by_type = defaultdict(lambda: [0, 0])
    for r in comp:
        by_type[r["svtype"]][0] += 1
        by_type[r["svtype"]][1] += r["clustered_site"] != "."
    tot_both = sum(r["carriers_both"] for r in m)
    tot_s = sum(r["sniffles_carriers"] for r in m)
    tot_c = sum(r["clustered_carriers"] for r in m)
    stats = [
        ("sniffles2_population_sites", len(comp)),
        ("clustered_sites", n_clustered),
        ("sniffles_sites_matched_to_clustered", len(m)),
        ("clustered_sites_matched_by_sniffles", len(matched_clustered)),
        ("matched_equal_carriers", eq),
        ("matched_sniffles_more_carriers", more),
        ("matched_sniffles_fewer_carriers", less),
        ("matched_carrier_pearson", f"{pearson(xs, ys):.4f}"),
        ("matched_carrier_spearman", f"{pearson(ranks(xs), ranks(ys)):.4f}"),
        ("matched_carrier_calls_both", tot_both),
        ("matched_carrier_calls_sniffles_total", tot_s),
        ("matched_carrier_calls_clustered_total", tot_c),
        ("matched_mean_abs_carrier_difference", f"{sum(abs(x - y) for x, y in zip(xs, ys)) / len(m):.4f}" if m else "NA"),
    ] + [(f"by_type_{t}_sniffles_sites_matched", f"{v[1]}/{v[0]}") for t, v in sorted(by_type.items())]
    with open(a.summary, "w") as fh:
        fh.write("statistic\tvalue\n")
        for k, v in stats:
            fh.write(f"{k}\t{v}\n")
    print("\n".join(f"{k}\t{v}" for k, v in stats), file=sys.stderr)


if __name__ == "__main__":
    main()
