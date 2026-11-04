#!/usr/bin/env python3
"""Flag patient CxSV loci that match a published pathogenic rearrangement (roadmap item 11).

Annotation only: no points, no effect on status, score or tier.

Inputs: the patient's junctions (junctions/junctions.tsv), its CxSV loci
(private/cxsv_panel_status.tsv) and the curated literature set
(literature_junctions.tsv + literature_events.tsv, acquisition/literature,
docs/LITERATURE_SET.md). A patient junction matches a literature junction when
the svtype and both contigs agree and each end is within the tolerance of the
literature junction's precision (--tolerance-bp for sequence-resolved junctions,
--tolerance-approximate for array/mate-pair/inverted-repeat positions,
--tolerance-region for ClinVar interval ends), either end order.

Per locus (--out):
  Literature_Match           event     the locus matches >= --min-event-fraction of one published
                                       event's junctions (the same rearrangement, as far as the
                                       published breakpoints resolve it)
                             junction  at least one locus junction matches a published junction
                             region    no junction match, but the locus and a published event's hs1
                                       region (e.g. a ClinVar region-level record) overlap reciprocally
                                       by >= --region-min-overlap (default 0.5) of each one's length
                             none
  Literature_Events          matched events as <event_id>(<junctions matched>/<event junctions>),
                             best first (at most 5)
  Literature_Junctions_Matched   locus junctions matching any literature junction / locus junctions
  Literature_Region_Events   number of events whose region overlaps the locus reciprocally
  Literature_Phenotypes      phenotype labels of the matched events (at most 3)
"""
import argparse
import csv
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import cxsv_junctions as cj  # noqa: E402

COLUMNS = ["Locus_ID", "Literature_Match", "Literature_Events", "Literature_Junctions_Matched",
           "Literature_Region_Events", "Literature_Phenotypes"]


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader([ln for ln in fh if not ln.startswith("#")], delimiter="\t"))


def near(c, p, c0, p0, tol):
    return c == c0 and abs(int(p) - int(p0)) <= tol


def junction_matches(j, lit, tol):
    """Literature junctions matching patient junction j (dict with svtype chrom1 pos1 chrom2 pos2)."""
    out = []
    for L in lit:
        if L["svtype"] != j["svtype"]:
            continue
        t = tol[L["precision"]]
        if (near(j["chrom1"], j["pos1"], L["chrom1"], L["pos1"], t) and near(j["chrom2"], j["pos2"], L["chrom2"],
                                                                            L["pos2"], t)) or \
           (near(j["chrom1"], j["pos1"], L["chrom2"], L["pos2"], t) and near(j["chrom2"], j["pos2"], L["chrom1"],
                                                                            L["pos1"], t)):
            out.append(L)
    return out


def parse_regions(s):
    out = []
    for part in (s or "").split(";"):
        if ":" in part and "-" in part:
            c, r = part.split(":")
            a, b = r.split("-")
            out.append((c, int(a), int(b)))
    return out


def reciprocal(s, e, a, b, f):
    """[s, e] and [a, b] overlap by >= f of each one's length (f = 0: any overlap)."""
    ov = min(e, b) - max(s, a)
    if ov < 0 or (ov == 0 and f > 0):
        return False
    return ov >= f * max(e - s, 1) and ov >= f * max(b - a, 1)


def annotate(loci, junctions, lit_junctions, lit_events, tol, min_fraction, region_min_overlap=0.5):
    by_locus = defaultdict(list)
    for j in junctions:
        for L in j["loci"]:
            by_locus[L].append(j)
    lit = [dict(r) for r in lit_junctions if r.get("lift_status", "lifted") == "lifted"]
    n_event = defaultdict(int)
    for r in lit:
        n_event[r["event_id"]] += 1
    ev = {e["event_id"]: e for e in lit_events}
    regions = {e: parse_regions(x.get("region_hs1")) for e, x in ev.items()}
    rows = []
    for r in loci:
        lid = f"{r['Chr']}:{r['Start']}-{r['End']}"
        s, e = int(r["Start"]), int(r["End"])
        js = by_locus.get(lid, [])
        hit_lit = defaultdict(set)
        k = 0
        for j in js:
            m = junction_matches(j, lit, tol)
            k += bool(m)
            for L in m:
                hit_lit[L["event_id"]].add(L["junction_id"])
        frac = {x: len(v) / n_event[x] for x, v in hit_lit.items()}
        reg = [x for x, ivs in regions.items()
               if any(c == r["Chr"] and reciprocal(s, e, a, b, region_min_overlap) for c, a, b in ivs)]
        if any(f >= min_fraction - 1e-9 for f in frac.values()):
            status = "event"
        elif hit_lit:
            status = "junction"
        elif reg:
            status = "region"
        else:
            status = "none"
        best = sorted(hit_lit, key=lambda x: (-frac[x], x))[:5]
        shown = best or sorted(reg)[:5]
        rows.append({"Locus_ID": lid, "Literature_Match": status,
                     "Literature_Events": ",".join(f"{x}({len(hit_lit[x])}/{n_event[x]})" for x in best) or
                     (",".join(sorted(reg)[:5]) if reg else "."),
                     "Literature_Junctions_Matched": f"{k}/{len(js)}",
                     "Literature_Region_Events": len(reg),
                     "Literature_Phenotypes": "; ".join(dict.fromkeys(ev[x].get("phenotype", ".") for x in shown[:3]
                                                                       if x in ev)) or "."})
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--junctions", required=True)
    ap.add_argument("--loci", required=True)
    ap.add_argument("--literature", required=True, help="literature set directory (literature_junctions.tsv, ...)")
    ap.add_argument("--tolerance-bp", type=int, default=500)
    ap.add_argument("--tolerance-approximate", type=int, default=10000)
    ap.add_argument("--tolerance-region", type=int, default=50000)
    ap.add_argument("--min-event-fraction", type=float, default=0.75)
    ap.add_argument("--region-min-overlap", type=float, default=0.5,
                    help="reciprocal overlap needed for a region-level match (0 = any overlap)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    tol = {"bp": a.tolerance_bp, "approximate": a.tolerance_approximate, "region": a.tolerance_region,
           "mixed": a.tolerance_approximate}
    rows = annotate(read_tsv(a.loci), cj.read_table(a.junctions), read_tsv(f"{a.literature}/literature_junctions.tsv"),
                    read_tsv(f"{a.literature}/literature_events.tsv"), tol, a.min_event_fraction,
                    a.region_min_overlap)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    n = defaultdict(int)
    for r in rows:
        n[r["Literature_Match"]] += 1
    print(f"{len(rows)} loci: " + ", ".join(f"{n[s]} {s}" for s in ("event", "junction", "region", "none")),
          file=sys.stderr)


if __name__ == "__main__":
    main()
