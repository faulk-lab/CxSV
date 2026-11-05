#!/usr/bin/env python3
"""Per-locus privacy and tier of one patient scored against two panel releases
(docs/LIBRARY_QC.md: catalog without the flagged libraries).

    privacy_compare.py --a A/cxsv_scores.tsv --b B/cxsv_scores.tsv --events truth_events.bed ... --out OUT.tsv

Prints counts of loci whose Panel_Status, panel_privacy points or Tier differ, and
every locus overlapping a truth event (BED: chrom start end name).
"""
import argparse
import csv


def t(p):
    return {r["Locus_ID"]: r for r in csv.DictReader(open(p), delimiter="\t")}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--a", required=True)
    ap.add_argument("--b", required=True)
    ap.add_argument("--events", action="append", default=[])
    ap.add_argument("--out", required=True)
    x = ap.parse_args()
    A, B = t(x.a), t(x.b)
    ev = []
    for p in x.events:
        for line in open(p):
            if line.strip() and not line.startswith("#"):
                f = line.split("\t")
                ev.append((f[0], int(f[1]), int(f[2]), f[3].strip()))
    keys = sorted(set(A) | set(B))
    cols = ["Locus_ID", "event", "Cmplx_Class", "Panel_Status_a", "Panel_Status_b", "Panel_Carriers_a",
            "Panel_Carriers_b", "privacy_points_a", "privacy_points_b", "Total_a", "Total_b", "Tier_a", "Tier_b"]
    n = {"loci_a": len(A), "loci_b": len(B), "status_changed": 0, "privacy_changed": 0, "tier_changed": 0}
    with open(x.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for k in keys:
            a, b = A.get(k, {}), B.get(k, {})
            c, s, e = (k.split(":")[0], *map(int, k.split(":")[1].split("-")))
            hit = ",".join(nm for ch, s0, e0, nm in ev if ch == c and s < e0 and e > s0)
            ch = [a.get("Panel_Status") != b.get("Panel_Status"), a.get("panel_privacy_points") != b.get("panel_privacy_points"),
                  a.get("Tier") != b.get("Tier")]
            for name, flag in zip(("status_changed", "privacy_changed", "tier_changed"), ch):
                n[name] += flag
            if hit or any(ch):
                w.writerow([k, hit or ".", (a or b).get("Cmplx_Class"), a.get("Panel_Status"), b.get("Panel_Status"),
                            a.get("Panel_Carriers"), b.get("Panel_Carriers"), a.get("panel_privacy_points"),
                            b.get("panel_privacy_points"), a.get("Total_Score"), b.get("Total_Score"), a.get("Tier"), b.get("Tier")])
    print("\t".join(f"{k}={v}" for k, v in n.items()))


if __name__ == "__main__":
    main()
