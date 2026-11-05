#!/usr/bin/env python3
"""Catalog impact of excluding library-QC-flagged genomes (docs/LIBRARY_QC.md).

    catalog_impact.py --old OLD_RELEASE --new NEW_RELEASE --flagged flagged.tsv \
        --truth LABEL=junctions.tsv ... --out-prefix OUT [--tolerance 500]

OLD/NEW_RELEASE: release directories (panel_junction_carriers.tsv, panel_cxsv_loci.tsv,
panel_cxsv_calls.tsv). Junction sites are compared by chrom1, pos1, chrom2, pos2, svtype;
an OLD site absent from NEW is "shifted" when a NEW site of the same type lies within
--tolerance at both ends, "removed" otherwise. Truth tables (E01-E13, L01-L04 truth
junctions, lit-v1.0 literature junctions) need columns chrom1, pos1, chrom2, pos2; a site
matches a truth junction when both ends lie within --tolerance (either orientation).
Writes OUT_summary.tsv, OUT_removed_sites.tsv, OUT_truth_hits.tsv, OUT_loci.tsv.
"""
import argparse
import bisect
import csv
from collections import defaultdict


def rows(path):
    with open(path) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


def key(r):
    return (r["chrom1"], int(r["pos1"]), r["chrom2"], int(r["pos2"]), r["svtype"])


class EndIndex:
    """Points (chrom, pos) -> payload, queried within a tolerance."""

    def __init__(self, items):
        self.d = defaultdict(list)
        for c, p, v in items:
            self.d[c].append((p, v))
        for c in self.d:
            self.d[c].sort(key=lambda x: x[0])
        self.pos = {c: [p for p, _ in v] for c, v in self.d.items()}

    def near(self, c, p, tol):
        ps = self.pos.get(c, [])
        i = bisect.bisect_left(ps, p - tol)
        out = []
        while i < len(ps) and ps[i] <= p + tol:
            out.append(self.d[c][i][1])
            i += 1
        return out


def pair_match(index, ends, a, b, tol):
    """ids of truth junctions whose two ends match (a, b) in either orientation."""
    def close(x, y):
        return x[0] == y[0] and abs(x[1] - y[1]) <= tol
    hits = set()
    for i in index.near(a[0], a[1], tol):
        e1, e2 = ends[i]
        if (close(e1, a) and close(e2, b)) or (close(e2, a) and close(e1, b)):
            hits.add(i)
    return hits


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--old", required=True)
    ap.add_argument("--new", required=True)
    ap.add_argument("--flagged", required=True)
    ap.add_argument("--truth", action="append", default=[])
    ap.add_argument("--tolerance", type=int, default=500)
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()
    tol = a.tolerance
    flagged = {(r["cohort"], r["sample_id"]) for r in rows(a.flagged)}
    old = {key(r): r for r in rows(f"{a.old}/panel_junction_carriers.tsv")}
    new = {key(r): r for r in rows(f"{a.new}/panel_junction_carriers.tsv")}
    nidx = EndIndex([(k[0], k[1], k) for k in new])
    removed, shifted = [], []
    for k, r in old.items():
        if k in new:
            continue
        cand = [n for n in nidx.near(k[0], k[1], tol) if n[4] == k[4] and n[2] == k[2] and abs(n[3] - k[3]) <= tol]
        (shifted if cand else removed).append(r)
    kept_changed = [(k, int(r["n_carriers"]), int(new[k]["n_carriers"])) for k, r in old.items()
                    if k in new and r["n_carriers"] != new[k]["n_carriers"]]
    # truth junctions
    truth_ends, truth_lab = {}, {}
    for spec in a.truth:
        lab, path = spec.split("=", 1)
        for i, r in enumerate(rows(path)):
            tid = (lab, i)
            truth_ends[tid] = ((r["chrom1"], int(r["pos1"])), (r["chrom2"], int(r["pos2"])))
            ev = r.get("event") or r.get("event_id") or r.get("event_uid") or lab
            truth_lab[tid] = f"{lab}:{ev}"
    tidx = EndIndex([(e[0], e[1], t) for t, (e, _) in truth_ends.items()] +
                    [(e[0], e[1], t) for t, (_, e) in truth_ends.items()])
    with open(f"{a.out_prefix}_truth_hits.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["change", "chrom1", "pos1", "chrom2", "pos2", "svtype", "n_carriers_old", "n_carriers_new",
                    "truth"])
        n_hit = defaultdict(int)
        changed_rows = [("removed", r, "0") for r in removed] + [("shifted", r, "NA") for r in shifted] + \
                       [("carriers_dropped", old[k], str(nn)) for k, o, nn in kept_changed if nn < o]
        for ch, r, nn in changed_rows:
            k = key(r)
            hits = pair_match(tidx, truth_ends, (k[0], k[1]), (k[2], k[3]), tol)
            if hits:
                n_hit[ch] += 1
                w.writerow([ch, *k, r["n_carriers"], nn, ",".join(sorted({truth_lab[h] for h in hits}))])
    with open(f"{a.out_prefix}_removed_sites.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["change", "chrom1", "pos1", "chrom2", "pos2", "svtype", "n_carriers", "carriers"])
        for ch, lst in (("removed", removed), ("shifted", shifted)):
            for r in lst:
                w.writerow([ch, *key(r), r["n_carriers"], r["carriers"]])
    # loci and calls
    lo = {(r["Chr"], r["Start"], r["End"]): r for r in rows(f"{a.old}/panel_cxsv_loci.tsv")}
    ln = {(r["Chr"], r["Start"], r["End"]): r for r in rows(f"{a.new}/panel_cxsv_loci.tsv")}
    calls = rows(f"{a.old}/panel_cxsv_calls.tsv")
    fcalls = [r for r in calls if (r["Cohort"], r["Sample_ID"]) in flagged]
    with open(f"{a.out_prefix}_loci.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["change", "Chr", "Start", "End", "N_Individuals_old", "N_Individuals_new"])
        for k, r in lo.items():
            if k not in ln:
                w.writerow(["absent_in_new", *k, r["N_Individuals"], 0])
        for k, r in ln.items():
            if k not in lo:
                w.writerow(["new_in_new", *k, 0, r["N_Individuals"]])
    summ = [("junction_sites_old", len(old)), ("junction_sites_new", len(new)),
            ("sites_removed", len(removed)), ("sites_shifted", len(shifted)),
            ("sites_carriers_dropped", sum(1 for _, o, n in kept_changed if n < o)),
            ("removed_sites_with_1_carrier", sum(1 for r in removed if r["n_carriers"] == "1")),
            ("cxsv_loci_old", len(lo)), ("cxsv_loci_new", len(ln)),
            ("loci_absent_in_new", sum(1 for k in lo if k not in ln)), ("loci_new_in_new", sum(1 for k in ln if k not in lo)),
            ("cxsv_calls_old", len(calls)), ("cxsv_calls_of_flagged_libraries", len(fcalls)),
            ("truth_hits_removed", n_hit["removed"]), ("truth_hits_shifted", n_hit["shifted"]),
            ("truth_hits_carriers_dropped", n_hit["carriers_dropped"])]
    with open(f"{a.out_prefix}_summary.tsv", "w") as fh:
        fh.write("metric\tvalue\n" + "".join(f"{k}\t{v}\n" for k, v in summ))
    for k, v in summ:
        print(f"{k}\t{v}")


if __name__ == "__main__":
    main()
