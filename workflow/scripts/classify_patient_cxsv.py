#!/usr/bin/env python3
"""Label each patient CxSV locus private, shared_junctions or public against the panel.

Matching is on breakpoints, not on region overlap (roadmap item 3). Each patient
CxSV locus is the set of its constituent junctions (junctions/junctions.tsv,
locus_id column). A junction matches a panel junction site when svtype and both
contigs agree and each end lies within --tolerance bp of the site's coordinate
envelope (INS: lengths within --ins-len-ratio). A panel individual carries the
locus when it carries at least --min-junction-fraction of the locus's junctions
(workflow/lib/cxsv_matching.py: classify_locus). Then

  public            >= --min-carriers panel individuals carry the rearrangement
  shared_junctions  at least one junction occurs in the panel, but not the
                    rearrangement (a different event, or a common simple SV, at
                    the same place): the locus is NOT public
  private           no junction of the locus occurs in the panel

Region overlap is kept as information only: Panel_Region_Carriers counts the
panel individuals with any CxSV call within --window bp of the locus (the
previous rule), and Panel_Hotspot = yes when that count reaches --min-carriers
although the locus is not public — the case where region overlap would have
masked a private event.

Outputs
-------
--out-all        master table + Panel_Junctions, Panel_Junctions_Seen, Panel_Carriers,
                 Panel_Carrier_Freq, Panel_Carriers_By_Platform,
                 Panel_Carriers_By_Superpopulation, Panel_Best_Junction_Fraction,
                 Panel_Region_Carriers, Panel_Hotspot, Panel_Status
--out-private    rows that are not public (private and shared_junctions)
--out-bed        those loci as BED (Chr, Start, End) for VCF extraction
--out-junctions  every patient junction (whole SV set) with its panel carriers,
                 carrier frequency and best-matching site's allele frequency
"""
import argparse
import bisect
import csv
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import cxsv_junctions as cj  # noqa: E402
import cxsv_matching as cm   # noqa: E402

STATUS_COLS = ["Panel_Junctions", "Panel_Junctions_Seen", "Panel_Carriers", "Panel_Carrier_Freq",
               "Panel_Carriers_By_Platform", "Panel_Carriers_By_Superpopulation",
               "Panel_Best_Junction_Fraction", "Panel_Region_Carriers", "Panel_Hotspot", "Panel_Status"]
JUNCTION_OUT = ["chrom1", "pos1", "chrom2", "pos2", "svtype", "svlen", "locus_id", "supp_vec", "dosage",
                "panel_sites", "panel_carriers", "panel_carrier_freq", "panel_af", "panel_status"]


def interval(row):
    start, end = max(0, int(row["Start"])), int(row["End"])
    if end < start:
        start, end = end, start
    if end <= start:
        end = start + 1
    return start, end


# ─── region overlap (informational) ──────────────────────────────────────────
def load_region_calls(calls_tsv):
    by_chr = defaultdict(list)
    with open(calls_tsv, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            by_chr[r["Chr"]].append((int(r["Start"]), int(r["End"]), r["Individual_ID"]))
    index = {}
    for chrom, rows in by_chr.items():
        rows.sort()
        index[chrom] = ([r[0] for r in rows], rows, max(r[1] - r[0] for r in rows))
    return index


def region_carriers(index, chrom, start, end, window):
    if chrom not in index:
        return set()
    starts, rows, maxlen = index[chrom]
    lo, hi = start - window, end + window
    i = bisect.bisect_left(starts, lo - maxlen)
    found = set()
    while i < len(rows) and rows[i][0] <= hi:
        s, e, ind = rows[i]
        if e >= lo:
            found.add(ind)
        i += 1
    return found


def load_panel_individuals(path):
    with open(path, newline="") as fh:
        return {r["individual_id"]: r for r in csv.DictReader(fh, delimiter="\t")}


def load_site_af(path):
    if not path:
        return {}
    out = {}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[(r["svtype"], r["chrom1"], int(r["pos1"]), r["chrom2"], int(r["pos2"]),
                 cm.parse_svlen(r["svlen"]))] = r["af"]
    return out


def breakdown(carriers, inds, field):
    c = defaultdict(int)
    for i in carriers:
        vals = (inds.get(i, {}).get(field) or "unknown").split(",")
        for v in vals:
            c[v] += 1
    return ";".join(f"{k}:{v}" for k, v in sorted(c.items())) or "none"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--master-table", required=True)
    ap.add_argument("--junctions", required=True, help="patient junctions/junctions.tsv")
    ap.add_argument("--panel-carriers", required=True)
    ap.add_argument("--panel-sites", default=None)
    ap.add_argument("--panel-individuals", required=True)
    ap.add_argument("--panel-calls", required=True, help="panel_cxsv_calls.tsv (region overlap, informational)")
    ap.add_argument("--tolerance", type=int, default=500)
    ap.add_argument("--ins-len-ratio", type=float, default=0.5)
    ap.add_argument("--min-junction-fraction", type=float, default=0.75)
    ap.add_argument("--min-carriers", type=int, default=1)
    ap.add_argument("--window", type=int, default=1000)
    ap.add_argument("--out-all", required=True)
    ap.add_argument("--out-private", required=True)
    ap.add_argument("--out-bed", required=True)
    ap.add_argument("--out-junctions", required=True)
    a = ap.parse_args()

    inds = load_panel_individuals(a.panel_individuals)
    n_panel = len(inds)
    index = cm.SiteIndex(cm.read_carriers(a.panel_carriers))
    site_af = load_site_af(a.panel_sites)
    region = load_region_calls(a.panel_calls)
    junctions = cj.read_table(a.junctions)
    by_locus = defaultdict(list)
    for j in junctions:
        for L in j["loci"]:
            by_locus[L].append(j)

    # every junction: panel carriers (for simple SVs and for the locus table)
    with open(a.out_junctions, "w") as fj:
        fj.write("\t".join(JUNCTION_OUT) + "\n")
        for j in junctions:
            carriers, hits = cm.junction_carriers(j, index, a.tolerance, a.ins_len_ratio)
            afs = [site_af.get((s["svtype"], s["chrom1"], s["pos1"], s["chrom2"], s["pos2"], s["svlen"]), ".")
                   for s in hits]
            afs = [float(x) for x in afs if x not in (".", "")]
            row = [j["chrom1"], j["pos1"], j["chrom2"], j["pos2"], j["svtype"], j["svlen"], j["locus_id"],
                   j["supp_vec"], j["dosage"], len(hits), len(carriers),
                   f"{len(carriers) / n_panel:.6f}" if n_panel else "NA",
                   f"{max(afs):.6f}" if afs else ".",
                   "public" if len(carriers) >= a.min_carriers else "private"]
            fj.write("\t".join(map(str, row)) + "\n")

    with open(a.master_table, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols = reader.fieldnames + STATUS_COLS
        rows = list(reader)

    counts = defaultdict(int)
    with open(a.out_all, "w", newline="") as fa, open(a.out_private, "w", newline="") as fp, \
            open(a.out_bed, "w") as fb:
        wa = csv.DictWriter(fa, cols, delimiter="\t", lineterminator="\n")
        wp = csv.DictWriter(fp, cols, delimiter="\t", lineterminator="\n")
        wa.writeheader()
        wp.writeheader()
        for r in rows:
            start, end = interval(r)
            locus_id = f"{r['Chr']}:{r['Start']}-{r['End']}"
            res = cm.classify_locus(by_locus.get(locus_id, []), index, a.tolerance, a.ins_len_ratio,
                                    a.min_junction_fraction, a.min_carriers)
            reg = region_carriers(region, r["Chr"], start, end, a.window)
            carriers = res["carriers"]
            r["Panel_Junctions"] = res["n_junctions"]
            r["Panel_Junctions_Seen"] = res["n_seen"]
            r["Panel_Carriers"] = len(carriers)
            r["Panel_Carrier_Freq"] = f"{len(carriers) / n_panel:.6f}" if n_panel else "NA"
            r["Panel_Carriers_By_Platform"] = breakdown(carriers, inds, "platforms")
            r["Panel_Carriers_By_Superpopulation"] = breakdown(carriers, inds, "superpopulation")
            r["Panel_Best_Junction_Fraction"] = f"{res['best_fraction']:.3f}"
            r["Panel_Region_Carriers"] = len(reg)
            r["Panel_Hotspot"] = "yes" if (len(reg) >= a.min_carriers and res["status"] != "public") else "no"
            r["Panel_Status"] = res["status"]
            if res["n_junctions"] == 0:
                print(f"WARNING: locus {locus_id} has no junctions in {a.junctions}; labelled private",
                      file=sys.stderr)
            counts[res["status"]] += 1
            wa.writerow(r)
            if res["status"] != "public":
                wp.writerow(r)
                fb.write(f"{r['Chr']}\t{start}\t{end}\n")
    print(f"{a.master_table}: {len(rows)} CxSV loci: " + ", ".join(f"{counts[s]} {s}" for s in cm.STATUSES)
          + f" (panel of {n_panel} individuals; tolerance {a.tolerance} bp, min junction fraction "
          f"{a.min_junction_fraction}, min carriers {a.min_carriers})", file=sys.stderr)


if __name__ == "__main__":
    main()
