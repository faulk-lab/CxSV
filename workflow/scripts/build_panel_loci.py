#!/usr/bin/env python3
"""Aggregate per-sample CxSV calls into the background panel.

Each background genome's CxSVs are discovered separately (rules/per_sample.smk).
This script stacks them into one table and counts, per locus, how many distinct
*individuals* carry a CxSV there. One person sequenced in two cohorts counts
once.

Inputs
------
--manifest   TSV: cohort, sample_id, individual_id, platform, coverage_x, master_table

Outputs
-------
--calls      panel_cxsv_calls.tsv: one row per (sample, CxSV locus) — the table
             patient loci are matched against
--loci       panel_cxsv_loci.tsv: overlapping calls merged across samples, with
             carrier counts (summary and plotting only: merging overlapping
             intervals can chain neighbouring loci, so matching uses --calls)
--sharing    4-column BED (Chr, Start, End, N_Individuals) for
             plot_cxsv_population_sharing.R
--summary    panel_summary.tsv: samples, individuals and calls per cohort,
             per platform and in total — the denominators for carrier frequency

Matching patient loci against the panel by breakpoints (rather than interval
overlap) is roadmap item 3; this table already carries what it needs.
"""
import argparse
import csv
import sys
from collections import defaultdict

CALL_COLS = ["Chr", "Start", "End", "Size", "SV_Types", "Cmplx_Class", "CxSV_Criteria",
             "Has_Inter_Chromosomal", "Cohort", "Sample_ID", "Individual_ID", "Platform"]


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def locus_interval(row):
    """Same clamping as annotate_cxsv.py applies to the loci BED."""
    start, end = max(0, int(row["Start"])), int(row["End"])
    if end < start:
        start, end = end, start
    if end <= start:
        end = start + 1
    return start, end


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--calls", required=True)
    ap.add_argument("--loci", required=True)
    ap.add_argument("--sharing", required=True)
    ap.add_argument("--summary", required=True)
    a = ap.parse_args()

    samples = read_tsv(a.manifest)
    calls = []
    for s in samples:
        for row in read_tsv(s["master_table"]):
            start, end = locus_interval(row)
            calls.append({
                "Chr": row["Chr"], "Start": start, "End": end, "Size": end - start,
                "SV_Types": row.get("SV_Types", ""), "Cmplx_Class": row.get("Cmplx_Class", ""),
                "CxSV_Criteria": row.get("CxSV_Criteria", ""),
                "Has_Inter_Chromosomal": row.get("Has_Inter_Chromosomal", ""),
                "Cohort": s["cohort"], "Sample_ID": s["sample_id"],
                "Individual_ID": s["individual_id"], "Platform": s["platform"],
            })
    calls.sort(key=lambda c: (c["Chr"], c["Start"], c["End"], c["Individual_ID"]))
    with open(a.calls, "w", newline="") as fh:
        w = csv.DictWriter(fh, CALL_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(calls)

    # merged loci with carrier counts
    loci = []
    cur = None
    for c in calls:
        if cur and c["Chr"] == cur["Chr"] and c["Start"] <= cur["End"]:
            cur["End"] = max(cur["End"], c["End"])
            cur["calls"].append(c)
        else:
            cur = {"Chr": c["Chr"], "Start": c["Start"], "End": c["End"], "calls": [c]}
            loci.append(cur)
    with open(a.loci, "w") as lf, open(a.sharing, "w") as sf:
        lf.write("Chr\tStart\tEnd\tN_Individuals\tN_Calls\tIndividuals_By_Platform\tCohorts\n")
        for L in loci:
            inds = {c["Individual_ID"] for c in L["calls"]}
            by_plat = defaultdict(set)
            for c in L["calls"]:
                by_plat[c["Platform"]].add(c["Individual_ID"])
            plat = ";".join(f"{p}:{len(v)}" for p, v in sorted(by_plat.items()))
            cohorts = ",".join(sorted({c["Cohort"] for c in L["calls"]}))
            lf.write(f"{L['Chr']}\t{L['Start']}\t{L['End']}\t{len(inds)}\t{len(L['calls'])}\t{plat}\t{cohorts}\n")
            sf.write(f"{L['Chr']}\t{L['Start']}\t{L['End']}\t{len(inds)}\n")

    # denominators
    n_calls = defaultdict(int)
    for c in calls:
        n_calls[c["Cohort"]] += 1
    groups = defaultdict(lambda: {"samples": set(), "individuals": set(), "calls": 0})
    for s in samples:
        for key in (("cohort", s["cohort"]), ("platform", s["platform"]), ("total", "all")):
            g = groups[key]
            g["samples"].add((s["cohort"], s["sample_id"]))
            g["individuals"].add(s["individual_id"])
    for c in calls:
        for key in (("cohort", c["Cohort"]), ("platform", c["Platform"]), ("total", "all")):
            groups[key]["calls"] += 1
    with open(a.summary, "w") as fh:
        fh.write("level\tname\tn_samples\tn_individuals\tn_cxsv_calls\n")
        for (level, name) in sorted(groups, key=lambda k: ({"total": 0, "platform": 1, "cohort": 2}[k[0]], k[1])):
            g = groups[(level, name)]
            fh.write(f"{level}\t{name}\t{len(g['samples'])}\t{len(g['individuals'])}\t{g['calls']}\n")

    tot = groups[("total", "all")]
    print(f"panel: {len(tot['samples'])} samples, {len(tot['individuals'])} individuals, "
          f"{len(calls)} CxSV calls, {len(loci)} merged loci", file=sys.stderr)


if __name__ == "__main__":
    main()
