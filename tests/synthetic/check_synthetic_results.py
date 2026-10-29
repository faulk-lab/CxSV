#!/usr/bin/env python3
"""Check the synthetic end-to-end run against expected.json.

For each planted CxSV (COMMON, RARE, PRIVATE) the patient must have at least
one CxSV locus overlapping it, and every such locus must carry the expected
panel carrier count and status (breakpoint-level matching). Every panel
individual must be callable at those loci (25x synthetic coverage) and the
release must have verified against its manifest. PRIVATE must also
show the HOTSPOT DECOY's region-overlap carriers (Panel_Region_Carriers,
Panel_Hotspot) while staying private. The SIMPLE deletion's junction must be
matched to its panel carriers (private/junction_panel_status.tsv). The private
VCF must contain records only near the PRIVATE event. Prints one line per check
and exits 1 on any failure.

Usage: check_synthetic_results.py <expected.json> <patient_results_root> <run> <sample>
"""
import csv
import json
import sys
from pathlib import Path

SLOP = 500      # breakpoint slop for the junction check (panel_match.tolerance)


def main(expected_json, patient_root, run, sample):
    exp = json.loads(Path(expected_json).read_text())
    root = Path(patient_root)
    sdir = root / "samples" / run / sample
    failures = []

    def check(ok, msg):
        print(("PASS " if ok else "FAIL ") + msg)
        if not ok:
            failures.append(msg)

    compat = (root / "panel_compatibility.txt").read_text()
    check("status: COMPATIBLE" in compat, "patient fingerprint matches the panel release")
    rcheck = (root / "panel_release_check.txt").read_text()
    check("status: VERIFIED" in rcheck, "panel release verified against its manifest")

    with open(sdir / "private" / "cxsv_panel_status.tsv", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    check(len(rows) > 0, f"patient has CxSV loci ({len(rows)})")

    contig = exp["contig"]
    for name, want in exp["patient_expectations"].items():
        ev = exp["events"][name]
        hits = [r for r in rows if r["Chr"] == contig
                and int(r["Start"]) <= ev["end"] and int(r["End"]) >= ev["start"]]
        check(len(hits) > 0, f"{name}: patient CxSV locus overlaps {contig}:{ev['start']}-{ev['end']} "
                             f"({len(hits)} locus/loci)")
        for r in hits:
            check(int(r["Panel_Carriers"]) == want["panel_carriers"] and r["Panel_Status"] == want["status"],
                  f"{name}: locus {r['Chr']}:{r['Start']}-{r['End']} carriers={r['Panel_Carriers']} "
                  f"status={r['Panel_Status']} (want {want['panel_carriers']}/{want['status']})")
            check(r["Panel_Callable"] == r["Panel_Individuals"] == str(exp.get("n_panel", 3)),
                  f"{name}: every panel individual callable at the locus's junctions "
                  f"(Panel_Callable={r['Panel_Callable']}, Panel_Individuals={r['Panel_Individuals']})")
            if "region_carriers" in want:
                check(int(r["Panel_Region_Carriers"]) == want["region_carriers"]
                      and r["Panel_Hotspot"] == want["hotspot"],
                      f"{name}: hotspot decoy seen by region overlap: region carriers="
                      f"{r['Panel_Region_Carriers']} hotspot={r['Panel_Hotspot']} "
                      f"(want {want['region_carriers']}/{want['hotspot']}) yet status {r['Panel_Status']}")

    with open(sdir / "private" / "junction_panel_status.tsv", newline="") as fh:
        jrows = list(csv.DictReader(fh, delimiter="\t"))
    for name, want in exp.get("junction_expectations", {}).items():
        m = [j for j in jrows if j["svtype"] == want["svtype"] and j["chrom1"] == contig
             and abs(int(j["pos1"]) - want["pos1"]) <= SLOP and abs(int(j["pos2"]) - want["pos2"]) <= SLOP]
        check(len(m) == 1 and int(m[0]["panel_carriers"]) == want["panel_carriers"],
              f"{name}: junction {want['svtype']} {want['pos1']}-{want['pos2']} matched on breakpoints, "
              f"panel carriers {[j['panel_carriers'] for j in m]} (want [{want['panel_carriers']}])")

    pv = exp["events"]["PRIVATE"]
    recs = [l.split("\t") for l in (sdir / "private" / "private_cxsv.vcf").read_text().splitlines()
            if l and not l.startswith("#")]
    check(len(recs) > 0, f"private VCF has records ({len(recs)})")
    outside = [f"{r[0]}:{r[1]}" for r in recs
               if not (r[0] == contig and pv["start"] - 60_000 <= int(r[1]) <= pv["end"] + 60_000)]
    check(not outside, f"private VCF records all lie near the PRIVATE event (outside: {outside[:5]})")

    print(f"\n{len(failures)} failure(s)")
    return 1 if failures else 0


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    sys.exit(main(*sys.argv[1:]))
