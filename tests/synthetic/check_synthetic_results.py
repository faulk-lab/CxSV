#!/usr/bin/env python3
"""Check the synthetic end-to-end run against expected.json.

For each planted CxSV (COMMON, RARE, PRIVATE) the patient must have at least
one CxSV locus overlapping it, and every such locus must carry the expected
panel carrier count and public/private status. The private VCF must contain
records only near the PRIVATE event. Prints one line per check and exits 1 on
any failure.

Usage: check_synthetic_results.py <expected.json> <patient_results_root> <run> <sample>
"""
import csv
import json
import sys
from pathlib import Path

WINDOW = 1000   # same as panel_match.window


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
