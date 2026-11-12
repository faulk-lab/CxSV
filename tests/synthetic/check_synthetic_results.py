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
VCF must contain records only near the non-public events (PRIVATE, PARTIAL). PARTIAL (junctions seen in
the panel, rearrangement carried by nobody) must be shared_junctions and get
exactly the graded privacy points of roadmap item 10 for its matched fraction. LONGINV, a 600 kb
inversion spanning RARE and PRIVATE, must not hide them in an oversized cluster, and its own
breakends must not be CxSV loci (roadmap item 13). Prints one line per check
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

    # roadmap item 13: a long span record is clustered at its breakends, so nothing it
    # spans is dropped as an oversized cluster and its breakends alone are not CxSVs
    with open(sdir / "cxsv" / "cxsv_summary.tsv", newline="") as fh:
        clusters = list(csv.DictReader(fh, delimiter="\t"))
    over = [f"{r['Chr']}:{r['Start']}-{r['End']}" for r in clusters if r["CxSV_Criteria"] == "FILTERED_oversized_cluster"]
    for name, want in exp.get("long_span_expectations", {}).items():
        check(not over, f"{name}: no oversized cluster in the patient (found: {over[:3]})")
        for b in want["breakends"]:
            at = [r for r in rows if r["Chr"] == contig and int(r["Start"]) - want["window"] <= b <= int(r["End"]) + want["window"]]
            check(bool(at) == want["cxsv_at_breakends"],
                  f"{name}: breakend {b} {'is' if at else 'is not'} in a CxSV locus "
                  f"(want {'yes' if want['cxsv_at_breakends'] else 'no'})")
        for inner in want["inside"]:
            ev = exp["events"][inner]
            hits = [r for r in rows if r["Chr"] == contig and int(r["Start"]) <= ev["end"] and int(r["End"]) >= ev["start"]]
            check(bool(hits) and all(int(r["End"]) - int(r["Start"]) < 100_000 for r in hits),
                  f"{name}: {inner} inside it is found as its own locus "
                  f"({[r['Chr'] + ':' + r['Start'] + '-' + r['End'] for r in hits]})")

    with open(sdir / "private" / "junction_panel_status.tsv", newline="") as fh:
        jrows = list(csv.DictReader(fh, delimiter="\t"))
    for name, want in exp.get("junction_expectations", {}).items():
        m = [j for j in jrows if j["svtype"] == want["svtype"] and j["chrom1"] == contig
             and abs(int(j["pos1"]) - want["pos1"]) <= SLOP and abs(int(j["pos2"]) - want["pos2"]) <= SLOP]
        check(len(m) == 1 and int(m[0]["panel_carriers"]) == want["panel_carriers"],
              f"{name}: junction {want['svtype']} {want['pos1']}-{want['pos2']} matched on breakpoints, "
              f"panel carriers {[j['panel_carriers'] for j in m]} (want [{want['panel_carriers']}])")

    # scoring (roadmap items 4-5)
    with open(sdir / "scoring" / "cxsv_scores.tsv", newline="") as fh:
        scores = list(csv.DictReader(fh, delimiter="\t"))
    check(len(scores) == len(rows), f"every CxSV locus scored ({len(scores)} of {len(rows)})")
    for r in scores:
        tot = sum(float(v) for k, v in r.items() if k.endswith("_points"))
        check(abs(tot - float(r["Total_Score"])) < 0.02,
              f"{r['Locus_ID']}: feature points add up to the total ({tot:.2f} vs {r['Total_Score']})")

    def event_scores(name):
        ev = exp["events"][name]
        return [r for r in scores if r["Chr"] == contig and int(r["Start"]) <= ev["end"] and int(r["End"]) >= ev["start"]]

    tier_rank = {t: i for i, t in enumerate(["Tier1_high_priority", "Tier2", "Tier3", "Benign_likely"])}
    for name, want in exp.get("scoring_expectations", {}).items():
        hits = event_scores(name)
        check(bool(hits), f"{name}: scored locus present")
        for r in hits:
            if "rank" in want:
                check(r["Rank"] == str(want["rank"]) and r["Tier"] == want["tier"],
                      f"{name}: rank {r['Rank']} tier {r['Tier']} total {r['Total_Score']} "
                      f"(want rank {want['rank']}, {want['tier']}); top: {r['Top_Features']}")
            if want.get("hpo_match"):
                check(float(r["hpo_match_points"]) > 0, f"{name}: HPO match {r['hpo_match_value']}")
            if want.get("graded_privacy"):
                P = json.loads((sdir / "scoring" / "scoring_params.json").read_text())
                pw = P["weights"]["panel_privacy"]
                g = pw.get("graded") or {}
                st = next(x for x in rows if f"{x['Chr']}:{x['Start']}-{x['End']}" == r["Locus_ID"])
                k, n = int(st["Panel_Junctions_Seen"]), int(st["Panel_Junctions"])
                f = k / n
                if g.get("shape", "ramp") == "ramp":
                    lo, hi = float(g["f_lo"]), float(g["f_hi"])
                    w = 0.0 if f <= lo + 1e-12 else 1.0 if f >= hi - 1e-12 else (f - lo) / (hi - lo)
                else:
                    w = f ** float(g["k"])
                base, pub = float(pw["shared_junctions"]), float(pw["public"])
                want_pts = round(base + (pub - base) * w, 4) if g.get("enabled") else base
                check(0 < k < n and abs(float(r["panel_privacy_points"]) - want_pts) < 0.005 + 1e-9  # TSV: 2 dp
                      and "junctions seen" in r["panel_privacy_value"],
                      f"{name}: graded privacy {r['panel_privacy_points']} for {k}/{n} junctions seen "
                      f"(want {want_pts:.4f}; {r['panel_privacy_value']})")
            if "catalogue" in want:
                check(want["catalogue"] in r["Catalogues_Containing"].split(","),
                      f"{name}: contained in {r['Catalogues_Containing']} (want {want['catalogue']}; "
                      f"matched {r.get(want['catalogue'] + '_Junctions_Matched')})")
            if "tier_below" in want:
                ref = event_scores(want["tier_below"])
                check(bool(ref) and all(float(r["Total_Score"]) < float(x["Total_Score"])
                                        and tier_rank[r["Tier"]] > tier_rank[x["Tier"]] for x in ref),
                      f"{name}: total {r['Total_Score']} ({r['Tier']}) below {want['tier_below']} "
                      f"({', '.join(x['Total_Score'] + ' ' + x['Tier'] for x in ref)})")
    # panel-v1.3 (roadmap item 16): sample QC section of the score report, clinician report rule
    md = (sdir / "scoring" / "cxsv_score_report.md").read_text()
    check("## Sample QC and sex check" in md and "sex_call " in md and "x_ratio " in md,
          "score report has the sample QC and sex check section")
    rep = sdir / "report"
    pdf, docx, png = (rep / f"{sample}.clinician_report.pdf", rep / f"{sample}.clinician_report.docx",
                      rep / f"{sample}.tier1_points.png")
    check(pdf.exists() and pdf.read_bytes()[:5] == b"%PDF-", f"clinician report PDF written ({pdf.name})")
    import zipfile
    ok = docx.exists() and zipfile.is_zipfile(docx) and "word/document.xml" in zipfile.ZipFile(docx).namelist()
    check(ok, f"clinician report DOCX written ({docx.name})")
    if ok:
        xml = zipfile.ZipFile(docx).read("word/document.xml").decode("utf-8", "replace")
        n1 = sum(1 for r in scores if r["Tier"] == "Tier1_high_priority")
        check("Sample quality" in xml and "What the tiers mean" in xml and "Sex-chromosome check" in xml
              and (n1 == 0 or "word/media/" in " ".join(zipfile.ZipFile(docx).namelist())),
              f"clinician report has its sections and, with {n1} Tier 1 loci, the Tier 1 points graphic")
    check(png.exists() and png.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n", f"Tier 1 points graphic written ({png.name})")
    chk = (sdir / "private" / "catalogue_check.tsv").read_text()
    check("VERIFIED" in chk and "MISMATCH" not in chk, "catalogue table verified against its manifest")

    # the private VCF holds every non-public locus: PRIVATE, and PARTIAL (shared_junctions)
    near = [exp["events"][e] for e in ("PRIVATE", "PARTIAL") if e in exp["events"]]
    recs = [l.split("\t") for l in (sdir / "private" / "private_cxsv.vcf").read_text().splitlines()
            if l and not l.startswith("#")]
    check(len(recs) > 0, f"private VCF has records ({len(recs)})")
    outside = [f"{r[0]}:{r[1]}" for r in recs
               if not (r[0] == contig and any(ev["start"] - 60_000 <= int(r[1]) <= ev["end"] + 60_000 for ev in near))]
    check(not outside, f"private VCF records all lie near the non-public events PRIVATE / PARTIAL "
                       f"(outside: {outside[:5]})")

    print(f"\n{len(failures)} failure(s)")
    return 1 if failures else 0


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    sys.exit(main(*sys.argv[1:]))
