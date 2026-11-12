#!/usr/bin/env python3
"""Checks of the synthetic intake run (tests/synthetic/check_intake.sh). Exit 0 = all pass."""
import json
import sys
from pathlib import Path

import pysam

W = Path(sys.argv[1])
REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "workflow" / "lib"))
import cxsv_inputs  # noqa: E402

D, R = W / "data", W / "results" / "intake"
expected = json.load(open(D / "expected.json"))
SOURCES = {"S_modbam": ["modbam"], "S_fastq": ["fastq"], "S_multi": ["other1", "other2"], "S_cram": ["othercram"]}
failures = []


def check(cond, msg):
    print(("PASS " if cond else "FAIL ") + msg)
    if not cond:
        failures.append(msg)


for sample, sources in SOURCES.items():
    exp = {k: v for s in sources for k, v in expected[s].items()}
    j = json.load(open(R / sample / "intake.json"))
    check(j["verdict"] == "PASS", f"{sample}: intake verdict {j['verdict']} {j.get('problems')}")
    check(j["reads_in"] == len(exp), f"{sample}: reads_in {j['reads_in']} == {len(exp)} input reads")
    cram = pysam.AlignmentFile(str(R / sample / f"{sample}.hs1.cram"), reference_filename=str(D / "ref" / "ref.fa"))
    hdr = cram.header.to_dict()
    check([s["SN"] for s in hdr["SQ"]] == ["chrA", "chrB"], f"{sample}: contigs are the intake reference's (chrA, chrB)")
    check([rg["SM"] for rg in hdr.get("RG", [])] == [sample], f"{sample}: one read group, SM = {sample}")
    seen, bad_tags, bad_rg, no_qs, n_unmapped, n_supp = {}, 0, 0, 0, 0, 0
    for a in cram.fetch(until_eof=True):
        if a.is_supplementary or a.is_secondary:
            n_supp += a.is_supplementary
            continue
        n_unmapped += a.is_unmapped
        seen[a.query_name] = seen.get(a.query_name, 0) + 1
        e = exp.get(a.query_name)
        if e is None or a.get_tag("MM") != e["mm"] or list(a.get_tag("ML")) != e["ml"]:
            bad_tags += 1
        if a.get_tag("RG") != sample:
            bad_rg += 1
        if not a.has_tag("qs"):
            no_qs += 1
    check(set(seen) == set(exp) and set(seen.values()) == {1}, f"{sample}: every input read exactly once as primary/unmapped ({len(seen)})")
    check(bad_tags == 0, f"{sample}: MM/ML identical to the input for every read ({bad_tags} differ)")
    check(bad_rg == 0, f"{sample}: per-read RG is the sample (dorado RG not copied; {bad_rg} wrong)")
    check(no_qs == 0, f"{sample}: qs tag kept ({no_qs} missing)")
    check(j["out_unmapped"] == n_unmapped, f"{sample}: intake.json unmapped {j['out_unmapped']} == {n_unmapped}")
    if sample in ("S_modbam", "S_cram"):
        check(n_unmapped > 0, f"{sample}: unmappable reads kept as unmapped records ({n_unmapped})")
    if sample == "S_multi":
        check(n_supp > 0, f"{sample}: chimeric reads give supplementary records ({n_supp}) yet count once")
    pg = [p for p in hdr["PG"] if p.get("PN") == "minimap2"][0]["CL"]
    check(all(t in pg.split() for t in ("-y", "-ax", "map-ont", "--eqx", "-k", "17")), f"{sample}: @PG has the recipe: {pg[:80]}")

table, notes = cxsv_inputs.build_sample_table({"intake_run": {"samples": str(R / "patient_sheet.tsv")}}, "bamfiles")
check(sorted(s for _, s in table) == sorted(SOURCES), f"patient_sheet.tsv loads as a patient run: {notes}")
check(table[("intake_run", "S_multi")].individual_id == "indiv_multi", "individual_id carried into the patient sheet")
summary = (R / "intake_summary.tsv").read_text().splitlines()
check(len(summary) == 1 + len(SOURCES) and all(l.endswith("PASS") for l in summary[1:]), "intake_summary.tsv: one PASS line per sample")
if failures:
    print(f"{len(failures)} intake check(s) FAILED")
    sys.exit(1)
print("all intake result checks passed")
