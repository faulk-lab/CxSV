#!/usr/bin/env python3
"""Helpers for the patient intake (workflow/intake/Snakefile, scripts/intake_align.sh).

    intake_tool.py fastq-tags --count-out N.txt IN.fastq[.gz] ...   FASTQ -> stdout, comments
                                                                    reduced to MM/ML/MN/qs tags
    intake_tool.py mmi-check  INDEX.mmi [--fai hs1.fa.fai]          k=17 w=10 map-ont index
    intake_tool.py verify     --sample S --reads-in N --flagstat F --header H
                              [--mod-in N --mod-out N] --out intake.json
    intake_tool.py sheets     --intake-sheet SHEET --root DIR --patient-sheet OUT --summary OUT

All checks live in workflow/lib/cxsv_intake.py (unit-tested in tests/unit/test_intake.py).
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import cxsv_intake as ci  # noqa: E402


def cmd_fastq_tags(a):
    n = ci.normalise_fastq(a.inputs, sys.stdout)
    sys.stdout.flush()
    with open(a.count_out, "w") as fh:
        fh.write(f"{n}\n")
    print(f"[intake_tool fastq-tags] processed {n} reads", file=sys.stderr)


def cmd_mmi_check(a):
    n = None
    if a.fai:
        with open(a.fai) as fh:
            n = sum(1 for l in fh if l.strip())
    h = ci.check_mmi(a.index, n)
    print(json.dumps(h))


def cmd_verify(a):
    flag = open(a.flagstat).read()
    header = open(a.header).read()
    prov = {}
    if a.provenance and os.path.exists(a.provenance):
        prov = dict(l.rstrip("\n").split("\t", 1) for l in open(a.provenance) if "\t" in l)
    try:
        res = ci.verify(a.sample, int(a.reads_in), flag, header, a.mod_in, a.mod_out)
        rc = 0
    except ci.IntakeCheckError as e:
        fs = ci.parse_flagstat(flag) if "primary" in flag else {}
        res = {"sample_id": a.sample, "reads_in": int(a.reads_in), **{f"out_{k}": v for k, v in fs.items()},
               "verdict": "FAIL", "problems": [str(e)]}
        rc = 6
    res["recipe"] = {k: (list(v) if isinstance(v, tuple) else v) for k, v in ci.RECIPE.items()}
    res["provenance"] = prov
    with open(a.out, "w") as fh:
        json.dump(res, fh, indent=2)
    print(f"{a.sample}: {res['verdict']} reads_in={res['reads_in']} "
          f"primary_mapped={res.get('out_primary_mapped')} unmapped={res.get('out_unmapped')}"
          + ("" if rc == 0 else f" -- {res['problems'][0]}"), file=sys.stderr)
    sys.exit(rc)


def cmd_sheets(a):
    recs = ci.read_intake_sheet(a.intake_sheet, check_files=False)
    root = os.path.abspath(a.root)
    rows = ci.patient_sheet_rows(recs, lambda s: os.path.join(root, s, f"{s}.hs1.cram"))
    ci.write_tsv(rows, a.patient_sheet, ["sample_id", "alignment", "index", "individual_id", "platform"])
    summ = []
    for r in recs:
        j = json.load(open(os.path.join(root, r.sample_id, "intake.json")))
        summ.append({"sample_id": r.sample_id, "inputs": ",".join(r.inputs), "formats": r.format_summary,
                     "reads_in": j["reads_in"], "out_primary_mapped": j["out_primary_mapped"],
                     "out_unmapped": j["out_unmapped"], "out_supplementary": j.get("out_supplementary"),
                     "mod_records_in_sample": j.get("mod_records_in_sample"),
                     "mod_records_out_sample": j.get("mod_records_out_sample"),
                     "cram_md5": j.get("provenance", {}).get("cram_md5", ""), "verdict": j["verdict"]})
    ci.write_tsv(summ, a.summary, list(summ[0]) if summ else ["sample_id"])
    print(f"patient sheet: {a.patient_sheet} ({len(rows)} samples)", file=sys.stderr)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("fastq-tags"); s.add_argument("--count-out", required=True); s.add_argument("inputs", nargs="+")
    s.set_defaults(func=cmd_fastq_tags)
    s = sub.add_parser("mmi-check"); s.add_argument("index"); s.add_argument("--fai"); s.set_defaults(func=cmd_mmi_check)
    s = sub.add_parser("verify")
    for o in ("--sample", "--reads-in", "--flagstat", "--header", "--out"):
        s.add_argument(o, required=True)
    s.add_argument("--mod-in", type=int); s.add_argument("--mod-out", type=int); s.add_argument("--provenance")
    s.set_defaults(func=cmd_verify)
    s = sub.add_parser("sheets")
    for o in ("--intake-sheet", "--root", "--patient-sheet", "--summary"):
        s.add_argument(o, required=True)
    s.set_defaults(func=cmd_sheets)
    a = p.parse_args(argv)
    try:
        a.func(a)
    except (ci.IntakeCheckError, ci.IntakeSheetError) as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(6)


if __name__ == "__main__":
    main()
