#!/usr/bin/env python3
"""Copy-number PoN stratum of every panel library (config cn_strata_table; rules/cn.smk).

The 1kGP panel sample sheets carry no chemistry column, so the stratum (the library-QC
group names of item 14) comes from the lab registry: the cohort fixes it for the 1kGP-1019
ONT and PacBio cohorts; for LRSC ONT the registry chemistry (R9 / R10) decides.

    python3 make_panel_strata.py --registry sample_registry.tsv \
        --sheet 1kGP_ONT_1019_CHM13=<sheet.tsv> --sheet 1kGP_LRSC_ONT_CHM13=<sheet.tsv> \
        --sheet 1kGP_PacBio_CHM13=<sheet.tsv> [--check lofo_r2_per_genome.tsv] --out panel_strata_v1.tsv

--check compares with the groups of the item-14 calibration table (columns cohort,
sample_id, group) and fails on any difference.
"""
import argparse
import csv
import os
import sys

COHORT_STRATUM = {"1kGP_ONT_1019_CHM13": "ONT_1kGP_1019", "1kGP_PacBio_CHM13": "PacBio_HiFi"}
CHEMISTRY_STRATUM = {"R9": "ONT_LRSC_R9", "R10": "ONT_LRSC_R10"}
CHECK_GROUP = {"ONT 1kGP-1019": "ONT_1kGP_1019", "ONT LRSC R9": "ONT_LRSC_R9",
               "ONT LRSC R10": "ONT_LRSC_R10", "PacBio HiFi": "PacBio_HiFi"}


def read_tsv(p):
    with open(p, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--registry", required=True)
    ap.add_argument("--sheet", action="append", required=True, help="COHORT=sample_sheet.tsv")
    ap.add_argument("--check", default=None)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    chem = {r["rel_path"]: (r.get("chemistry") or "").strip() for r in read_tsv(a.registry)}
    rows = []
    for spec in a.sheet:
        cohort, path = spec.split("=", 1)
        for r in read_tsv(path):
            if cohort in COHORT_STRATUM:
                st = COHORT_STRATUM[cohort]
            else:
                aln = r["alignment"]
                rel = "/".join(aln.split("/")[-2:])
                c = chem.get(rel, "")
                if c not in CHEMISTRY_STRATUM:
                    sys.exit(f"{cohort}/{r['sample_id']}: registry chemistry {c!r} for {rel} has no stratum")
                st = CHEMISTRY_STRATUM[c]
            rows.append((cohort, r["sample_id"], st))
    if a.check:
        want = {(r["cohort"], r["sample_id"]): CHECK_GROUP[r["group"]] for r in read_tsv(a.check)}
        got = {(c, s): st for c, s, st in rows}
        diff = sorted(k for k in set(want) | set(got) if want.get(k) != got.get(k))
        if diff:
            sys.exit(f"{len(diff)} libraries differ from {a.check}: {diff[:10]}")
        print(f"check: {len(got)} libraries agree with {os.path.basename(a.check)}", file=sys.stderr)
    with open(a.out, "w") as fh:
        fh.write("cohort\tsample_id\tcn_stratum\n")
        fh.writelines(f"{c}\t{s}\t{st}\n" for c, s, st in sorted(rows))
    print(f"{len(rows)} libraries -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
