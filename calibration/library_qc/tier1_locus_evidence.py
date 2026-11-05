#!/usr/bin/env python3
"""Per-locus evidence of scored loci for the conservative-mode prototype
(docs/LIBRARY_QC.md): junction caller support and copy-number change.

    tier1_locus_evidence.py --out OUT.tsv [--tiers Tier1_high_priority,Tier2] \
        LABEL=SCORES_TSV=JUNCTIONS_TSV ...

SCORES_TSV is a score_cxsv.py table (cxsv_scores.tsv or a LOFO scores_LOFO.tsv);
JUNCTIONS_TSV the sample's junctions/junctions.tsv. For each locus of the listed
tiers: junctions (n), with >= 2 callers (n_multi), breakend-type junctions
(BND/TRA/INV; n_bnd), of them single-caller (n_bnd_single) and SVIM-only
(n_bnd_svim_only), and CN-changed segments parsed from ct_cn_oscillation_value.
"""
import argparse
import csv
import re
import sys
from collections import defaultdict

sys.path.insert(0, __file__.rsplit("/calibration/", 1)[0] + "/workflow/lib")
import cxsv_junctions as cj  # noqa: E402

BNDLIKE = {"BND", "TRA", "INV"}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tiers", default="Tier1_high_priority")
    ap.add_argument("specs", nargs="+")
    a = ap.parse_args()
    tiers = set(a.tiers.split(","))
    cols = ["label", "Locus_ID", "Cmplx_Class", "Panel_Status", "Total_Score", "Tier", "n", "n_multi",
            "n_bnd", "n_bnd_single", "n_bnd_svim_only", "cn_changed", "caller_agreement_value"]
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for spec in a.specs:
            label, scores, junc = spec.split("=", 2)
            rows = [r for r in csv.DictReader(open(scores), delimiter="\t") if r["Tier"] in tiers]
            if not rows:
                continue
            by = defaultdict(list)
            for j in cj.read_table(junc):
                for L in j["loci"]:
                    by[L].append(j)
            for r in rows:
                js = by.get(r["Locus_ID"], [])
                bnd = [j for j in js if j["svtype"] in BNDLIKE]
                m = re.search(r"(\d+) changed", r.get("ct_cn_oscillation_value", ""))
                w.writerow([label, r["Locus_ID"], r["Cmplx_Class"], r["Panel_Status"], r["Total_Score"], r["Tier"],
                            len(js), sum(j["supp_vec"].count("1") >= 2 for j in js), len(bnd),
                            sum(j["supp_vec"].count("1") == 1 for j in bnd),
                            sum(j["supp_vec"] == "010" for j in bnd), m.group(1) if m else "NA",
                            r.get("caller_agreement_value", "NA")])


if __name__ == "__main__":
    main()
