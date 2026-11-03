#!/usr/bin/env python3
"""Leave-out panel status of every panel genome's own CxSV loci (item 10 calibration).

For each panel sample and each mode (LOIO: its individual removed; LOFO: its
whole panel family removed; NONE: nothing removed, used only to check this
script against classify_patient_cxsv.py), write what classify_patient_cxsv.py
would write for that sample against the reduced panel:

  <out>/<cohort>/<sample>/status_<MODE>.tsv     = private/cxsv_panel_status.tsv
  <out>/<cohort>/<sample>/jstatus_<MODE>.tsv    = private/junction_panel_status.tsv

The matching code is the production code: cxsv_matching.classify_locus and
junction_carriers on the released panel_junction_carriers.tsv, through
gp_common.ExcludingIndex, which drops the left-out individuals from every
matched site's carriers (the site envelopes of the release are kept; see
validate_rebuild.py for the check against true rebuilds). Panel_Callable comes
from locus_callability.py (same rules, leave-out applied); Panel_Individuals is
the panel size minus the left-out individuals; region carriers (informational)
lose the left-out individuals too.

Columns not used by the scorer and not recomputed here: junction-level
panel_callable and panel_af ('NA' / '.'). Two extra columns are appended to the
locus table: Panel_Matched_Fraction (= Panel_Junctions_Seen / Panel_Junctions)
and Leaveout_Excluded (number of panel individuals removed).
"""
import argparse
import gc
import gzip
import importlib.util
import multiprocessing as mp
import os
import sys
import time
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gp_common import REPO, ExcludingIndex, excluded, load_families, read_tsv  # noqa: E402
import cxsv_junctions as cj  # noqa: E402
import cxsv_matching as cm   # noqa: E402

spec = importlib.util.spec_from_file_location(
    "classify_patient_cxsv", os.path.join(REPO, "workflow", "scripts", "classify_patient_cxsv.py"))
cp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cp)

EXTRA = ["Panel_Matched_Fraction", "Leaveout_Excluded"]
G = {}


def load_callable(path, wanted):
    out = defaultdict(dict)
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            key = (f[ix["cohort"]], f[ix["sample_id"]])
            if wanted is not None and key not in wanted:
                continue
            out[key][f[ix["locus_id"]]] = {m: f[ix[f"callable_{m}"]] for m in ("NONE", "LOIO", "LOFO")}
    return out


def one_sample(task):
    cohort, sample, individual = task
    a = G["args"]
    t0 = time.time()
    S = f"{a.samples_root}/{cohort}/{sample}"
    od = f"{a.out}/{cohort}/{sample}"
    os.makedirs(od, exist_ok=True)
    inds, n_panel = G["inds"], len(G["inds"])
    junctions = cj.read_table(f"{S}/junctions/junctions.tsv")
    by_locus = defaultdict(list)
    for j in junctions:
        for L in j["loci"]:
            by_locus[L].append(j)
    with open(f"{S}/cxsv/cxsv_master_table.tsv", newline="") as fh:
        import csv
        reader = csv.DictReader(fh, delimiter="\t")
        cols = reader.fieldnames + cp.STATUS_COLS + EXTRA
        rows = list(reader)
    callab = G["callable"].get((cohort, sample), {})
    summary = []
    for mode in a.modes:
        excl = excluded(individual, G["families"], mode) & set(inds)
        idx = ExcludingIndex(G["index"], excl)
        n_eff = n_panel - len(excl)
        tmpj = f"{od}/jstatus_{mode}.tsv.tmp"
        with open(tmpj, "w") as fj:
            fj.write("\t".join(cp.JUNCTION_OUT) + "\n")
            for j in junctions:
                carriers, hits = cm.junction_carriers(j, idx, a.tolerance, a.ins_len_ratio)
                row = [j["chrom1"], j["pos1"], j["chrom2"], j["pos2"], j["svtype"], j["svlen"], j["locus_id"],
                       j["supp_vec"], j["dosage"], len(hits), len(carriers), "NA",
                       f"{len(carriers) / n_eff:.6f}" if n_eff else "NA", ".",
                       "public" if len(carriers) >= a.min_carriers else "private"]
                fj.write("\t".join(map(str, row)) + "\n")
        os.replace(tmpj, f"{od}/jstatus_{mode}.tsv")
        counts = defaultdict(int)
        tmps = f"{od}/status_{mode}.tsv.tmp"
        import csv
        with open(tmps, "w", newline="") as fa:
            wa = csv.DictWriter(fa, cols, delimiter="\t", lineterminator="\n")
            wa.writeheader()
            for r0 in rows:
                r = dict(r0)
                start, end = cp.interval(r)
                lid = f"{r['Chr']}:{r['Start']}-{r['End']}"
                res = cm.classify_locus(by_locus.get(lid, []), idx, a.tolerance, a.ins_len_ratio,
                                        a.min_junction_fraction, a.min_carriers)
                reg = cp.region_carriers(G["region"], r["Chr"], start, end, a.window) - excl
                carriers = res["carriers"]
                r["Panel_Junctions"] = res["n_junctions"]
                r["Panel_Junctions_Seen"] = res["n_seen"]
                r["Panel_Carriers"] = len(carriers)
                c = callab.get(lid, {}).get(mode, "NA")
                r["Panel_Callable"] = c if by_locus.get(lid) else "NA"
                r["Panel_Individuals"] = n_eff
                r["Panel_Carrier_Freq"] = f"{len(carriers) / n_eff:.6f}" if n_eff else "NA"
                r["Panel_Carriers_By_Platform"] = cp.breakdown(carriers, inds, "platforms")
                r["Panel_Carriers_By_Superpopulation"] = cp.breakdown(carriers, inds, "superpopulation")
                r["Panel_Best_Junction_Fraction"] = f"{res['best_fraction']:.3f}"
                r["Panel_Region_Carriers"] = len(reg)
                r["Panel_Hotspot"] = "yes" if (len(reg) >= a.min_carriers and res["status"] != "public") else "no"
                r["Panel_Status"] = res["status"]
                r["Panel_Matched_Fraction"] = (f"{res['n_seen'] / res['n_junctions']:.4f}"
                                               if res["n_junctions"] else "0.0000")
                r["Leaveout_Excluded"] = len(excl)
                counts[res["status"]] += 1
                wa.writerow(r)
        os.replace(tmps, f"{od}/status_{mode}.tsv")
        summary.append(f"{mode}:" + "/".join(str(counts[s]) for s in cm.STATUSES) + f"(-{len(excl)})")
    return f"{cohort}/{sample} {len(rows)} loci {' '.join(summary)} {time.time() - t0:.1f}s"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--release-dir", required=True)
    ap.add_argument("--sample-manifest", required=True)
    ap.add_argument("--samples-root", required=True)
    ap.add_argument("--families", required=True)
    ap.add_argument("--locus-callable", required=True, help="locus_callability.py output")
    ap.add_argument("--samples", default=None, help="optional file: cohort<TAB>sample_id per line (subset)")
    ap.add_argument("--modes", default="LOIO,LOFO")
    ap.add_argument("--tolerance", type=int, default=500)
    ap.add_argument("--ins-len-ratio", type=float, default=0.5)
    ap.add_argument("--min-junction-fraction", type=float, default=0.75)
    ap.add_argument("--min-carriers", type=int, default=1)
    ap.add_argument("--window", type=int, default=1000)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    a.modes = a.modes.split(",")
    t0 = time.time()
    man = read_tsv(a.sample_manifest)
    tasks = [(r["cohort"], r["sample_id"], r["individual_id"]) for r in man]
    if a.samples:
        want = {tuple(l.split("\t")[:2]) for l in open(a.samples).read().splitlines() if l.strip()}
        tasks = [t for t in tasks if t[:2] in want]
    R = a.release_dir
    G["args"] = a
    G["families"] = load_families(a.families)
    G["inds"] = cp.load_panel_individuals(f"{R}/panel_individuals.tsv")
    G["index"] = cm.SiteIndex(cm.read_carriers(f"{R}/panel_junction_carriers.tsv"))
    G["region"] = cp.load_region_calls(f"{R}/panel_cxsv_calls.tsv")
    G["callable"] = load_callable(a.locus_callable, {t[:2] for t in tasks})
    print(f"loaded panel ({len(G['inds'])} individuals) in {time.time() - t0:.0f} s; {len(tasks)} samples",
          file=sys.stderr, flush=True)
    gc.collect()
    gc.freeze()
    with mp.get_context("fork").Pool(a.threads) as pool:
        for msg in pool.imap_unordered(one_sample, tasks, chunksize=1):
            print(msg, file=sys.stderr, flush=True)
    print(f"done {time.time() - t0:.0f} s", file=sys.stderr)


if __name__ == "__main__":
    main()
