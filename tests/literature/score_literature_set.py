#!/usr/bin/env python3
"""Annotation / privacy / scoring test of the literature set (roadmap item 11, test A).

Every curated literature event is turned into a synthetic patient whose calls
are exactly the published junctions (no reads), and the patient workflow's own
scripts are run on it, rule by rule, with the same parameters:

  merged_sv.vcf        one record per lifted junction, SUPP_VEC 111 (all three callers),
                       ALT in breakend notation carrying the junction's sides, GT 0/1
                       (1/1 when the carrier is homozygous or hemizygous)
  vcf_to_bed.py -> bedtools sort/cluster -d cluster_distance -> classify_cxsv.py
                       (twice: --max-cluster-size as on main, and 0 = filter off, so an event
                       the oversized-cluster filter drops is still scored; item-13 dependency)
  extract_cxsv_vcf.py -> annotate_cxsv.py -> extract_junctions.py
  classify_patient_cxsv.py   against the panel-v1.0 release (verified separately);
                       Panel_Callable is filled from the release's callability masks,
                       computed once for all literature junctions (cxsv_callability)
  match_catalogues.py  the four external catalogues
  locus_segments.py    then read depth SYNTHESISED from the curated copy-number segments
                       (baseline 30x, copies from literature_cn_segments.tsv, carrier sex)
  score_cxsv.py        weights/tiers from config/scoring.yaml, no HPO terms, no parents

So the caller-evidence axes are idealised (all callers, no read-support
penalty); privacy, genes, class, size, artefact context and the chromothripsis
criteria are computed exactly as for a patient.

Separately, every lifted junction is matched against the panel junction sites and
each catalogue with the item-3 matcher (junction_privacy.tsv), whether or not
discovery makes a CxSV locus of it.

Outputs in --work:
  events/<event>/<mode>/...           all intermediate files
  literature_scores.tsv               one row per (event, mode, locus): class, status, axis
                                      points (privacy, genes, event, evidence, chromothripsis), total, tier
  literature_event_summary.tsv        one row per event: discovery (main), best tier main / filter-off
  junction_privacy.tsv                per junction: panel sites / carriers, catalogue hits + max AF

Usage (freya, env cxsv-pipeline):
  python tests/literature/score_literature_set.py --literature <lit dir> --repo <repo> \\
      --release /mnt/sargasso1/CxSV/build/panel-v1.0/results/panel/release --work <dir> [--jobs 6]
"""
import argparse
import concurrent.futures as cf
import csv
import json
import os
import subprocess
import sys
from collections import defaultdict

import yaml

MODES = {"main": None, "filter_off": 0}       # None = calling.yaml max_cluster_size
GROUPS = ("privacy", "genes", "event", "evidence", "chromothripsis", "phenotype", "inheritance")
ALT_FMT = {("+", "-"): "N[{c}:{p}[", ("+", "+"): "N]{c}:{p}]", ("-", "-"): "[{c}:{p}[N", ("-", "+"): "]{c}:{p}]N"}


def read_tsv(path):
    with open(path, newline="") as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def sh(cmd, log):
    with open(log, "a") as fh:
        fh.write(f"+ {cmd}\n")
        fh.flush()
        subprocess.run(cmd, shell=True, check=True, executable="/bin/bash", stdout=fh, stderr=fh)


# ─── synthetic inputs ────────────────────────────────────────────────────────
def write_vcf(path, junctions, contigs, gt):
    """SURVIVOR-like merged VCF: three sample columns (cuteSV, SVIM, Sniffles2)."""
    with open(path, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n##source=CxSV literature set (tests/literature/score_literature_set.py)\n")
        for c, n in contigs.items():
            fh.write(f"##contig=<ID={c},length={n}>\n")
        for k in ("SVTYPE", "END", "SVLEN", "CHR2", "SUPP_VEC"):
            fh.write(f'##INFO=<ID={k},Number=1,Type=String,Description="{k}">\n')
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write('##FORMAT=<ID=ID,Number=1,Type=String,Description="Caller record id">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tcuteSV\tSVIM\tSniffles2\n")
        for i, j in enumerate(sorted(junctions, key=lambda j: (j["chrom1"], int(j["pos1"]))), 1):
            rid = f"lit{i}"
            alt = ALT_FMT[(j["side1"], j["side2"])].format(c=j["chrom2"], p=j["pos2"])
            info = [f"SVTYPE={j['svtype']}", "SUPP_VEC=111"]
            if j["svtype"] == "BND":
                info.append(f"CHR2={j['chrom2']}")
            else:
                info += [f"END={j['pos2']}", f"SVLEN={abs(int(j['pos2']) - int(j['pos1']))}"]
            samples = "\t".join(f"{gt}:{rid}_{c}" for c in ("cuteSV", "SVIM", "Sniffles2"))
            fh.write(f"{j['chrom1']}\t{j['pos1']}\t{rid}\tN\t{alt}\t60\tPASS\t{';'.join(info)}\tGT:ID\t{samples}\n")


def ploidy(chrom, sex):
    if chrom == "chrY":
        return 1 if sex == "male" else 0
    if chrom == "chrX":
        return 1 if sex == "male" else 2
    return 2


def synth_depth(bed, out, cn_rows, sex, base=30.0):
    """locus_segments.bed -> bedcov-style rows with depth from the curated copy changes."""
    by = defaultdict(list)
    for r in cn_rows:
        by[r["chrom"]].append((int(r["start"]), int(r["end"]), float(r["delta"])))
    with open(bed) as fi, open(out, "w") as fo:
        for line in fi:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            c, s, e = f[0], int(f[1]), int(f[2])
            pl = ploidy(c, sex)
            if pl == 0 or e <= s:
                d = 0.0
            else:
                ov = sum((min(e, b) - max(s, a)) * dl for a, b, dl in by.get(c, ()) if min(e, b) > max(s, a))
                d = base / 2 * max(0.0, pl + ov / (e - s))
            fo.write("\t".join(f) + f"\t{int(round(d * (e - s)))}\n")


def scoring_params(cfg, release_meta, sex):
    sc, pat, cal = cfg["scoring"], cfg["patient"], cfg["calling"]
    m = pat.get("panel_match", {})
    return {"weights": sc["weights"], "tiers": sc["tiers"], "copy_number": sc["copy_number"],
            "clustering": sc["clustering"], "foldback": sc["foldback"],
            "scoring_resources": {k: v for k, v in sc["scoring_resources"].items() if v},
            "panel_match": {"tolerance": m.get("tolerance", 500), "ins_len_ratio": m.get("ins_len_ratio", 0.5),
                            "min_junction_fraction": m.get("min_junction_fraction", 0.75)},
            "collapse_distance": cal["junction_params"]["collapse_distance"],
            "composite_fraction": cal["cxsv_params"]["composite_fraction"],
            "panel_release_id": release_meta.get("release_id"),
            "patient": {"sex": sex, "hpo_terms": [], "parents_affected": False}}


# ─── one event ───────────────────────────────────────────────────────────────
def run_event(task):
    ev, junctions, cn_rows, mode, a, cfg, meta, callable_at = task
    S = os.path.join(a.work, "events", ev["event_id"], mode)
    for d in ("cxsv", "junctions", "private", "scoring"):
        os.makedirs(f"{S}/{d}", exist_ok=True)
    log = f"{S}/run.log"
    open(log, "w").close()
    W = os.path.join(a.repo, "workflow", "scripts")
    cal, sc = cfg["calling"], cfg["scoring"]
    cp, jp = cal["cxsv_params"], cal["junction_params"]
    ann = cfg["annotation"]
    sex = ev["carrier_sex"] if ev["carrier_sex"] in ("male", "female") else "female"
    gt = "1/1" if ev["zygosity"] in ("homozygous", "hemizygous") else "0/1"
    contigs = {ln.split("\t")[0]: int(ln.split("\t")[1]) for ln in open(a.fai)}
    write_vcf(f"{S}/merged_sv.vcf", junctions, contigs, gt)
    max_cluster = cp["max_cluster_size"] if MODES[mode] is None else MODES[mode]
    sh(f"python {W}/vcf_to_bed.py --vcf {S}/merged_sv.vcf --out {S}/cxsv/merged_sv.bed --window {cp['bnd_window']}", log)
    sh(f"bedtools sort -i {S}/cxsv/merged_sv.bed > {S}/cxsv/merged_sv_sorted.bed && "
       f"bedtools cluster -d {cp['cluster_distance']} -i {S}/cxsv/merged_sv_sorted.bed > {S}/cxsv/clusters.bed", log)
    sh(f"python {W}/classify_cxsv.py --bed {S}/cxsv/clusters.bed --vcf {S}/merged_sv.vcf --fai {a.fai} "
       f"--output {S}/cxsv/cxsv_summary.tsv --collapse-distance {jp['collapse_distance']} "
       f"--ins-len-ratio {jp['ins_len_ratio']} --composite-fraction {cp['composite_fraction']} "
       f"--min-bp {cp['min_breakpoints']} --max-window {cp['cluster_distance']} "
       f"--chromoanag-min {cp['chromoanag_min_window']} --chromoanag-max {cp['chromoanag_max_window']} "
       f"--chromoanag-min-sv {cp['chromoanag_min_svs']} --extreme-min-breakends {cp['extreme_min_breakends']} "
       f"--min-interleaved {cp['min_interleaved']} --interleave-alt-min-intra {cp['interleave_alt_min_intra']} "
       f"--interleave-alt-min-inter {cp['interleave_alt_min_inter']} --max-cluster-size {max_cluster}", log)
    sh(f"python {W}/extract_cxsv_vcf.py --vcf {S}/merged_sv.vcf --summary {S}/cxsv/cxsv_summary.tsv "
       f"--output {S}/cxsv/cxsv_only.vcf", log)
    sh(f"python {W}/annotate_cxsv.py --summary {S}/cxsv/cxsv_summary.tsv --vcf {S}/cxsv/cxsv_only.vcf "
       f"--genes {ann['genes_bed']} --exons {ann['exons_bed']} --repeats {ann['repeats_bed']} "
       f"--segdups {ann['segdup_bed']} --out-table {S}/cxsv/cxsv_master_table.tsv --out-bed {S}/cxsv/cxsv_loci.bed", log)
    sh(f"python {W}/extract_junctions.py --vcf {S}/merged_sv.vcf --clusters {S}/cxsv/clusters.bed "
       f"--summary {S}/cxsv/cxsv_summary.tsv --fai {a.fai} --collapse-distance {jp['collapse_distance']} "
       f"--ins-len-ratio {jp['ins_len_ratio']} --out {S}/junctions/junctions.tsv", log)
    summary = read_tsv(f"{S}/cxsv/cxsv_summary.tsv")
    out = {"event_id": ev["event_id"], "mode": mode, "clusters": len(summary),
           "oversized": sum(r["CxSV_Criteria"] == "FILTERED_oversized_cluster" for r in summary),
           "cxsv_loci": sum(r["IsCxSV"] == "YES" for r in summary),
           "max_cluster_bp": max((int(r["Size"]) for r in summary), default=0), "loci": []}
    if not out["cxsv_loci"]:
        return out
    R = a.release
    m = cfg["patient"].get("panel_match", {})
    sh(f"python {W}/classify_patient_cxsv.py --master-table {S}/cxsv/cxsv_master_table.tsv "
       f"--junctions {S}/junctions/junctions.tsv --panel-carriers {R}/panel_junction_carriers.tsv "
       f"--panel-sites {R}/panel_junction_sites.tsv --panel-individuals {R}/panel_individuals.tsv "
       f"--panel-calls {R}/panel_cxsv_calls.tsv --tolerance {m.get('tolerance', 500)} "
       f"--ins-len-ratio {m.get('ins_len_ratio', 0.5)} --min-junction-fraction {m.get('min_junction_fraction', 0.75)} "
       f"--min-carriers {m.get('min_carriers', 1)} --window {m.get('hotspot_window', 1000)} "
       f"--out-all {S}/private/cxsv_panel_status.tsv --out-private {S}/private/private_cxsv_master_table.tsv "
       f"--out-bed {S}/private/private_cxsv_loci.bed --out-junctions {S}/private/junction_panel_status.tsv", log)
    fill_callable(S, callable_at)
    cats = " ".join(f"--catalogue {n}:{v}:{t}:{mf}" for n, v, t, mf in cfg["catalogues"])
    cm_ = sc["external_catalogues"].get("match", {})
    sh(f"python {W}/match_catalogues.py --junctions {S}/junctions/junctions.tsv --loci {S}/private/cxsv_panel_status.tsv "
       f"{cats} --tolerance {cm_.get('tolerance', 500)} --ins-len-ratio {cm_.get('ins_len_ratio', 0.5)} "
       f"--min-junction-fraction {cm_.get('min_junction_fraction', 0.75)} --out-loci {S}/private/cxsv_catalogue_hits.tsv "
       f"--out-junctions {S}/private/junction_catalogue_hits.tsv --out-check {S}/private/catalogue_check.tsv", log)
    cn = sc["copy_number"]
    sh(f"python {W}/locus_segments.py --junctions {S}/junctions/junctions.tsv --loci {S}/private/cxsv_panel_status.tsv "
       f"--fai {a.fai} --collapse-distance {jp['collapse_distance']} --composite-fraction {cp['composite_fraction']} "
       f"--flank {cn['flank_bp']} --gap {cn['gap_bp']} --margin {cn['margin_bp']} --out {S}/scoring/locus_segments.bed", log)
    synth_depth(f"{S}/scoring/locus_segments.bed", f"{S}/scoring/locus_depth.tsv", cn_rows, sex)
    with open(f"{S}/scoring/scoring_params.json", "w") as fh:
        json.dump(scoring_params(cfg, meta, sex), fh, indent=2, sort_keys=True)
    for c in ("sniffles", "cutesv", "svim"):
        open(f"{S}/{c}_empty.vcf", "w").close()
    sh(f"python {W}/score_cxsv.py --status {S}/private/cxsv_panel_status.tsv --catalogue-hits {S}/private/cxsv_catalogue_hits.tsv "
       f"--catalogue-check {S}/private/catalogue_check.tsv --junctions {S}/junctions/junctions.tsv "
       f"--junction-status {S}/private/junction_panel_status.tsv --vcf {S}/cxsv/cxsv_only.vcf "
       f"--caller-vcf Sniffles2={S}/sniffles_empty.vcf --caller-vcf cuteSV={S}/cutesv_empty.vcf "
       f"--caller-vcf SVIM={S}/svim_empty.vcf --depth {S}/scoring/locus_depth.tsv --params {S}/scoring/scoring_params.json "
       f"--genes {ann['genes_bed']} --exons {ann['exons_bed']} --repeats {ann['repeats_bed']} --segdups {ann['segdup_bed']} "
       f"--fai {a.fai} --sample {ev['event_id']} --out-tsv {S}/scoring/cxsv_scores.tsv "
       f"--out-features {S}/scoring/cxsv_score_features.tsv --out-report {S}/scoring/cxsv_score_report.md", log)
    pts = defaultdict(lambda: defaultdict(float))
    for r in read_tsv(f"{S}/scoring/cxsv_score_features.tsv"):
        pts[r["Locus_ID"]][r["group"]] += float(r["points"])
    status = {f"{s['Chr']}:{s['Start']}-{s['End']}": s for s in read_tsv(f"{S}/private/cxsv_panel_status.tsv")}
    for r in read_tsv(f"{S}/scoring/cxsv_scores.tsv"):
        row = {k: r.get(k, "NA") for k in ("Locus_ID", "Size", "Cmplx_Class", "SV_Types", "N_Junctions", "N_Breakends",
                                          "Max_Interleaved", "Panel_Status", "Panel_Carriers", "Panel_Callable",
                                          "Panel_Individuals", "Catalogues_Containing", "Catalogue_Max_AF",
                                          "Genes_Affected", "Total_Score", "Tier", "Top_Features")}
        row.update({f"pts_{g}": round(pts[r["Locus_ID"]][g], 2) for g in GROUPS})
        st = status.get(r["Locus_ID"], {})
        row["Panel_Junctions_Seen"] = st.get("Panel_Junctions_Seen", "NA")
        row["Panel_Junctions"] = st.get("Panel_Junctions", "NA")
        out["loci"].append(row)
    return out


def fill_callable(S, callable_at):
    """Panel_Callable per locus = individuals callable at every junction (classify_patient_cxsv.n_callable)."""
    by = defaultdict(list)
    for j in read_tsv(f"{S}/junctions/junctions.tsv"):
        key = (j["chrom1"], int(j["pos1"]), j["chrom2"], int(j["pos2"]))
        for L in (j["locus_id"].split(",") if j["locus_id"] != "." else []):
            by[L].append(callable_at.get(key))
    path = f"{S}/private/cxsv_panel_status.tsv"
    rows = read_tsv(path)
    for r in rows:
        sets = by.get(f"{r['Chr']}:{r['Start']}-{r['End']}", [])
        r["Panel_Callable"] = len(set.intersection(*sets)) if sets and all(s is not None for s in sets) else "NA"
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


# ─── junction-level privacy ──────────────────────────────────────────────────
def junction_privacy(junctions, a, cfg):
    sys.path.insert(0, os.path.join(a.repo, "workflow", "lib"))
    import cxsv_matching as cm
    tol = 500
    index = cm.SiteIndex(cm.read_carriers(f"{a.release}/panel_junction_carriers.tsv"))
    q = defaultdict(list)
    for j in junctions:
        q[j["chrom1"]].append(int(j["pos1"]))
    for v in q.values():
        v.sort()
    cat_idx = {n: cm.SiteIndex(cm.read_catalogue_sites(t, q, tol)) for n, v, t, mf in cfg["catalogues"]}
    rows = []
    for j in junctions:
        jj = dict(j, pos1=int(j["pos1"]), pos2=int(j["pos2"]),
                  svlen=int(j["svlen"]) if j["svlen"] not in ("", ".") else 0)
        car, hits = cm.junction_carriers(jj, index, tol, 0.5)
        r = {k: j[k] for k in ("event_id", "junction_id", "svtype", "chrom1", "pos1", "chrom2", "pos2", "precision")}
        r.update(panel_sites=len(hits), panel_carriers=len(car))
        for n, idx in cat_idx.items():
            h = idx.query(jj, tol, 0.5)
            afs = [s["af"] for s in h if s.get("af") is not None]
            r[f"{n}_sites"] = len(h)
            r[f"{n}_max_af"] = f"{max(afs):.4g}" if afs else "."
        rows.append(r)
    return rows


def compute_callable(junctions, a):
    sys.path.insert(0, os.path.join(a.repo, "workflow", "lib"))
    import cxsv_callability as cc
    meta = json.load(open(f"{a.release}/release.json"))
    flank = int(meta["callability"]["flank_bp"])
    keys = sorted({(j["chrom1"], int(j["pos1"]), j["chrom2"], int(j["pos2"])) for j in junctions})
    sets = cc.callable_individuals(f"{a.release}/callability",
                                   [dict(chrom1=k[0], pos1=k[1], chrom2=k[2], pos2=k[3]) for k in keys], flank)
    return dict(zip(keys, sets)), meta


def load_config(a):
    C = os.path.join(a.repo, "config")
    cal = yaml.safe_load(open(f"{C}/calling.yaml"))
    sc = yaml.safe_load(open(f"{C}/scoring.yaml"))
    pat = yaml.safe_load(open(f"{C}/patient.yaml"))
    ann = {k: os.path.join(a.resources, v) for k, v in (("genes_bed", "genes.final.bed"),
                                                         ("exons_bed", "exons.filtered.bed"),
                                                         ("repeats_bed", "repeats.final.bed"),
                                                         ("segdup_bed", "segdup.final.bed"))}
    root = sc["external_catalogues"]["root"]
    cats = [(t["name"], t["version"],
             f"{root}/{t['name']}/{t['version']}/{t['name']}.{t['version']}.junction_sites.tsv.gz",
             f"{root}/{t['name']}/{t['version']}/{t['name']}.{t['version']}.manifest.tsv")
            for t in sc["external_catalogues"]["tables"]]
    return {"calling": cal, "scoring": sc, "patient": pat, "annotation": ann, "catalogues": cats}


def write(path, rows):
    cols = []
    for r in rows:
        for k in r:
            if k not in cols:
                cols.append(k)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n", restval="NA")
        w.writeheader()
        w.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--literature", required=True, help="literature set directory (literature_*.tsv)")
    ap.add_argument("--repo", required=True)
    ap.add_argument("--release", required=True)
    ap.add_argument("--resources", default="/mnt/sargasso1/CxSV/resources")
    ap.add_argument("--fai", default="/mnt/sargasso1/CxSV/refs/hs1.fa.fai")
    ap.add_argument("--work", required=True)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--events", default=None, help="comma-separated subset")
    ap.add_argument("--sources", default="nazaryan2018,collins2017,eisfeldt2024,bahrambeigi2019,grochowski2024",
                    help="sources with junction-level events to score (ClinVar is region level)")
    a = ap.parse_args()
    os.makedirs(a.work, exist_ok=True)
    cfg = load_config(a)
    events = [e for e in read_tsv(f"{a.literature}/literature_events.tsv") if e["source_id"] in a.sources.split(",")]
    if a.events:
        keep = set(a.events.split(","))
        events = [e for e in events if e["event_id"] in keep]
    junctions = [j for j in read_tsv(f"{a.literature}/literature_junctions.tsv") if j["lift_status"] == "lifted"]
    cn = [r for r in read_tsv(f"{a.literature}/literature_cn_segments.tsv") if r["lift_status"] == "lifted"]
    by_j, by_cn = defaultdict(list), defaultdict(list)
    for j in junctions:
        by_j[j["event_id"]].append(j)
    for r in cn:
        by_cn[r["event_id"]].append(r)
    ids = {e["event_id"] for e in events}
    mine = [j for j in junctions if j["event_id"] in ids]
    print(f"{len(events)} events, {len(mine)} lifted junctions", file=sys.stderr)
    write(f"{a.work}/junction_privacy.tsv", junction_privacy(mine, a, cfg))
    callable_at, meta = compute_callable(mine, a)
    tasks = [(e, by_j[e["event_id"]], by_cn[e["event_id"]], m, a, cfg, meta, callable_at)
             for e in events if by_j[e["event_id"]] for m in MODES]
    results = []
    with cf.ProcessPoolExecutor(a.jobs) as ex:
        for res in ex.map(run_event, tasks):
            results.append(res)
            print(f"{res['event_id']} {res['mode']}: {res['cxsv_loci']} CxSV loci, {res['oversized']} oversized",
                  file=sys.stderr)
    rows, summ = [], {}
    for res in results:
        for L in res["loci"] or [{}]:
            rows.append(dict({k: res[k] for k in ("event_id", "mode", "clusters", "oversized", "cxsv_loci",
                                                 "max_cluster_bp")}, **L))
        s = summ.setdefault(res["event_id"], {"event_id": res["event_id"]})
        best = max(res["loci"], key=lambda L: float(L["Total_Score"]), default=None)
        s[f"{res['mode']}_cxsv_loci"] = res["cxsv_loci"]
        s[f"{res['mode']}_oversized_clusters"] = res["oversized"]
        s[f"{res['mode']}_best_score"] = best["Total_Score"] if best else "NA"
        s[f"{res['mode']}_best_tier"] = best["Tier"] if best else "not_called"
        s[f"{res['mode']}_best_class"] = best["Cmplx_Class"] if best else "NA"
        s[f"{res['mode']}_best_status"] = best["Panel_Status"] if best else "NA"
    write(f"{a.work}/literature_scores.tsv", rows)
    write(f"{a.work}/literature_event_summary.tsv", list(summ.values()))


if __name__ == "__main__":
    main()
