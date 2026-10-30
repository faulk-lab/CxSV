#!/usr/bin/env python3
"""Score every patient CxSV locus with the configurable points system (roadmap item 5).

Reads the patient's CxSV loci with their panel status (private/cxsv_panel_status.tsv),
catalogue hits (private/cxsv_catalogue_hits.tsv), junctions, merged-VCF records,
caller VCFs (supporting reads), locus read depth (scoring/locus_depth.tsv) and the
annotation / scoring resources, computes every feature of
workflow/lib/cxsv_scoring.py and writes

  --out-tsv       one row per locus, ranked: total, tier, top-3 contributions and,
                  for every feature, <feature>_value and <feature>_points
  --out-features  long format: Locus_ID Rank feature group value points rule
  --out-report    human-readable per-call report (markdown)

Weights, tiers and patient information come from --params (JSON written by the
workflow from config/scoring.yaml). Thresholds are operational choices, not
validated cutoffs (docs/SCORING.md); the report says so.
"""
import argparse
import csv
import json
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import cxsv_classes as cc      # noqa: E402
import cxsv_junctions as cj    # noqa: E402
import cxsv_matching as cm     # noqa: E402
import cxsv_scoring as sc      # noqa: E402

CALLERS = ("cuteSV", "SVIM", "Sniffles2")      # SURVIVOR sample-column order


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def info_dict(s):
    d = {}
    for kv in s.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
    return d


def read_merged_vcf(path):
    """{record id: {chrom, pos, alt, svtype, caller_ids: {caller: [ids]}}}"""
    recs = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            info = info_dict(f[7])
            ids = {}
            if len(f) > 9:
                keys = f[8].split(":")
                if "ID" in keys:
                    k = keys.index("ID")
                    for ci, col in enumerate(f[9:12]):
                        v = col.split(":")
                        if k < len(v) and v[k] not in ("NaN", "NA", "."):
                            ids[CALLERS[ci]] = v[k].split(",")
            recs[f[2]] = {"chrom": f[0], "pos": int(f[1]), "alt": f[4],
                          "svtype": info.get("SVTYPE", "."), "caller_ids": ids}
    return recs


def read_support(paths):
    """{caller: {record id: supporting reads}} from the callers' filtered VCFs."""
    keys = {"Sniffles2": ("SUPPORT",), "cuteSV": ("RE",), "SVIM": ("SUPPORT",)}
    out = {}
    for caller, path in paths.items():
        d = {}
        if path and os.path.exists(path):
            with open(path) as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    f = line.split("\t", 8)
                    if len(f) < 8:
                        continue
                    info = info_dict(f[7])
                    for k in keys[caller]:
                        if k in info:
                            try:
                                d[f[2]] = int(float(info[k]))
                            except ValueError:
                                pass
                            break
        out[caller] = d
    return out


def load_dosage(path):
    if not path:
        return {}
    return {r["gene_symbol"]: {"hi": r["hi_score"] or None, "ts": r["ts_score"] or None} for r in read_tsv(path)}


def load_regions(path):
    """ClinGen curated regions lifted to hs1 (scoring/clingen_region_dosage.tsv,
    acquisition/resources/build_resources.py): [{chrom, start, end, isca_id,
    name, hi, ts}], 0-based half-open."""
    if not path:
        return []
    return [{"chrom": r["chrom"], "start": int(r["start"]), "end": int(r["end"]), "isca_id": r["isca_id"],
             "name": r["region_name"], "hi": r["hi_score"] or None, "ts": r["ts_score"] or None}
            for r in read_tsv(path)]


def load_constraint(path):
    if not path:
        return {}
    out = {}
    for r in read_tsv(path):
        try:
            out[r["gene_symbol"]] = float(r["loeuf"])
        except ValueError:
            pass
    return out


def load_hpo(path):
    g = defaultdict(set)
    if path:
        for r in read_tsv(path):
            g[r["gene_symbol"]].add(r["hpo_id"])
    return g


def load_omim(path):
    """Genes with an OMIM phenotype (genemap2.txt, mapping key 3) — licensed, optional."""
    genes = set()
    if not path:
        return genes
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) > 12 and "(3)" in f[12]:
                for g in (f[8], *f[6].split(", ")):
                    if g:
                        genes.add(g.strip())
    return genes


def resource_versions(manifest):
    if not manifest or not os.path.exists(manifest):
        return []
    return [(r["file"], r["source"], r["source_version"], r["sha256"][:12], r["retrieved_utc"])
            for r in read_tsv(manifest)]


def junction_obj(r):
    return cj.Junction(r["chrom1"], r["pos1"], r["chrom2"], r["pos2"], r["svtype"], r["svlen"],
                       r["supp_vec"], r["dosage"], r["record_ids"].split(","))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    for k in ("status", "catalogue-hits", "junctions", "junction-status", "vcf", "depth", "params",
              "genes", "exons", "repeats", "segdups", "fai", "out-tsv", "out-features", "out-report"):
        ap.add_argument(f"--{k}", required=True)
    ap.add_argument("--caller-vcf", action="append", default=[], help="CALLER=path (Sniffles2, cuteSV, SVIM)")
    ap.add_argument("--parent", action="append", default=[], help="ROLE=junctions.tsv of a sequenced parent")
    ap.add_argument("--catalogue-check", default=None)
    ap.add_argument("--sample", default="patient")
    a = ap.parse_args()

    P = json.load(open(a.params))
    W, res = P["weights"], P.get("scoring_resources", {})
    tol = int(P["collapse_distance"])
    cls_p = {**cc.DEFAULTS, "collapse_distance": tol, "composite_fraction": float(P["composite_fraction"])}
    order = cj.contig_order(a.fai)
    chrom_len = {}
    for line in open(a.fai):
        f = line.split("\t")
        chrom_len[f[0]] = int(f[1])

    status = read_tsv(a.status)
    cat_rows = {r["Locus_ID"]: r for r in read_tsv(a.catalogue_hits)}
    cat_names = [c[:-len("_Contains")] for c in (next(iter(cat_rows.values())).keys() if cat_rows else [])
                 if c.endswith("_Contains")]
    by_locus = defaultdict(list)
    for r in cj.read_table(a.junctions):
        for L in r["loci"]:
            by_locus[L].append(r)
    jstatus = read_tsv(a.junction_status)
    recs = read_merged_vcf(a.vcf)
    support = read_support(dict(x.split("=", 1) for x in a.caller_vcf))
    depth = defaultdict(list)
    for line in open(a.depth):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 7:
            s, e = int(f[1]), int(f[2])
            depth[f[3]].append((f[0], s, e, f[4], float(f[6]) / max(1, e - s)))
    genes = sc.IntervalIndex(sc.read_bed(a.genes, split_labels=True))
    exons = sc.IntervalIndex(sc.read_bed(a.exons))
    lowc_classes = set(W["artefact_context"]["low_complexity_classes"])
    reps = sc.IntervalIndex([r for r in sc.read_bed(a.repeats) if r[3] in lowc_classes])
    segdup_rows = sc.read_bed(a.segdups)
    segdups = sc.IntervalIndex(segdup_rows)
    par = [(c, s, e) for c, s, e, _ in sc.read_bed(res["par_bed"])] if res.get("par_bed") else []
    dosage = load_dosage(res.get("clingen_dosage"))
    constraint = load_constraint(res.get("constraint"))
    regions = load_regions(res.get("clingen_regions"))
    region_index = sc.IntervalIndex([(r["chrom"], r["start"], r["end"], i) for i, r in enumerate(regions)])
    hpo_genes = load_hpo(res.get("hpo_annotations")) if P["patient"].get("hpo_terms") else {}
    hpo_parents = sc.read_obo_parents(res["hpo_ontology"]) if hpo_genes and res.get("hpo_ontology") else {}
    omim_genes = load_omim(res.get("omim_genemap2"))
    parents = {}
    for spec in a.parent:
        role, path = spec.split("=", 1)
        parents[role] = cm.SiteIndex(cm.junction_point_sites(cj.read_table(path), role))

    # sex: given, or inferred from the chrX / autosome depth of the locus flanks
    flanks = [(c, s, e, d) for rows in depth.values() for c, s, e, k, d in rows if k.startswith("flank")]
    sex = (P["patient"].get("sex") or "").lower()
    if sex in ("male", "female"):
        sex_note = f"sex {sex} (given)"
    else:
        sex, ratio = sc.infer_sex(flanks, par, P["copy_number"]["sex_ratio_male"], P["copy_number"]["sex_ratio_female"])
        sex_note = f"sex {sex} (inferred, chrX/autosome depth {sc.fnum(ratio)})"

    # chromosome breakend sets for criterion 7: non-public junctions, segdups masked
    eff = {}
    mask = defaultdict(list)
    if P["clustering"].get("mask_segdups", True):
        for c, s, e, _ in segdup_rows:
            mask[c].append((s, e))
    chrom_bk = defaultdict(list)
    for r in jstatus:
        if r.get("panel_carriers", "0") not in ("0", "NA"):
            continue
        for c, p in ((r["chrom1"], int(r["pos1"])), (r["chrom2"], int(r["pos2"]))):
            chrom_bk[c].append(p)

    def eff_of(c):
        if c not in eff:
            eff[c] = sc.EffectiveCoords(chrom_len.get(c, 0), mask.get(c, []))
        return eff[c]

    def merge(ps):
        out = []
        for p in sorted(ps):
            if not out or p - out[-1] > tol:
                out.append(p)
        return out

    loci = []
    for row in status:
        L = f"{row['Chr']}:{row['Start']}-{row['End']}"
        jr = by_locus.get(L, [])
        kept, comp = cc.distinct_junctions([junction_obj(r) for r in jr], order, cls_p)
        bk = cc.breakends(kept, tol)
        orients = []
        reads = []
        for j in kept:
            o, best = set(), None
            for rid in j.record_ids:
                rec = recs.get(rid)
                if not rec:
                    continue
                o |= sc.record_orientation(rec["chrom"], rec["pos"], rec["alt"], rec["svtype"])
                for caller, ids in rec["caller_ids"].items():
                    for i in ids:
                        v = support.get(caller, {}).get(i)
                        if v is not None:
                            best = v if best is None else max(best, v)
            orients.append(o)
            reads.append(best)
        segs = sc.cn_segments(depth.get(L, []), sex, par, P["copy_number"])
        gains = [(c, s, e) for c, s, e, cn, st, pl, ch in (segs or []) if ch and st > pl]
        losses = [(c, s, e) for c, s, e, cn, st, pl, ch in (segs or []) if ch and st < pl]
        # affected genes: wholly lost / gained (depth), else broken by a breakend
        affected = {}
        for mode, ivs in (("lost", losses), ("gained", gains)):
            for c, s, e in ivs:
                for gs, ge, g in genes.overlapping(c, s, e):
                    if gs >= s - P["copy_number"]["margin_bp"] and ge <= e + P["copy_number"]["margin_bp"]:
                        affected.setdefault(g, mode)
        for c, p in bk:
            for gs, ge, g in genes.overlapping(c, p - 1, p):
                affected.setdefault(g, "broken")
        # curated ClinGen regions (item 9): runs of lost / gained segments, same depth calls
        margin = int(P["copy_number"]["margin_bp"])
        reg_hits = [h for mode, ivs in (("lost", losses), ("gained", gains))
                    for h in sc.region_hits(region_index, regions, sc.merge_runs(ivs, 2 * margin), mode, margin)]
        # exons: in the rearranged interval of each chromosome, or hit by a breakend
        ex = set()
        per_chrom = defaultdict(list)
        for c, p in bk:
            per_chrom[c].append(p)
        for c, ps in per_chrom.items():
            for s, e, g in exons.overlapping(c, min(ps) - 1, max(ps)):
                ex.add((c, s, e))
        n_genes = len(affected)
        # artefact context
        wbp = int(W["artefact_context"]["window_bp"])
        f_sd = sum(1 for c, p in bk if segdups.overlapping(c, p - wbp, p + wbp)) / len(bk) if bk else 0.0
        f_lc = sum(1 for c, p in bk if reps.overlapping(c, p - wbp, p + wbp)) / len(bk) if bk else 0.0
        # criterion 7
        chrom = row["Chr"]
        E = eff_of(chrom)
        lpos = [p for c, p in bk if c == chrom and not E.masked(p)]
        cpos = merge([p for p in chrom_bk.get(chrom, []) if not E.masked(p)] + lpos)
        clus = sc.clustering_test(lpos, cpos, E, int(P["clustering"]["n_simulations"]),
                                  sc.seed_of(L, int(P["clustering"]["seed"]))) if len(lpos) >= 3 else None
        # criterion 8
        intra = [(j, o) for j, o in zip(kept, orients) if j.chrom1 == j.chrom2]
        counts, p_bal = sc.join_balance([o for _, o in intra])
        fb = sc.foldbacks([j for j, _ in intra], [o for _, o in intra], tol,
                          int(P["foldback"]["max_span"]), gains)
        n_inv = sum(1 for _, o in intra if o & {("+", "+"), ("-", "-")})
        start, end = int(row["Start"]), int(row["End"])
        subtel = min(start, chrom_len.get(chrom, end) - end) < int(P["foldback"]["subtelomeric_bp"])
        # inheritance
        mjf = float(P["panel_match"]["min_junction_fraction"])
        par_res = {}
        for role, idx in parents.items():
            k, n, _ = cm.match_fraction(jr, idx, int(P["panel_match"]["tolerance"]),
                                        float(P["panel_match"]["ins_len_ratio"]))
            par_res[role] = bool(n) and k / n >= mjf - 1e-9
        loci.append(dict(L=L, row=row, kept=kept, comp=comp, bk=bk, orients=orients, reads=reads, segs=segs,
                         affected=affected, reg_hits=reg_hits, n_exons=len(ex), n_genes=n_genes, f_sd=f_sd,
                         f_lc=f_lc, clus=clus,
                         counts=counts, p_bal=p_bal, fb=fb, n_inv=n_inv, subtel=subtel, par_res=par_res))

    qs = sc.bh_qvalues([x["clus"]["p_scan"] if x["clus"] else None for x in loci])
    for x, q in zip(loci, qs):
        row = x["row"]
        feats = []
        pp = sc.f_panel_privacy(row["Panel_Status"], row["Panel_Carriers"], row["Panel_Individuals"],
                                row.get("Panel_Callable", "NA"), W["panel_privacy"])
        feats.append(pp)
        cr = cat_rows.get(x["L"], {})
        x["cat"] = cr
        cats = []
        for n in cat_names:
            af = cr.get(f"{n}_Max_AF", ".")
            cats.append((n, cr.get(f"{n}_Junctions_Matched", "0/0"), None if af in (".", "") else float(af),
                         cr.get(f"{n}_Contains") == "yes"))
        feats.append(sc.f_catalogue_frequency(cats, pp.points, W["catalogue_frequency"]))
        feats.append(sc.f_gene_dosage(x["affected"], dosage, constraint, W["gene_dosage"]))
        feats.append(sc.f_region_dosage(x["reg_hits"], feats[-1].points, W["region_dosage"], bool(regions)))
        feats.append(sc.f_exons(x["n_exons"], W["exons"]))
        feats.append(sc.f_gene_count(x["n_genes"], W["gene_count"]))
        feats.append(sc.f_no_genes(x["n_genes"], x["n_exons"], W["no_genes"]))
        feats.append(sc.f_class(row["Cmplx_Class"], W["complexity_class"]))
        feats.append(sc.f_size(int(row["End"]) - int(row["Start"]), W["size"]))
        feats.append(sc.f_caller_agreement([j.supp_vec for j in x["kept"]], W["caller_agreement"]))
        known = [r for r in x["reads"] if r is not None]
        feats.append(sc.f_read_support(min(known) if known else None, W["read_support"]))
        feats.append(sc.f_artefact(x["f_sd"], x["f_lc"], W["artefact_context"]))
        feats.append(sc.f_interleaving(int(row.get("Max_Interleaved", 0) or 0), W["ct_interleaving"]))
        feats.append(sc.f_breakpoint_tier(row.get("Complexity_Tier", "NA"), row.get("N_Breakends", "NA"),
                                          W["ct_breakpoint_tier"]))
        feats.append(sc.f_cn_oscillation(x["segs"], W["ct_cn_oscillation"], sex_note))
        feats.append(sc.f_clustering(x["clus"], q, W["ct_clustering"]))
        ct_sub = sum(f.points for f in feats if f.name in sc.CT_FEATURES)
        feats.append(sc.f_bfb(ct_sub, x["counts"], x["p_bal"], len(x["fb"]), sum(1 for _, g in x["fb"] if g),
                              x["n_inv"], x["subtel"], P["foldback"], W["ct_bfb_modifier"]))
        terms = P["patient"].get("hpo_terms") or []
        gene_terms = set().union(*(hpo_genes.get(g, set()) for g in x["affected"])) if hpo_genes else set()
        feats.append(sc.f_hpo(terms, sc.hpo_match(terms, gene_terms, hpo_parents) if terms else set(),
                              W["hpo_match"]))
        feats.append(sc.f_omim({g for g in x["affected"] if g in omim_genes}, bool(omim_genes), W["omim"]))
        feats.append(sc.f_inheritance(x["par_res"], bool(P["patient"].get("parents_affected")), W["inheritance"]))
        x["feats"] = feats
        x["ct"] = round(sum(f.points for f in feats if f.group == "chromothripsis"), 4)
        x["total"] = round(sum(f.points for f in feats), 4)
        x["tier"] = sc.tier_of(x["total"], P["tiers"])
        x["q"] = q

    loci.sort(key=lambda x: (-x["total"], x["L"]))
    write_outputs(a, P, loci, cat_names, sex_note)


def top3(feats):
    nz = sorted((f for f in feats if f.points), key=lambda f: (-abs(f.points), sc.FEATURE_ORDER.index(f.name)))
    return "; ".join(f"{f.name} {f.points:+.2f}" for f in nz[:3]) or "none"


def write_outputs(a, P, loci, cat_names, sex_note):
    base = ["Rank", "Locus_ID", "Chr", "Start", "End", "Size", "Cmplx_Class", "SV_Types", "CxSV_Criteria",
            "Complexity_Tier", "N_Junctions", "N_Breakends", "Max_Interleaved", "Panel_Status", "Panel_Carriers",
            "Panel_Callable", "Panel_Individuals", "Catalogues_Containing", "Catalogue_Max_AF"]
    cats = [f"{n}_{s}" for n in cat_names for s in ("Junctions_Matched", "Max_AF")]
    cols = base + cats + ["Genes_Affected", "Total_Score", "Tier", "Chromothripsis_Points", "Top_Features"]
    cols += [f"{n}_{k}" for n in sc.FEATURE_ORDER for k in ("value", "points")]
    with open(a.out_tsv, "w", newline="") as fh, open(a.out_features, "w", newline="") as ff:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        wf = csv.writer(ff, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        wf.writerow(["Locus_ID", "Rank", "feature", "group", "value", "points", "rule"])
        for rank, x in enumerate(loci, 1):
            x["rank"] = rank
            fm = {f.name: f for f in x["feats"]}
            out = {"Rank": rank, "Locus_ID": x["L"]}
            out.update({c: x["row"].get(c, "NA") for c in base[2:17]})
            out["Catalogues_Containing"] = x["cat"].get("Catalogues_Containing", "NA")
            out["Catalogue_Max_AF"] = x["cat"].get("Catalogue_Max_AF", "NA")
            out.update({c: x["cat"].get(c, "NA") for c in cats})
            out["Genes_Affected"] = ",".join(f"{g}:{m}" for g, m in sorted(x["affected"].items())) or "none"
            out["Total_Score"] = f"{x['total']:.2f}"
            out["Tier"] = x["tier"]
            out["Chromothripsis_Points"] = f"{x['ct']:.2f}"
            out["Top_Features"] = top3(x["feats"])
            for n in sc.FEATURE_ORDER:
                out[f"{n}_value"] = fm[n].value
                out[f"{n}_points"] = f"{fm[n].points:.2f}"
            w.writerow([out.get(c, "NA") for c in cols])
            for n in sc.FEATURE_ORDER:
                f = fm[n]
                wf.writerow([x["L"], rank, f.name, f.group, f.value, f"{f.points:.2f}", f.rule])
    write_report(a, P, loci, sex_note)


def write_report(a, P, loci, sex_note):
    L = []
    L.append(f"# CxSV score report: {a.sample}\n")
    L.append("Explicit points system modelled on ACMG/ClinGen CNV scoring plus computed chromothripsis "
             "criteria 1, 2, 5, 7 and 8 (docs/SCORING.md). **Every threshold and weight is an operational, "
             "cohort-specific choice, not a validated cutoff.** This report is research output for review by "
             "qualified clinical scientists; it is not a diagnosis.\n")
    L.append(f"- loci scored: {len(loci)}; {sex_note}")
    L.append(f"- panel release: {P.get('panel_release_id', 'NA')}; matching tolerance "
             f"{P['panel_match']['tolerance']} bp, min junction fraction {P['panel_match']['min_junction_fraction']}")
    vers = resource_versions(P.get("scoring_resources", {}).get("manifest"))
    if vers:
        L.append("- resources: " + "; ".join(f"{f} ({s} {v}, sha256 {h})" for f, s, v, h, _ in vers
                                             if not f.endswith(".bak")))
    if a.catalogue_check and os.path.exists(a.catalogue_check):
        cats = read_tsv(a.catalogue_check)
        L.append("- external catalogues: " + ("; ".join(f"{r['catalogue']} {r['version']} ({r['status']})"
                                                        for r in cats) or "none configured"))
    tiers = defaultdict(int)
    for x in loci:
        tiers[x["tier"]] += 1
    L.append("- tiers: " + ", ".join(f"{t['name']} {tiers[t['name']]}" for t in P["tiers"]) + "\n")
    L.append("## Ranking\n")
    L.append("| Rank | Locus | Class | Panel | Genes | Total | Tier | Top contributions |")
    L.append("|---:|---|---|---|---|---:|---|---|")
    for x in loci:
        g = ",".join(sorted(x["affected"])) or "-"
        L.append(f"| {x['rank']} | {x['L']} | {x['row']['Cmplx_Class']} | {x['row']['Panel_Status']} "
                 f"({x['row']['Panel_Carriers']}/{x['row']['Panel_Individuals']}) | {g} | {x['total']:.2f} | "
                 f"{x['tier']} | {top3(x['feats'])} |")
    for x in loci:
        r = x["row"]
        L.append(f"\n## {x['rank']}. {x['L']}: {r['Cmplx_Class']}, total {x['total']:.2f} ({x['tier']})\n")
        L.append(f"- junctions {r.get('N_Junctions', 'NA')} ({len(x['comp'])} composite records not counted), "
                 f"breakends {r.get('N_Breakends', 'NA')} ({r.get('Complexity_Tier', 'NA')}), "
                 f"types {r.get('SV_Types', 'NA')}, criteria {r.get('CxSV_Criteria', 'NA')}")
        L.append("- genes: " + (", ".join(f"{g} ({m})" for g, m in sorted(x["affected"].items())) or "none")
                 + f"; exons disrupted: {x['n_exons']}")
        if x["segs"]:
            L.append("- copy number: " + "; ".join(f"{c}:{s:,}-{e:,} CN {cn:.2f}" + (" (changed)" if ch else "")
                                                   for c, s, e, cn, st, pl, ch in x["segs"]))
        L.append("")
        L.append("| Feature | Value | Points | Rule |")
        L.append("|---|---|---:|---|")
        for f in x["feats"]:
            L.append(f"| {f.name} | {str(f.value).replace('|', '/')} | {f.points:+.2f} | {f.rule} |")
        L.append(f"| **total** | | **{x['total']:+.2f}** | chromothripsis points {x['ct']:+.2f} |")
    with open(a.out_report, "w") as fh:
        fh.write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
