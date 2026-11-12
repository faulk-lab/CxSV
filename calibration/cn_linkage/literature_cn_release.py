#!/usr/bin/env python3
"""CN literature test A through the release pipeline (panel-v1.3 item d; docs/CN_LINKAGE.md
section 6a). Replaces the Phase 1 prototype path of literature_cn.py (calibration
cnlib/cnlink against phase1/seg and its 110-genome PoN) with the pipeline's own scripts and a
panel release's copy-number components.

Every lit-v1.0 junction-level event (63; 26 DUP-TRP/INV-DUP) is a synthetic patient, as in
tests/literature/score_literature_set.py (item 11, test A): its calls are exactly the published
junctions. Per event:

  junction track  score_literature_set.run_event, mode "main" (discovery as calling.yaml sets
                  it; its scores are the junction-level result)
  depth           the whole-genome 1 kb depth of a real healthy donor of the carrier's sex
                  (--donor-male / --donor-female, cn_depth_v1 tracks), multiplied inside the
                  event by (ploidy + delta) / ploidy from literature_cn_segments.tsv
                  (ploidy: cxsv_cn.ploidy_track with cn_params.par), written as the
                  mosdepth regions file the pipeline reads
  segments        workflow/scripts/cn_segment.py against <release>/cn/pon (the donor's
                  stratum from <release>/cn/pon/pon_members.tsv; sex called from the track,
                  sex-stratified chrX / chrY reference), cn_params of --repo (plus --cn-set)
  linkage         workflow/scripts/cn_linkage.py with <release>/cn/panel_cn_segments.tsv and
                  <release>/cn/pon/pon_members.tsv, --individual = the donor (its own panel
                  segments are left out, as for a patient who is not in the panel)
  loci, scoring   extract_junctions.py, annotate_cxsv.py, then the patient chain as
                  rules/patient_scoring.smk runs it (classify_patient_cxsv, match_catalogues,
                  locus_segments --cn-loci, score_cxsv --cn-only --cn-stats); locus depth is
                  synthesised from the curated copy number (score_literature_set.synth_depth),
                  as in item 11 test A and the prototype

Variants (--variants; segmentation is shared, linkage and scoring are per variant):
  v13      cn_params of --repo as committed
  v12rule  the same, with the panel-v1.2-r2 chrX/chrY public rule (public_min_genomes_sexchrom
           1, no sex-discordance exclusion), to separate D2 from the other changes

Outputs in --work: events/<event>/{main,release/<variant>}/..., literature_cn_release_scores.tsv
(per locus), literature_cn_release_summary.tsv (per event and variant; columns as
literature_cn.py's summary), literature_cn_release_counts.tsv (called / Tier 1 counts, next
to --prototype-summary when given).

Differences from the prototype that remain: the donors are members of the release PoN
(the prototype left the donor's family out of its PoN; one member of the donor stratum's
PoN), and only the donor's own panel libraries are left out of the public filter (the
prototype left out its family).
"""
import argparse
import concurrent.futures as cf
import csv
import gzip
import importlib.util
import json
import os
import shutil
import subprocess
import sys
from collections import defaultdict

G = {}
TIER1 = "Tier1_high_priority"


def load_mod(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def run(cmd, log):
    with open(log, "a") as fh:
        fh.write("$ " + " ".join(map(str, cmd)) + "\n")
        fh.flush()
        subprocess.run([str(c) for c in cmd], check=True, stdout=fh, stderr=fh)


def write_modelled_depth(src, out, cn_rows, sex, par, cn):
    """Donor regions.bed.gz -> out, depth x (1 + frac * delta / ploidy) per 1 kb bin overlapping
    a curated copy-number segment (frac: overlapped fraction of the bin)."""
    ev = defaultdict(list)
    for r in cn_rows:
        ev[r["chrom"]].append((int(r["start"]), int(r["end"]), float(r["delta"])))
    ploidy = {}
    with gzip.open(src, "rt") as fi, gzip.open(out, "wt", compresslevel=1) as fo:
        for line in fi:
            c = line.split("\t", 1)[0]
            if c not in ev:
                fo.write(line)
                continue
            _, s, e, d = line.rstrip("\n").split("\t")
            s, e, d = int(s), int(e), float(d)
            b = s // 1000
            if c not in ploidy:
                ploidy[c] = cn.ploidy_track(c, 400_000, sex, par, binsize=1000)
            pl = int(ploidy[c][min(b, len(ploidy[c]) - 1)])
            f = 1.0
            for a0, a1, dl in ev[c]:
                ov = min(e, a1) - max(s, a0)
                if ov > 0 and pl > 0:
                    f *= 1 + (ov / 1000.0) * dl / pl
            fo.write(f"{c}\t{s}\t{e}\t{d * f:.2f}\n" if f != 1.0 else line)


def score_chain(S, ev, cn_rows, sex, a, cfg, meta, callable_at):
    """rules/patient_scoring.smk on directory S (privacy -> catalogues -> locus depth -> score)."""
    slit = G["slit"]
    W = os.path.join(a.repo, "workflow", "scripts")
    cal, sc = cfg["calling"], cfg["scoring"]
    cp, jp = cal["cxsv_params"], cal["junction_params"]
    ann = cfg["annotation"]
    log = f"{S}/run.log"
    for d in ("private", "scoring"):
        os.makedirs(f"{S}/{d}", exist_ok=True)
    R = a.release
    m = cfg["patient"].get("panel_match", {})
    run([sys.executable, f"{W}/classify_patient_cxsv.py", "--master-table", f"{S}/cxsv/cxsv_master_table.tsv",
         "--junctions", f"{S}/junctions/junctions.tsv", "--panel-carriers", f"{R}/panel_junction_carriers.tsv",
         "--panel-sites", f"{R}/panel_junction_sites.tsv", "--panel-individuals", f"{R}/panel_individuals.tsv",
         "--panel-calls", f"{R}/panel_cxsv_calls.tsv", "--tolerance", m.get("tolerance", 500),
         "--ins-len-ratio", m.get("ins_len_ratio", 0.5), "--min-junction-fraction", m.get("min_junction_fraction", 0.75),
         "--min-carriers", m.get("min_carriers", 1), "--window", m.get("hotspot_window", 1000),
         "--out-all", f"{S}/private/cxsv_panel_status.tsv", "--out-private", f"{S}/private/private_cxsv_master_table.tsv",
         "--out-bed", f"{S}/private/private_cxsv_loci.bed", "--out-junctions", f"{S}/private/junction_panel_status.tsv"], log)
    slit.fill_callable(S, callable_at)
    cm_ = sc["external_catalogues"].get("match", {})
    cats = [x for n, v, t, mf in cfg["catalogues"] for x in ("--catalogue", f"{n}:{v}:{t}:{mf}")]
    run([sys.executable, f"{W}/match_catalogues.py", "--junctions", f"{S}/junctions/junctions.tsv",
         "--loci", f"{S}/private/cxsv_panel_status.tsv", *cats, "--tolerance", cm_.get("tolerance", 500),
         "--ins-len-ratio", cm_.get("ins_len_ratio", 0.5), "--min-junction-fraction", cm_.get("min_junction_fraction", 0.75),
         "--out-loci", f"{S}/private/cxsv_catalogue_hits.tsv", "--out-junctions", f"{S}/private/junction_catalogue_hits.tsv",
         "--out-check", f"{S}/private/catalogue_check.tsv"], log)
    cn = sc["copy_number"]
    run([sys.executable, f"{W}/locus_segments.py", "--junctions", f"{S}/junctions/junctions.tsv",
         "--loci", f"{S}/private/cxsv_panel_status.tsv", "--fai", a.fai, "--collapse-distance", jp["collapse_distance"],
         "--composite-fraction", cp["composite_fraction"], "--flank", cn["flank_bp"], "--gap", cn["gap_bp"],
         "--margin", cn["margin_bp"], "--cn-loci", f"{S}/cn/cn_loci.tsv", "--out", f"{S}/scoring/locus_segments.bed"], log)
    slit.synth_depth(f"{S}/scoring/locus_segments.bed", f"{S}/scoring/locus_depth.tsv", cn_rows, sex)
    with open(f"{S}/scoring/scoring_params.json", "w") as fh:
        json.dump(slit.scoring_params(cfg, meta, sex), fh, indent=2, sort_keys=True)
    for c in ("sniffles", "cutesv", "svim"):
        open(f"{S}/{c}_empty.vcf", "w").close()
    run([sys.executable, f"{W}/score_cxsv.py", "--status", f"{S}/private/cxsv_panel_status.tsv",
         "--catalogue-hits", f"{S}/private/cxsv_catalogue_hits.tsv", "--catalogue-check", f"{S}/private/catalogue_check.tsv",
         "--junctions", f"{S}/junctions/junctions.tsv", "--junction-status", f"{S}/private/junction_panel_status.tsv",
         "--vcf", f"{S}/cxsv/cxsv_only.vcf", "--caller-vcf", f"Sniffles2={S}/sniffles_empty.vcf",
         "--caller-vcf", f"cuteSV={S}/cutesv_empty.vcf", "--caller-vcf", f"SVIM={S}/svim_empty.vcf",
         "--depth", f"{S}/scoring/locus_depth.tsv", "--params", f"{S}/scoring/scoring_params.json",
         "--genes", ann["genes_bed"], "--exons", ann["exons_bed"], "--repeats", ann["repeats_bed"],
         "--segdups", ann["segdup_bed"], "--fai", a.fai, "--cn-only", f"{S}/cn/cn_only.tsv",
         "--cn-stats", f"{S}/cn/stats.json", "--sample", ev["event_id"], "--out-tsv", f"{S}/scoring/cxsv_scores.tsv",
         "--out-features", f"{S}/scoring/cxsv_score_features.tsv", "--out-report", f"{S}/scoring/cxsv_score_report.md"], log)
    return slit.read_tsv(f"{S}/scoring/cxsv_scores.tsv")


def variant_params(base, variant):
    p = json.loads(json.dumps(base))
    if variant == "v12rule":
        p["cn_params"]["linkage"].update(public_min_genomes_sexchrom=1, sexchrom_discordance_x=1e9,
                                         sexchrom_discordance_y=1e9)
    elif variant != "v13":
        raise ValueError(f"unknown variant {variant!r}")
    return p


def run_event(task):
    ev, junctions, cn_rows = task
    slit, a, cfg, meta, callable_at = G["slit"], G["a"], G["cfg"], G["meta"], G["callable_at"]
    W = os.path.join(a.repo, "workflow", "scripts")
    cal = cfg["calling"]
    jp, cp = cal["junction_params"], cal["cxsv_params"]
    main = slit.run_event((ev, junctions, cn_rows, "main", a, cfg, meta, callable_at))
    S0 = os.path.join(a.work, "events", ev["event_id"], "main")
    S1 = os.path.join(a.work, "events", ev["event_id"], "release")
    os.makedirs(f"{S1}/cn", exist_ok=True)
    log = f"{S1}/run.log"
    open(log, "w").close()
    sex = ev["carrier_sex"] if ev["carrier_sex"] in ("male", "female") else "female"
    donor = G["donor"][sex]
    write_modelled_depth(donor["depth"], f"{S1}/cn/depth.regions.bed.gz", cn_rows, sex, G["par"], G["cn"])
    with open(f"{S1}/cn_params.json", "w") as fh:
        json.dump(G["params"], fh, sort_keys=True)
    run([sys.executable, f"{W}/cn_segment.py", "--depth", f"{S1}/cn/depth.regions.bed.gz", "--pon", f"{a.release}/cn/pon",
         "--stratum", donor["stratum"], "--params", f"{S1}/cn_params.json",
         *[x for s in a.cn_set for x in ("--set", s)],
         "--out", f"{S1}/cn/segments.tsv", "--stats", f"{S1}/cn/stats.json"], log)
    os.remove(f"{S1}/cn/depth.regions.bed.gz")                 # 1,619-genome-sized intermediates are not kept
    contigs = {j["chrom1"] for j in junctions} | {j["chrom2"] for j in junctions} | {r["chrom"] for r in cn_rows}
    segs = [r for r in slit.read_tsv(f"{S1}/cn/segments.tsv") if r["chrom"] in contigs]
    stats = json.load(open(f"{S1}/cn/stats.json"))
    out = dict(event=ev, main=main, variants={}, stats=stats, contigs=sorted(contigs),
               segments=[f"{s['chrom']}:{s['start']}-{s['end']}:{s['state']}/{s['ploidy']}" for s in segs])
    for var in a.variants:
        S = f"{S1}/{var}"
        for d in ("cxsv", "junctions", "cn"):
            os.makedirs(f"{S}/{d}", exist_ok=True)
        for f in ("segments.tsv", "stats.json"):           # one segmentation, shared by the variants
            shutil.copyfile(f"{S1}/cn/{f}", f"{S}/cn/{f}")
        vlog = f"{S}/run.log"
        open(vlog, "w").close()
        with open(f"{S}/cn_params.json", "w") as fh:
            json.dump(variant_params(G["params"], var), fh, sort_keys=True)
        run([sys.executable, f"{W}/cn_linkage.py", "--vcf", f"{S0}/merged_sv.vcf", "--clusters", f"{S0}/cxsv/clusters.bed",
             "--summary", f"{S0}/cxsv/cxsv_summary.tsv", "--cxsv-vcf", f"{S0}/cxsv/cxsv_only.vcf",
             "--single", f"{S0}/cxsv/single_caller_filtered.tsv", "--segments", f"{S}/cn/segments.tsv",
             "--panel-segments", f"{a.release}/cn/panel_cn_segments.tsv",
             "--pon-members", f"{a.release}/cn/pon/pon_members.tsv" if var == "v13" else "none",
             "--individual", donor["individual"], "--stats", f"{S}/cn/stats.json", "--fai", a.fai,
             "--calling", f"{S}/cn_params.json", "--scripts", W, "--work", f"{S}/cn/work",
             "--out-clusters", f"{S}/cxsv/clusters.bed", "--out-summary", f"{S}/cxsv/cxsv_summary.tsv",
             "--out-vcf", f"{S}/cxsv/cxsv_only.vcf", "--out-candidates", f"{S}/cn/cn_candidates.tsv",
             "--out-loci", f"{S}/cn/cn_loci.tsv", "--out-cn-only", f"{S}/cn/cn_only.tsv"], vlog)
        ann = cfg["annotation"]
        run([sys.executable, f"{W}/extract_junctions.py", "--vcf", f"{S0}/merged_sv.vcf", "--clusters", f"{S}/cxsv/clusters.bed",
             "--summary", f"{S}/cxsv/cxsv_summary.tsv", "--fai", a.fai, "--collapse-distance", jp["collapse_distance"],
             "--ins-len-ratio", jp["ins_len_ratio"], "--out", f"{S}/junctions/junctions.tsv"], vlog)
        run([sys.executable, f"{W}/annotate_cxsv.py", "--summary", f"{S}/cxsv/cxsv_summary.tsv", "--vcf", f"{S}/cxsv/cxsv_only.vcf",
             "--genes", ann["genes_bed"], "--exons", ann["exons_bed"], "--repeats", ann["repeats_bed"],
             "--segdups", ann["segdup_bed"], "--out-table", f"{S}/cxsv/cxsv_master_table.tsv",
             "--out-bed", f"{S}/cxsv/cxsv_loci.bed"], vlog)
        n_loci = sum(r["IsCxSV"] == "YES" for r in slit.read_tsv(f"{S}/cxsv/cxsv_summary.tsv"))
        scores = score_chain(S, ev, cn_rows, sex, a, cfg, meta, callable_at) if n_loci else []
        loci = slit.read_tsv(f"{S}/cn/cn_loci.tsv")
        cands = slit.read_tsv(f"{S}/cn/cn_candidates.tsv")
        out["variants"][var] = dict(scores=scores, cn_ids={L["locus_id"] for L in loci if L["written"] == "yes"},
                                    cands=cands, loci=loci)
    return out


def summarise(results, variants):
    rows, summ = [], []
    for r in results:
        ev = r["event"]
        main_loci = r["main"]["loci"]
        best_main = max(main_loci, key=lambda L: float(L["Total_Score"]), default=None)
        for var in variants:
            v = r["variants"][var]
            sc = v["scores"]
            best = max(sc, key=lambda L: float(L["Total_Score"]), default=None)
            cnl = [L for L in sc if L["Locus_ID"] in v["cn_ids"]]
            best_cnl = max(cnl, key=lambda L: float(L["Total_Score"]), default=None)
            for L in sc:
                rows.append(dict(event_id=ev["event_id"], variant=var, author_event_type=ev["author_event_type"],
                                 cn_locus=L["Locus_ID"] in v["cn_ids"],
                                 **{k: L.get(k) for k in ("Locus_ID", "Size", "Cmplx_Class", "Panel_Status", "Genes_Affected",
                                                          "Total_Score", "Tier", "Top_Features")}))
            summ.append(dict(event_id=ev["event_id"], variant=var, author_event_type=ev["author_event_type"],
                             carrier_sex=ev["carrier_sex"], dup_trp=ev["author_event_type"] == "DUP-TRP/INV-DUP",
                             donor_sex_call=r["stats"].get("sex_call"), stratum=r["stats"].get("stratum"),
                             pon_sex_reference=r["stats"].get("pon_sex_reference"),
                             main_called=bool(main_loci), main_best_tier=best_main["Tier"] if best_main else "not_called",
                             main_best_score=best_main["Total_Score"] if best_main else "NA",
                             cn_called=bool(sc), cn_best_tier=best["Tier"] if best else "not_called",
                             cn_best_score=best["Total_Score"] if best else "NA",
                             cn_loci=len(cnl), cn_locus_classes=";".join(sorted({L["Cmplx_Class"] for L in cnl})),
                             cn_locus_best_tier=best_cnl["Tier"] if best_cnl else "NA",
                             cn_locus_best_score=best_cnl["Total_Score"] if best_cnl else "NA",
                             candidates=";".join(f"{c['decision']}:{c['cls']}:{c['n_segments']}seg:{c['n_records']}rec:"
                                                 f"{c['n_long']}long:pub={c['public_counts']}" for c in v["cands"]
                                                 if any(p.rsplit(":", 1)[0] in r["contigs"] for p in c["parts"].split(";"))),
                             segments=";".join(r["segments"])))
    return rows, summ


def counts(summ, label):
    def n(rows, f):
        return sum(1 for x in rows if f(x))
    out = []
    groups = {"all": summ, "dup_trp": [x for x in summ if str(x["dup_trp"]) == "True"],
              "other": [x for x in summ if str(x["dup_trp"]) != "True"]}
    for g, rs in groups.items():
        out.append(dict(source=label, group=g, events=len(rs),
                        junction_called=n(rs, lambda x: str(x["main_called"]) == "True"),
                        junction_tier1=n(rs, lambda x: x["main_best_tier"] == TIER1),
                        cn_called=n(rs, lambda x: str(x["cn_called"]) == "True"),
                        cn_tier1=n(rs, lambda x: x["cn_best_tier"] == TIER1),
                        cn_locus_written=n(rs, lambda x: int(x["cn_loci"] or 0) > 0),
                        cn_class_dup_trp=n(rs, lambda x: "DUP-TRP" in (x["cn_locus_classes"] or ""))))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--literature", required=True)
    ap.add_argument("--repo", required=True)
    ap.add_argument("--release", required=True, help="panel release (cn/pon, cn/panel_cn_segments.tsv, ...)")
    ap.add_argument("--resources", default="/mnt/sargasso1/CxSV/resources")
    ap.add_argument("--fai", default="/mnt/sargasso1/CxSV/refs/hs1.fa.fai")
    ap.add_argument("--depth-root", default="/mnt/sargasso1/CxSV/build/cn_depth_v1",
                    help="<root>/<cohort>/<sample>/depth.regions.bed.gz")
    ap.add_argument("--donor-male", default="1kGP_LRSC_ONT_CHM13/HG02026")
    ap.add_argument("--donor-female", default="1kGP_LRSC_ONT_CHM13/HG02181")
    ap.add_argument("--variants", default="v13,v12rule")
    ap.add_argument("--cn-set", action="append", default=[], help="cn_segment.py --set override (e.g. "
                    "segmentation.nested_min_z=5.5); repeatable")
    ap.add_argument("--prototype-summary", default=None, help="literature_cn.py summary to report next to")
    ap.add_argument("--work", required=True)
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--events", default=None)
    ap.add_argument("--sources", default="nazaryan2018,collins2017,eisfeldt2024,bahrambeigi2019,grochowski2024")
    a = ap.parse_args()
    a.variants = a.variants.split(",")
    os.makedirs(a.work, exist_ok=True)
    slit = load_mod("slit", f"{a.repo}/tests/literature/score_literature_set.py")
    sys.path.insert(0, f"{a.repo}/workflow/lib")
    import cxsv_cn as cn
    cfg = slit.load_config(a)
    cal = cfg["calling"]
    params = {"cn_params": cal["cn_params"], "cxsv_params": cal["cxsv_params"], "junction_params": cal["junction_params"]}
    members = {r["key"]: r for r in slit.read_tsv(f"{a.release}/cn/pon/pon_members.tsv")}
    donor = {}
    for sex, key in (("male", a.donor_male), ("female", a.donor_female)):
        m = members[key]
        donor[sex] = dict(key=key, individual=m["individual_id"], stratum=m["stratum"], in_pon=m["in_pon"],
                          sex_call=m["sex_call"], depth=f"{a.depth_root}/{key}/depth.regions.bed.gz")
        print(f"donor {sex}: {key} stratum {m['stratum']} sex_call {m['sex_call']} in_pon {m['in_pon']}", file=sys.stderr)
    events = [e for e in slit.read_tsv(f"{a.literature}/literature_events.tsv") if e["source_id"] in a.sources.split(",")]
    if a.events:
        events = [e for e in events if e["event_id"] in set(a.events.split(","))]
    junctions = [j for j in slit.read_tsv(f"{a.literature}/literature_junctions.tsv") if j["lift_status"] == "lifted"]
    cnr = [r for r in slit.read_tsv(f"{a.literature}/literature_cn_segments.tsv") if r["lift_status"] == "lifted"]
    by_j, by_cn = defaultdict(list), defaultdict(list)
    for j in junctions:
        by_j[j["event_id"]].append(j)
    for r in cnr:
        by_cn[r["event_id"]].append(r)
    mine = [j for j in junctions if j["event_id"] in {e["event_id"] for e in events}]
    callable_at, meta = slit.compute_callable(mine, a)
    G.update(slit=slit, a=a, cfg=cfg, meta=meta, callable_at=callable_at, donor=donor, cn=cn, params=params,
             par=cn.par_from_params(cal["cn_params"].get("par")))
    tasks = [(e, by_j[e["event_id"]], by_cn[e["event_id"]]) for e in events if by_j[e["event_id"]]]
    print(f"{len(tasks)} events, variants {a.variants}, release {meta.get('release_id')}", file=sys.stderr)
    results = []
    with cf.ProcessPoolExecutor(a.jobs) as ex:
        for r in ex.map(run_event, tasks):
            results.append(r)
            print(r["event"]["event_id"], {v: len(x["scores"]) for v, x in r["variants"].items()}, file=sys.stderr)
    rows, summ = summarise(results, a.variants)
    slit.write(f"{a.work}/literature_cn_release_scores.tsv", rows)
    slit.write(f"{a.work}/literature_cn_release_summary.tsv", summ)
    cnt = []
    if a.prototype_summary:
        cnt += counts(slit.read_tsv(a.prototype_summary), "prototype (literature_cn.py)")
    for var in a.variants:
        cnt += counts([x for x in summ if x["variant"] == var], f"release pipeline, {var}")
    slit.write(f"{a.work}/literature_cn_release_counts.tsv", cnt)
    with open(f"{a.work}/run_info.json", "w") as fh:
        json.dump(dict(release_id=meta.get("release_id"), release=a.release, repo=a.repo, donors=donor,
                       variants=a.variants, cn_set=a.cn_set, events=len(tasks)), fh, indent=1, default=str)
    for c in cnt:
        print("\t".join(str(c[k]) for k in c), file=sys.stderr)


if __name__ == "__main__":
    main()
