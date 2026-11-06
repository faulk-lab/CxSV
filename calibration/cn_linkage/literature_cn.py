#!/usr/bin/env python3
"""Test A: the junction-level literature events with copy number modelled
(roadmap item 15, Phase 1; docs/CN_LINKAGE_PHASE1.md).

Every lit-v1.0 junction-level event (63; 26 DUP-TRP/INV-DUP) is a synthetic patient as in
tests/literature/score_literature_set.py (item 11, test A): its calls are exactly the
published junctions. Mode "main" is that test's discovery on main (calling.yaml as is)
and is re-run here against panel-v1.1 (baseline). Mode "cn" adds the copy-number track:

  depth      a real healthy genome of the carrier's sex (DONOR_<sex>, an unflagged
             LRSC R10 subset genome near 30x) gives the 1 kb depth track: its noise,
             biases and own CNVs; inside the event the depth is multiplied by
             (ploidy + delta) / ploidy from literature_cn_segments.tsv
  segments   cnlib against the donor stratum's PoN with the donor's family left out
             (only the event's contigs are segmented)
  linkage    cnlink.candidates on the event's records; panel frequency of the segments
             against the other-family subset genomes (segments_all.tsv)
  loci       augment.augment, then the same privacy / catalogue / depth / scoring chain as
             test A (depth synthesised from the curated CN, panel-v1.1 release)

Outputs in --work: events/<event>/{main,cn}/..., literature_cn_scores.tsv (per locus),
literature_cn_summary.tsv (per event: main vs cn, called, best tier/score/class, CN loci,
candidate decisions).
"""
import argparse
import concurrent.futures as cf
import csv
import importlib.util
import json
import os
import sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import augment  # noqa: E402
import cnlib  # noqa: E402
import cnlink  # noqa: E402

G = {}


def load_slit(repo):
    spec = importlib.util.spec_from_file_location("slit", f"{repo}/tests/literature/score_literature_set.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def modelled_track(fine, cn_rows, sex, par):
    """Donor 1 kb depth x (ploidy + delta) / ploidy over the curated CN segments."""
    out = {c: a.copy() for c, a in fine.items()}
    for r in cn_rows:
        c = r["chrom"]
        if c not in out:
            continue
        s, e, d = int(r["start"]), int(r["end"]), float(r["delta"])
        b0, b1 = s // 1000, -(-e // 1000)
        n = len(out[c])
        pl = cnlib.ploidy_track(c, n, sex, par, binsize=1000).astype(float)
        for b in range(max(0, b0), min(n, b1)):
            frac = (min(e, (b + 1) * 1000) - max(s, b * 1000)) / 1000.0
            if pl[b] > 0 and frac > 0:
                out[c][b] *= 1 + frac * d / pl[b]
    return out


def score_chain(S, ev, cn_rows, sex, slit, a, cfg, meta, callable_at):
    """The second half of slit.run_event (privacy -> score) on directory S."""
    W = os.path.join(a.repo, "workflow", "scripts")
    cal, sc = cfg["calling"], cfg["scoring"]
    cp, jp = cal["cxsv_params"], cal["junction_params"]
    ann = cfg["annotation"]
    log = f"{S}/run.log"
    for d in ("private", "scoring"):
        os.makedirs(f"{S}/{d}", exist_ok=True)
    R = a.release
    m = cfg["patient"].get("panel_match", {})
    slit.sh(f"python {W}/classify_patient_cxsv.py --master-table {S}/cxsv/cxsv_master_table.tsv "
            f"--junctions {S}/junctions/junctions.tsv --panel-carriers {R}/panel_junction_carriers.tsv "
            f"--panel-sites {R}/panel_junction_sites.tsv --panel-individuals {R}/panel_individuals.tsv "
            f"--panel-calls {R}/panel_cxsv_calls.tsv --tolerance {m.get('tolerance', 500)} "
            f"--ins-len-ratio {m.get('ins_len_ratio', 0.5)} --min-junction-fraction {m.get('min_junction_fraction', 0.75)} "
            f"--min-carriers {m.get('min_carriers', 1)} --window {m.get('hotspot_window', 1000)} "
            f"--out-all {S}/private/cxsv_panel_status.tsv --out-private {S}/private/private_cxsv_master_table.tsv "
            f"--out-bed {S}/private/private_cxsv_loci.bed --out-junctions {S}/private/junction_panel_status.tsv", log)
    slit.fill_callable(S, callable_at)
    cats = " ".join(f"--catalogue {n}:{v}:{t}:{mf}" for n, v, t, mf in cfg["catalogues"])
    cm_ = sc["external_catalogues"].get("match", {})
    slit.sh(f"python {W}/match_catalogues.py --junctions {S}/junctions/junctions.tsv --loci {S}/private/cxsv_panel_status.tsv "
            f"{cats} --tolerance {cm_.get('tolerance', 500)} --ins-len-ratio {cm_.get('ins_len_ratio', 0.5)} "
            f"--min-junction-fraction {cm_.get('min_junction_fraction', 0.75)} --out-loci {S}/private/cxsv_catalogue_hits.tsv "
            f"--out-junctions {S}/private/junction_catalogue_hits.tsv --out-check {S}/private/catalogue_check.tsv", log)
    cn = sc["copy_number"]
    slit.sh(f"python {W}/locus_segments.py --junctions {S}/junctions/junctions.tsv --loci {S}/private/cxsv_panel_status.tsv "
            f"--fai {a.fai} --collapse-distance {jp['collapse_distance']} --composite-fraction {cp['composite_fraction']} "
            f"--flank {cn['flank_bp']} --gap {cn['gap_bp']} --margin {cn['margin_bp']} --out {S}/scoring/locus_segments.bed", log)
    slit.synth_depth(f"{S}/scoring/locus_segments.bed", f"{S}/scoring/locus_depth.tsv", cn_rows, sex)
    with open(f"{S}/scoring/scoring_params.json", "w") as fh:
        json.dump(slit.scoring_params(cfg, meta, sex), fh, indent=2, sort_keys=True)
    for c in ("sniffles", "cutesv", "svim"):
        open(f"{S}/{c}_empty.vcf", "w").close()
    slit.sh(f"python {W}/score_cxsv.py --status {S}/private/cxsv_panel_status.tsv --catalogue-hits {S}/private/cxsv_catalogue_hits.tsv "
            f"--catalogue-check {S}/private/catalogue_check.tsv --junctions {S}/junctions/junctions.tsv "
            f"--junction-status {S}/private/junction_panel_status.tsv --vcf {S}/cxsv/cxsv_only.vcf "
            f"--caller-vcf Sniffles2={S}/sniffles_empty.vcf --caller-vcf cuteSV={S}/cutesv_empty.vcf "
            f"--caller-vcf SVIM={S}/svim_empty.vcf --depth {S}/scoring/locus_depth.tsv --params {S}/scoring/scoring_params.json "
            f"--genes {ann['genes_bed']} --exons {ann['exons_bed']} --repeats {ann['repeats_bed']} --segdups {ann['segdup_bed']} "
            f"--fai {a.fai} --sample {ev['event_id']} --out-tsv {S}/scoring/cxsv_scores.tsv "
            f"--out-features {S}/scoring/cxsv_score_features.tsv --out-report {S}/scoring/cxsv_score_report.md", log)
    return slit.read_tsv(f"{S}/scoring/cxsv_scores.tsv")


def run_event(task):
    ev, junctions, cn_rows = task
    slit, a, cfg, meta, callable_at = G["slit"], G["a"], G["cfg"], G["meta"], G["callable_at"]
    res = slit.run_event((ev, junctions, cn_rows, "main", a, cfg, meta, callable_at))
    S0 = os.path.join(a.work, "events", ev["event_id"], "main")
    S1 = os.path.join(a.work, "events", ev["event_id"], "cn")
    os.makedirs(S1, exist_ok=True)
    sex = ev["carrier_sex"] if ev["carrier_sex"] in ("male", "female") else "female"
    donor = G["donor"][sex]
    contigs = sorted({j["chrom1"] for j in junctions} | {j["chrom2"] for j in junctions} | {r["chrom"] for r in cn_rows})
    track = modelled_track(donor["fine"], cn_rows, sex, G["par"])
    segs, _ = cnlib.segment_genome(track, sex, G["par"], G["pon"], family=donor["family"],
                                   contigs=[c for c in contigs if c in cnlib.CONTIGS])
    recs = cnlink.read_records(f"{S0}/merged_sv.vcf")
    others = G["others"][donor["family"]]
    cands = cnlink.candidates(segs, recs, others)
    rep = augment.augment(S0, S1, f"{S0}/merged_sv.vcf", cands, a.repo, a.fai, f"{S1}/augment.log")
    n_cn = sum(r["written"] == "yes" for r in rep)
    scores = []
    if res["cxsv_loci"] or n_cn:
        open(f"{S1}/run.log", "w").close()
        scores = score_chain(S1, ev, cn_rows, sex, slit, a, cfg, meta, callable_at)
    cn_ids = {r["locus_id"] for r in rep if r["written"] == "yes"}
    return dict(event=ev, main=res, cn_scores=scores, cn_ids=cn_ids,
                cands=[dict(decision=c["decision"], cls=c["cls"], n_segments=len(c["segs"]), n_records=len(c["records"]),
                            n_long=c["n_long"], public=",".join(map(str, c["public_counts"]))) for c in cands],
                segments=[f"{s['chrom']}:{s['start']}-{s['end']}:{s['state']}/{s['ploidy']}" for s in segs])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--literature", required=True)
    ap.add_argument("--repo", required=True)
    ap.add_argument("--release", required=True)
    ap.add_argument("--resources", default="/mnt/sargasso1/CxSV/resources")
    ap.add_argument("--fai", default="/mnt/sargasso1/CxSV/refs/hs1.fa.fai")
    ap.add_argument("--seg", required=True, help="segment_all.py output (pon/, segments_all.tsv)")
    ap.add_argument("--depth-root", required=True)
    ap.add_argument("--donor-male", required=True, help="key of the male donor genome")
    ap.add_argument("--donor-female", required=True)
    ap.add_argument("--genomes", required=True, help="genomes.tsv (segment_all.py)")
    ap.add_argument("--families", required=True)
    ap.add_argument("--par", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--jobs", type=int, default=3)
    ap.add_argument("--events", default=None)
    ap.add_argument("--sources", default="nazaryan2018,collins2017,eisfeldt2024,bahrambeigi2019,grochowski2024")
    a = ap.parse_args()
    os.makedirs(a.work, exist_ok=True)
    slit = load_slit(a.repo)
    cfg = slit.load_config(a)
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
    genomes = {g["key"]: g for g in slit.read_tsv(a.genomes)}
    fam = {r["individual_id"]: r["family_id"] for r in slit.read_tsv(a.families)}
    par = cnlib.read_par(a.par)
    strata = {genomes[a.donor_male]["stratum"], genomes[a.donor_female]["stratum"]}
    assert len(strata) == 1, "donors must share a stratum"
    pon = cnlib.PanelOfNormals.from_npz(f"{a.seg}/pon/{strata.pop()}.npz")
    donor = {}
    for sex, key in (("male", a.donor_male), ("female", a.donor_female)):
        donor[sex] = dict(key=key, fine=cnlib.read_bins(f"{a.depth_root}/{key}/depth.regions.bed.gz"),
                          family=fam.get(genomes[key]["individual_id"], genomes[key]["individual_id"]))
    allseg = defaultdict(list)
    for r in slit.read_tsv(f"{a.seg}/segments_all.tsv"):
        if genomes.get(r["key"], {}).get("role") == "panel":
            allseg[(r["key"], r["family"])].append(dict(chrom=r["chrom"], start=int(r["start"]), end=int(r["end"]),
                                                        state=int(r["state"]), ploidy=int(r["ploidy"])))
    others = {d["family"]: [(k, f, s) for (k, f), s in allseg.items() if f != d["family"]] for d in donor.values()}
    G.update(slit=slit, a=a, cfg=cfg, meta=meta, callable_at=callable_at, pon=pon, donor=donor, par=par, others=others)
    tasks = [(e, by_j[e["event_id"]], by_cn[e["event_id"]]) for e in events if by_j[e["event_id"]]]
    with cf.ProcessPoolExecutor(a.jobs) as ex:
        results = list(ex.map(run_event, tasks))
    rows, summ = [], []
    for r in results:
        ev = r["event"]
        main_loci = r["main"]["loci"]
        best_main = max(main_loci, key=lambda L: float(L["Total_Score"]), default=None)
        cn_loci = r["cn_scores"]
        best_cn = max(cn_loci, key=lambda L: float(L["Total_Score"]), default=None)
        cn_only = [L for L in cn_loci if L["Locus_ID"] in r["cn_ids"]]
        best_cnl = max(cn_only, key=lambda L: float(L["Total_Score"]), default=None)
        for L in cn_loci:
            rows.append(dict(event_id=ev["event_id"], author_event_type=ev["author_event_type"], cn_locus=L["Locus_ID"] in r["cn_ids"],
                             **{k: L.get(k) for k in ("Locus_ID", "Cmplx_Class", "Panel_Status", "Genes_Affected", "Total_Score", "Tier", "Top_Features")}))
        summ.append(dict(event_id=ev["event_id"], author_event_type=ev["author_event_type"], carrier_sex=ev["carrier_sex"],
                         dup_trp=ev["author_event_type"] == "DUP-TRP/INV-DUP",
                         main_called=bool(main_loci), main_best_tier=best_main["Tier"] if best_main else "not_called",
                         main_best_score=best_main["Total_Score"] if best_main else "NA",
                         cn_called=bool(cn_loci), cn_best_tier=best_cn["Tier"] if best_cn else "not_called",
                         cn_best_score=best_cn["Total_Score"] if best_cn else "NA",
                         cn_loci=len(cn_only), cn_locus_classes=";".join(sorted({L["Cmplx_Class"] for L in cn_only})),
                         cn_locus_best_tier=best_cnl["Tier"] if best_cnl else "NA",
                         cn_locus_best_score=best_cnl["Total_Score"] if best_cnl else "NA",
                         candidates=";".join(f"{c['decision']}:{c['cls']}:{c['n_segments']}seg:{c['n_records']}rec:{c['n_long']}long"
                                             for c in r["cands"]),
                         segments=";".join(r["segments"])))
    slit.write(f"{a.work}/literature_cn_scores.tsv", rows)
    slit.write(f"{a.work}/literature_cn_summary.tsv", summ)


if __name__ == "__main__":
    main()
