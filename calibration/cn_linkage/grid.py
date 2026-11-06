#!/usr/bin/env python3
"""Test B: read-level simulation grid of DUP-TRP/INV-DUP (roadmap item 15, Phase 1).

Structure: the positive-control E01 template (tests/positive_control/events.yaml):
p B C D q -> p B+ C+ D+ C- B+ C+ D+ q (B, D duplicated, C triplicated, middle copy
inverted; junctions D+|C- and C-|B+). C (the triplicated segment, size T) is centred on
the locus gene; B and D are (S - T) / 2 each (S = total span).

Grid (one simulation per cell, coverages by down-sampling):
  loci       MECP2 (chrX, male host), PLP1 (chrX, male host), LMNB1 (chr5, female host)
  platforms  R10 (LRSC R10), R9 (LRSC R9), PacBio (HiFi): one male and one female host each
  S x T      S 0.2 / 0.5 / 1 / 2 Mb x T 50 / 100 / 200 / 500 kb with T <= S / 2
  coverage   10, 20, 30x (genome-median equivalent)

  grid.py design   --out GRID                   -> GRID/cells.tsv, GRID/cells/<cell>/events.yaml
  grid.py simulate --out GRID --cell CELL ...   (freya, env with badread) build_haplotypes.py,
                   badread at the host's local depth x 30 / host genome median (male chrX:
                   derivative alone; female: derivative + reference haplotype, half each),
                   minimap2 like the panel, keep reads overlapping the replace windows
  grid.py run      --out GRID --cell CELL --cov X ...   (any node, pipeline image):
                   host reads of the window (event +/- 1.5 Mb) minus reads overlapping the
                   replace windows, down-sampled to X / host median; simulated reads
                   down-sampled to X / 30; merged BAM; the patient workflow's discovery
                   (snakemake, targets cxsv + junctions); mosdepth 1 kb; segmentation of the
                   window against the host stratum's PoN (host family left out); cnlink;
                   augment; patient_rescore.py (panel-v1.1 release)
  grid.py collect  --out GRID                   -> GRID/grid_results.tsv
"""
import argparse
import csv
import gzip
import json
import os
import random
import shutil
import subprocess
import sys
import time

import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

LOCI = {"MECP2": ("chrX", 152333314, "male"), "PLP1": ("chrX", 102231523, "male"), "LMNB1": ("chr5", 127327059, "female")}
HOSTS = {("R10", "male"): "1kGP_LRSC_ONT_CHM13/HG03745", ("R10", "female"): "1kGP_LRSC_ONT_CHM13/HG02772",
         ("R9", "male"): "1kGP_LRSC_ONT_CHM13/NA18865", ("R9", "female"): "1kGP_LRSC_ONT_CHM13/HG03861",
         ("PacBio", "male"): "1kGP_PacBio_CHM13/HG00477", ("PacBio", "female"): "1kGP_PacBio_CHM13/HG00555"}
STRATUM = {"R10": "LRSC_R10", "R9": "LRSC_R9", "PacBio": "PacBio"}
SPANS = [200_000, 500_000, 1_000_000, 2_000_000]
TRIPS = [50_000, 100_000, 200_000, 500_000]
COVS = [10, 20, 30]
WINDOW = 1_500_000
FLANK_SIM, FLANK_REPLACE = 150_000, 50_000


def sh(cmd, log=None):
    if log:
        with open(log, "a") as fh:
            fh.write("$ " + cmd + "\n")
            fh.flush()
            subprocess.run(cmd, shell=True, check=True, stdout=fh, stderr=fh, executable="/bin/bash")
    else:
        subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")


def sh_out(cmd):
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True, executable="/bin/bash").stdout


def cells():
    out = []
    for locus in LOCI:
        for plat in ("R10", "R9", "PacBio"):
            for S in SPANS:
                for T in TRIPS:
                    if T <= S // 2:
                        out.append(dict(cell=f"{locus}_{plat}_S{S // 1000}_T{T // 1000}", locus=locus, platform=plat,
                                        span=S, trip=T))
    return out


def event_yaml(c):
    chrom, center, sex = LOCI[c["locus"]]
    T, F = c["trip"], (c["span"] - c["trip"]) // 2
    c0, c1 = center - T // 2, center + T // 2
    seg = {"p": [chrom, c0 - F - FLANK_SIM, c0 - F], "B": [chrom, c0 - F, c0], "C": [chrom, c0, c1],
           "D": [chrom, c1, c1 + F], "q": [chrom, c1 + F, c1 + F + FLANK_SIM]}
    ev = {"id": "G", "name": "DUP-TRP/INV-DUP grid", "gene": c["locus"], "segments": seg,
          "paths": {"der": ["p+", "B+", "C+", "D+", "C-", "B+", "C+", "D+", "q+"]},
          "zygosity": "hemizygous" if sex == "male" else "heterozygous"}
    return {"version": 1, "reference": "T2T-CHM13v2 (hs1)", "host": "grid", "seed": 15,
            "flank_sim": FLANK_SIM, "flank_replace": FLANK_REPLACE, "events": [ev]}


def design(a):
    os.makedirs(f"{a.out}/cells", exist_ok=True)
    rows = cells()
    with open(f"{a.out}/cells.tsv", "w") as fh:
        fh.write("cell\tlocus\tplatform\tspan\ttrip\tchrom\tB_start\tC_start\tC_end\tD_end\thost\tsex\n")
        for c in rows:
            y = event_yaml(c)
            d = f"{a.out}/cells/{c['cell']}"
            os.makedirs(d, exist_ok=True)
            yaml.safe_dump(y, open(f"{d}/events.yaml", "w"), sort_keys=False)
            s = y["events"][0]["segments"]
            chrom, _, sex = LOCI[c["locus"]]
            fh.write(f"{c['cell']}\t{c['locus']}\t{c['platform']}\t{c['span']}\t{c['trip']}\t{chrom}\t{s['B'][1]}\t"
                     f"{s['C'][1]}\t{s['C'][2]}\t{s['D'][2]}\t{HOSTS[(c['platform'], sex)]}\t{sex}\n")
    print(f"{len(rows)} cells x {len(COVS)} coverages = {len(rows) * len(COVS)} runs", file=sys.stderr)


def read_cells(out):
    with open(f"{out}/cells.tsv") as fh:
        return {r["cell"]: r for r in csv.DictReader(fh, delimiter="\t")}


def host_info(a, host):
    """alignment path, genome median depth and 1 kb track of the host (depth pass)."""
    sheet = {}
    cohort, sample = host.split("/")
    with open(f"{a.sheets}/{cohort}.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            sheet[r["sample_id"]] = r["alignment"]
    import cnlib
    fine = cnlib.read_bins(f"{a.depth_root}/{host}/depth.regions.bed.gz")
    gmed = cnlib.genome_median(cnlib.to_coarse(fine))
    return sheet[sample], gmed, fine


def local_depth(fine, windows):
    vals = []
    for c, s, e in windows:
        vals.append(fine[c][s // 1000:e // 1000])
    v = np.concatenate(vals)
    return float(np.median(v[v > 0]))


def simulate(a):
    import simulate_and_splice as ss       # tests/positive_control (on sys.path via --pc-code)
    c = read_cells(a.out)[a.cell]
    d = f"{a.out}/cells/{a.cell}"
    log = f"{d}/simulate.log"
    if os.path.exists(f"{d}/sim.ok"):
        return
    open(log, "w").close()
    t0 = time.time()
    sh(f"python {a.pc_code}/build_haplotypes.py {d}/events.yaml {a.ref} {d}/build --check", log)
    aln, gmed, fine = host_info(a, c["host"])
    win = ss.read_bed(f"{d}/build/G/replace.bed")
    loc = local_depth(fine, win)
    scale = 30.0 / gmed
    model = a.models[c["platform"]]
    stats = ss.read_stats(f"{model}/host_read_stats.tsv")
    fq = f"{d}/sim.fq"
    for p in (fq, f"{d}/origins.tsv"):
        if os.path.exists(p):
            os.remove(p)
    with open(f"{d}/origins.tsv", "w") as orig:
        if c["sex"] == "male":
            n_alt = ss.simulate(f"{d}/build/G/alt.fa", loc * scale, model, stats, 15_000 + hash(a.cell) % 1000, "G", "alt", fq, orig)
            n_ref = 0
        else:
            n_alt = ss.simulate(f"{d}/build/G/alt.fa", loc * scale / 2, model, stats, 15_000 + hash(a.cell) % 1000, "G", "alt", fq, orig)
            n_ref = ss.simulate(f"{d}/build/G/ref.fa", loc * scale / 2, model, stats, 25_000 + hash(a.cell) % 1000, "G", "ref", fq, orig)
    preset = "map-hifi" if c["platform"] == "PacBio" else "map-ont"
    index = a.index_hifi if c["platform"] == "PacBio" else a.index_ont
    bed = f"{d}/build/G/replace.bed"
    sh(f"minimap2 -t {a.threads} -x {preset} -a --eqx -K 4g -R '@RG\\tID:sim\\tSM:sim' {index} {fq} 2> {d}/minimap2.log "
       f"| samtools sort -@ 4 -o {d}/sim_all.bam - && samtools index {d}/sim_all.bam", log)
    sh(f"samtools view -L {bed} {d}/sim_all.bam | cut -f1 | sort -u > {d}/sim_keep.txt", log)
    sh(f"samtools view -b -N {d}/sim_keep.txt -o {d}/sim.bam {d}/sim_all.bam && samtools index {d}/sim.bam && "
       f"rm -f {d}/sim_all.bam {d}/sim_all.bam.bai {fq}", log)
    json.dump(dict(cell=a.cell, host=c["host"], host_aln=aln, host_genome_median=gmed, local_depth=loc,
                   sim_depth_alt=loc * scale if c["sex"] == "male" else loc * scale / 2, alt_reads=n_alt, ref_reads=n_ref,
                   kept_reads=int(sh_out(f"wc -l < {d}/sim_keep.txt").strip()), seconds=round(time.time() - t0, 1)),
              open(f"{d}/sim.json", "w"), indent=1)
    open(f"{d}/sim.ok", "w").close()


def run(a):
    import augment
    import cnlib
    import cnlink
    c = read_cells(a.out)[a.cell]
    d = f"{a.out}/cells/{a.cell}"
    R_out = f"{d}/cov{a.cov}"
    if os.path.exists(f"{R_out}/run.ok"):
        return
    # work on local disk: snakemake's mtime checks trip over CIFS timestamp skew
    R = f"{a.local}/{a.cell}/cov{a.cov}" if a.local else R_out
    shutil.rmtree(R, ignore_errors=True)
    os.makedirs(R)
    log = f"{R}/run.log"
    t0 = time.time()
    sim = json.load(open(f"{d}/sim.json"))
    aln, gmed = sim["host_aln"], sim["host_genome_median"]
    chrom = c["chrom"]
    lo, hi = max(0, int(c["B_start"]) - WINDOW), int(c["D_end"]) + WINDOW
    region = f"{chrom}:{lo + 1}-{hi}"
    bed = f"{d}/build/G/replace.bed"
    f_host = min(1.0, a.cov / gmed)
    f_sim = min(1.0, a.cov / 30.0)
    seed = 15 + a.cov
    X = a.exec_prefix
    sh(f"{X} samtools view -b --reference {a.ref} -o {R}/host_win.bam {aln} {region} && {X} samtools index {R}/host_win.bam", log)
    sh(f"{X} samtools view -L {bed} {R}/host_win.bam | cut -f1 | sort -u > {R}/host_drop.txt", log)
    sh(f"{X} samtools view -b -N ^{R}/host_drop.txt --subsample {f_host:.4f} --subsample-seed {seed} -o {R}/host_keep.bam {R}/host_win.bam", log)
    sh(f"{X} samtools view -b --subsample {f_sim:.4f} --subsample-seed {seed} -o {R}/sim_keep.bam {d}/sim.bam", log)
    # the simulated BAM's header lists contigs in the index order: merge, then sort
    sh(f"{X} samtools merge -f -o {R}/window_unsorted.bam {R}/host_keep.bam {R}/sim_keep.bam && "
       f"{X} samtools sort -@ 2 -o {R}/window.bam {R}/window_unsorted.bam && {X} samtools index {R}/window.bam && "
       f"rm -f {R}/window_unsorted.bam && "
       f"rm -f {R}/host_win.bam {R}/host_win.bam.bai {R}/host_keep.bam {R}/sim_keep.bam", log)
    # discovery: the patient workflow's rules up to junctions (callers, filters, merge, discovery)
    W = f"{R}/wf"
    os.makedirs(f"{W}/config", exist_ok=True)
    plat = "PacBio" if c["platform"] == "PacBio" else "ONT"
    with open(f"{W}/sheet.tsv", "w") as fh:
        fh.write(f"sample_id\talignment\tindividual_id\tplatform\nG\t{R}/window.bam\tG\t{plat}\n")
    with open(f"{W}/config/run.yaml", "w") as fh:
        # results under the work directory (the bundle site config points at /work/...)
        fh.write(f"patient_runs:\n  grid: {{samples: \"{W}/sheet.tsv\"}}\nthreads: {{samtools: 2, sv_callers: {a.threads}}}\n"
                 f"patient_results_root: results/patient\n")
    rel = "results/patient/samples/grid/G"          # snakemake targets relative to --directory
    S = f"{W}/{rel}"
    # the bundle's site config with /bundle replaced by the bundle's host path (no container mounts)
    with open(f"{a.bundle}/config/bundle_site.yaml") as fi, open(f"{W}/config/site.yaml", "w") as fo:
        fo.write(fi.read().replace("/bundle/", f"{a.bundle}/"))
    sh(f"{X} snakemake -s {a.code}/workflow/patient/Snakefile --directory {W} --cores {a.threads} "
       f"{rel}/junctions/junctions.tsv {rel}/cxsv/cxsv_master_table.tsv "
       f"--configfile {W}/config/site.yaml {W}/config/run.yaml > {R}/snakemake.log 2>&1", log)
    # depth of the window
    sh(f"{X} mosdepth -t 2 -n -x --mapq 20 --by 1000 -c {chrom} -f {a.ref} {R}/depth {R}/window.bam", log)
    fine = cnlib.read_bins(f"{R}/depth.regions.bed.gz", contigs=[chrom])
    par = cnlib.read_par(a.par)
    pon = cnlib.PanelOfNormals.from_npz(f"{a.seg}/pon/{STRATUM[c['platform']]}.npz")
    fam = a.host_family[c["host"]]
    segs = segment_window(fine[chrom], chrom, lo, hi, gmed * f_host, c["sex"], par, pon, fam)
    recs = cnlink.read_records(f"{S}/merged_sv.vcf")
    others = a.others(fam)
    cands = cnlink.candidates(segs, recs, others)
    A = f"{R}/aug"
    os.makedirs(f"{A}/callers", exist_ok=True)
    for f in ("sniffles_filtered.vcf", "cutesv_filtered.vcf", "svim_filtered.vcf"):
        shutil.copy(f"{S}/callers/{f}", f"{A}/callers/{f}")
    rep = augment.augment(S, A, f"{S}/merged_sv.vcf", cands, a.code, f"{a.ref}.fai", f"{R}/augment.log")
    # scoring parameters as the workflow writes them (patient: sex given, no HPO)
    params = make_params(a, c["sex"])
    os.makedirs(f"{S}/scoring", exist_ok=True)
    json.dump(params, open(f"{S}/scoring/scoring_params.json", "w"), indent=2, sort_keys=True)
    sh(f"python {HERE}/patient_rescore.py --run {S} --aug {A} --aln {R}/window.bam --release {a.release} "
       f"--code {a.code} --ref {a.ref} --sample G > {R}/rescore.json", log)
    res = outcome(c, A, rep, cands, segs)
    res.update(cell=a.cell, cov=a.cov, wall_s=round(time.time() - t0, 1), host_fraction=round(f_host, 4))
    json.dump(res, open(f"{R}/result.json", "w"), indent=1)
    for f in ("window.bam", "window.bam.bai"):
        if not a.keep_bam and os.path.exists(f"{R}/{f}"):
            os.remove(f"{R}/{f}")
    if R != R_out:
        # keep the small outputs with the cell; the local work tree is removed
        os.makedirs(f"{R_out}/aug", exist_ok=True)
        for f in ("result.json", "rescore.json", "run.log", "snakemake.log", "augment.log"):
            if os.path.exists(f"{R}/{f}"):
                shutil.copy(f"{R}/{f}", f"{R_out}/{f}")
        for f in ("scoring/cxsv_scores.tsv", "cxsv/clusters.bed", "cn_candidates.tsv"):
            if os.path.exists(f"{A}/{f}"):
                shutil.copy(f"{A}/{f}", f"{R_out}/aug/{os.path.basename(f)}")
        if os.path.exists(f"{S}/merged_sv.vcf"):
            shutil.copy(f"{S}/merged_sv.vcf", f"{R_out}/merged_sv.vcf")
        if not a.keep_bam:
            shutil.rmtree(R, ignore_errors=True)
    open(f"{R_out}/run.ok", "w").close()


def segment_window(track, chrom, lo, hi, gmed, sex, par, pon, fam):
    """Segments inside [lo, hi) of one contig; noise from the window's own bins."""
    import cnlib
    P = cnlib.DEFAULTS
    n = len(track)
    a0, a1 = lo // 1000, min(n, -(-hi // 1000))
    fine = np.full(n, np.nan, np.float32)
    fine[a0:a1] = track[a0:a1]
    coarse = cnlib.to_coarse({chrom: fine})[chrom]
    med, rel, cnt = pon.reference(chrom, sex, fam)
    pl = cnlib.ploidy_track(chrom, len(coarse), sex, par)
    with np.errstate(all="ignore"):
        r = coarse / gmed / med / (np.maximum(pl, 1) / 2.0)
    ok = np.isfinite(r) & (med >= P["min_pon_median"]) & (rel <= P["max_pon_relmad"])
    v = r[ok]
    sigma = float(1.4826 * np.median(np.abs(v - np.median(v)))) if len(v) else 0.1
    segs, _ = cnlib.segment_contig(chrom, coarse, fine, gmed, pl, med, rel, cnt, sigma, P)
    return segs


def make_params(a, sex):
    spec = __import__("importlib.util").util.spec_from_file_location("slit", f"{a.code}/tests/literature/score_literature_set.py")
    slit = __import__("importlib.util").util.module_from_spec(spec)
    spec.loader.exec_module(slit)

    class A:
        repo = a.code
        resources = "/mnt/sargasso1/CxSV/resources"
    cfg = slit.load_config(A)
    meta = json.load(open(f"{a.release}/release.json"))
    return slit.scoring_params(cfg, meta, sex)


def outcome(c, A, rep, cands, segs):
    """Called? (any CxSV locus overlapping the event); CN locus overlapping the event; tiers."""
    chrom, b0, d1 = c["chrom"], int(c["B_start"]), int(c["D_end"])
    rows = []
    with open(f"{A}/scoring/cxsv_scores.tsv") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    cn_ids = {r["locus_id"] for r in rep if r["written"] == "yes"}
    hit = [r for r in rows if r["Chr"] == chrom and int(r["End"]) >= b0 - 20000 and int(r["Start"]) <= d1 + 20000]
    std = [r for r in hit if r["Locus_ID"] not in cn_ids]
    cnl = [r for r in hit if r["Locus_ID"] in cn_ids]
    best = lambda L: max(L, key=lambda r: float(r["Total_Score"])) if L else None
    bs, bc = best(std), best(cnl)
    evseg = [s for s in segs if s["chrom"] == chrom and s["end"] >= b0 - 50000 and s["start"] <= d1 + 50000]
    off_t1 = [r for r in rows if r["Tier"].startswith("Tier1") and r not in hit]
    cand_ev = [k for k in cands if any(s["chrom"] == chrom and s["end"] >= b0 and s["start"] <= d1 for s in k["segs"])]
    return dict(std_called=bool(std), std_tier=bs["Tier"] if bs else "not_called", std_score=bs["Total_Score"] if bs else "NA",
                std_class=bs["Cmplx_Class"] if bs else "NA",
                cn_called=bool(cnl), cn_tier=bc["Tier"] if bc else "not_called", cn_score=bc["Total_Score"] if bc else "NA",
                cn_class=bc["Cmplx_Class"] if bc else "NA", cn_status=bc["Panel_Status"] if bc else "NA",
                event_segments=";".join(f"{s['start']}-{s['end']}:{s['state']}/{s['ploidy']}" for s in evseg),
                candidate_decisions=";".join(f"{k['decision']}:{k['cls']}:{len(k['records'])}rec" for k in cand_ev),
                any_called=bool(hit), best_tier=best(hit)["Tier"] if hit else "not_called",
                off_event_tier1=len(off_t1))


def collect(a):
    out = []
    for cell, c in read_cells(a.out).items():
        for cov in COVS:
            p = f"{a.out}/cells/{cell}/cov{cov}/result.json"
            if os.path.exists(p):
                r = json.load(open(p))
                r.update({k: c[k] for k in ("locus", "platform", "span", "trip", "sex")})
                out.append(r)
    cols = ["cell", "locus", "platform", "sex", "span", "trip", "cov", "std_called", "std_tier", "std_score", "std_class",
            "cn_called", "cn_tier", "cn_score", "cn_class", "cn_status", "any_called", "best_tier", "event_segments",
            "candidate_decisions", "off_event_tier1", "host_fraction", "wall_s"]
    with open(f"{a.out}/grid_results.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(out)
    print(f"{len(out)} runs collected", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("design", "simulate", "run", "collect"))
    ap.add_argument("--out", required=True)
    ap.add_argument("--cell")
    ap.add_argument("--cov", type=int)
    ap.add_argument("--ref", default="/mnt/sargasso1/CxSV/refs/hs1.fa")
    ap.add_argument("--index-ont", default="/mnt/sargasso1/CxSV/refs/hs1.map-ont.k17.mmi")
    ap.add_argument("--index-hifi", default="/mnt/sargasso1/CxSV/calibration/cn_linkage/refs/hs1.map-hifi.mmi")
    ap.add_argument("--models", default="")
    ap.add_argument("--pc-code", default="")
    ap.add_argument("--code", default="")
    ap.add_argument("--bundle", default="/mnt/sargasso1/CxSV/bundles/cxsv-bundle-20261105-8fc5bef725bd")
    ap.add_argument("--release", default="/mnt/sargasso1/CxSV/build/panel-v1.1-r2/results/panel/release")
    ap.add_argument("--sheets", default="/mnt/sargasso1/CxSV/datasets/registry/sample_sheets")
    ap.add_argument("--depth-root", default="/mnt/sargasso1/CxSV/calibration/cn_linkage/depth")
    ap.add_argument("--seg", default="")
    ap.add_argument("--par", default="/mnt/sargasso1/CxSV/resources/scoring/par.bed")
    ap.add_argument("--families", default="/mnt/sargasso1/CxSV/calibration/item10/families.tsv")
    ap.add_argument("--genomes", default="")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--exec-prefix", default="")
    ap.add_argument("--keep-bam", action="store_true")
    ap.add_argument("--local", default="", help="local work directory for run (default: in the cell directory)")
    a = ap.parse_args()
    if a.mode == "design":
        return design(a)
    if a.mode == "collect":
        return collect(a)
    if a.mode == "simulate":
        sys.path.insert(0, a.pc_code)
        a.models = dict(kv.split("=", 1) for kv in a.models.split(","))
        return simulate(a)
    fam = {}
    with open(a.families) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            fam[r["individual_id"]] = r["family_id"]
    a.host_family = {k: fam.get(k.split("/")[1], k.split("/")[1]) for k in HOSTS.values()}
    segs = {}
    with open(f"{a.seg}/segments_all.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["family"].startswith("patient:"):
                continue
            segs.setdefault((r["key"], r["family"]), []).append(dict(chrom=r["chrom"], start=int(r["start"]), end=int(r["end"]),
                                                                    state=int(r["state"]), ploidy=int(r["ploidy"])))
    a.others = lambda f: [(k, ff, s) for (k, ff), s in segs.items() if ff != f]
    run(a)


if __name__ == "__main__":
    main()
