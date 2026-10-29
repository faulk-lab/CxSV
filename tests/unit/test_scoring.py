"""Scoring features (workflow/lib/cxsv_scoring.py), catalogue matching and the scorer CLI."""
import gzip
import json
import subprocess
import sys

import yaml

import cxsv_junctions as cj
import cxsv_matching as cm
import cxsv_scoring as sc
from conftest import REPO

W = yaml.safe_load((REPO / "config/scoring.yaml").read_text())["weights"]


# ─── privacy ─────────────────────────────────────────────────────────────────
def test_panel_privacy_is_weighted_by_callability():
    full = sc.f_panel_privacy("private", "0", "10", "10", W["panel_privacy"])
    half = sc.f_panel_privacy("private", "0", "10", "5", W["panel_privacy"])
    unknown = sc.f_panel_privacy("private", "0", "10", "NA", W["panel_privacy"])
    none = sc.f_panel_privacy("private", "0", "10", "0", W["panel_privacy"])
    assert full.points == 0.30 and half.points == 0.15 and unknown.points == 0.15 and none.points == 0.0
    assert sc.f_panel_privacy("public", "3", "10", "10", W["panel_privacy"]).points == -1.0
    assert sc.f_panel_privacy("shared_junctions", "0", "10", "10", W["panel_privacy"]).points == 0.0


def test_catalogue_most_severe_counts_never_summed():
    w = W["catalogue_frequency"]
    cats = [("A", "3/3", 0.05, True), ("B", "3/3", 0.002, True), ("C", "1/3", 0.4, False)]
    f = sc.f_catalogue_frequency(cats, 0.30, w)
    assert f.points == -1.30                       # common in A: -1, and the private bonus is cancelled
    f = sc.f_catalogue_frequency(cats, -1.0, w)    # already public in the panel: privacy = min(-1, -1)
    assert f.points == 0.0
    f = sc.f_catalogue_frequency([("B", "3/3", 0.002, True)], 0.30, w)
    assert f.points == -0.80
    f = sc.f_catalogue_frequency([("C", "1/3", 0.4, False)], 0.30, w)   # partial match: not contained
    assert f.points == 0.0 and "C 1/3" in f.value


# ─── genes ───────────────────────────────────────────────────────────────────
def test_gene_dosage_uses_clingen_then_loeuf():
    dosage = {"HIG": {"hi": "3", "ts": "0"}, "TSG": {"hi": "0", "ts": "3"}, "AR": {"hi": "30", "ts": None},
              "UNL": {"hi": "40", "ts": None}}
    loeuf = {"AR": 0.2, "NEW": 0.5, "UNL": 0.1}
    w = W["gene_dosage"]
    assert sc.f_gene_dosage({"HIG": "broken"}, dosage, loeuf, w).points == 0.90
    assert sc.f_gene_dosage({"HIG": "lost"}, dosage, loeuf, w).points == 1.00
    assert sc.f_gene_dosage({"TSG": "gained"}, dosage, loeuf, w).points == 1.00
    assert sc.f_gene_dosage({"TSG": "broken"}, dosage, loeuf, w).points == 0.0
    assert sc.f_gene_dosage({"AR": "broken"}, dosage, loeuf, w).points == 0.30      # LOEUF < 0.35
    assert sc.f_gene_dosage({"NEW": "broken"}, dosage, loeuf, w).points == 0.15     # LOEUF < 0.6
    assert sc.f_gene_dosage({"UNL": "broken"}, dosage, loeuf, w).points == 0.0      # curated 40 wins
    f = sc.f_gene_dosage({"HIG": "broken", "NEW": "broken"}, dosage, loeuf, w)
    assert f.points == 0.90 and f.rule.startswith("HIG")                           # max, not sum


def test_gene_count_exons_and_intergenic():
    assert sc.f_gene_count(3, W["gene_count"]).points == 0.0
    assert sc.f_gene_count(30, W["gene_count"]).points == 0.45
    assert sc.f_gene_count(40, W["gene_count"]).points == 0.90
    assert sc.f_exons(2, W["exons"]).points == 0.15 and sc.f_exons(0, W["exons"]).points == 0.0
    assert sc.f_no_genes(0, 0, W["no_genes"]).points == -0.60
    assert sc.f_no_genes(1, 0, W["no_genes"]).points == 0.0


def test_class_size_callers_reads_artefacts():
    assert sc.f_class("Chromothripsis_like", W["complexity_class"]).points == 0.30
    assert sc.f_class("Complex_other", W["complexity_class"]).points == 0.05
    assert sc.f_size(50_000, W["size"]).points == 0.0 and sc.f_size(2_000_000, W["size"]).points == 0.15
    assert sc.f_caller_agreement(["101", "011"], W["caller_agreement"]).points == 0.10
    assert sc.f_caller_agreement(["001", "001"], W["caller_agreement"]).points == -0.10
    assert sc.f_caller_agreement(["001", "001", "111"], W["caller_agreement"]).points == 0.0
    assert sc.f_read_support(3, W["read_support"]).points == -0.15
    assert sc.f_read_support(12, W["read_support"]).points == 0.0
    assert sc.f_read_support(None, W["read_support"]).points == 0.0
    assert sc.f_artefact(1.0, 0.0, W["artefact_context"]).points == -0.30
    assert sc.f_artefact(0.5, 0.5, W["artefact_context"]).points == -0.225


# ─── chromothripsis criteria ─────────────────────────────────────────────────
def test_interleaving_and_tier_points():
    assert sc.f_interleaving(2, W["ct_interleaving"]).points == 0.0
    assert sc.f_interleaving(3, W["ct_interleaving"]).points == 0.20
    assert sc.f_interleaving(7, W["ct_interleaving"]).points == 0.30
    assert sc.f_breakpoint_tier("extreme", 12, W["ct_breakpoint_tier"]).points == 0.25
    assert sc.f_breakpoint_tier("intermediate", 4, W["ct_breakpoint_tier"]).points == 0.05


def test_record_orientation_from_alt_brackets():
    assert sc.record_orientation("chr1", 100, "A[chr1:500[", "BND") == {("+", "-")}       # DEL-like
    assert sc.record_orientation("chr1", 100, "A]chr1:500]", "INV") == {("+", "+")}
    assert sc.record_orientation("chr1", 100, "[chr1:500[A", "INV") == {("-", "-")}
    assert sc.record_orientation("chr1", 500, "]chr1:100]A", "BND") == {("+", "-")}       # seen from the far end
    assert sc.record_orientation("chr1", 100, "<INV>", "INV") == {("+", "+"), ("-", "-")}
    assert sc.record_orientation("chr1", 100, "<DUP>", "DUP") == {("-", "+")}
    assert sc.record_orientation("chr1", 100, "ACGT", "INS") == set()


def test_join_balance_and_foldbacks():
    counts, p = sc.join_balance([{("+", "-")}, {("-", "+")}, {("+", "+")}, {("-", "-")}])
    assert counts == {"DEL-like": 1, "DUP-like": 1, "h2hINV": 1, "t2tINV": 1} and p > 0.99
    counts, p = sc.join_balance([{("+", "+")}] * 6)
    assert p < 0.01
    J = lambda a, b: cj.Junction("chr1", a, "chr1", b, "INV", b - a, "001", 1, [f"{a}"])
    js = [J(1000, 7000), J(20000, 30000), J(40000, 45000), J(40001, 45002)]
    orients = [{("-", "-")}, {("+", "+")}, {("+", "+")}, {("-", "-")}]
    fb = sc.foldbacks(js, orients, 300, 30000, [("chr1", 2000, 6000)])
    # the last two are each other's reciprocal (a balanced inversion): not fold-backs
    assert [(j.pos1, g) for j, g in fb] == [(1000, True), (20000, False)]
    w = W["ct_bfb_modifier"]
    f = sc.f_bfb(0.40, counts, p, 2, 1, 2, False, {"min_fraction": 0.5}, w)
    assert f.points == -0.20
    f = sc.f_bfb(0.40, counts, p, 2, 1, 2, True, {"min_fraction": 0.5}, w)
    assert f.points == -0.30
    assert sc.f_bfb(0.40, counts, p, 2, 0, 2, True, {"min_fraction": 0.5}, w).points == 0.0   # no CN gain


def test_ploidy_sex_and_par():
    par = [("chrX", 0, 2394410), ("chrX", 153925834, 154259566)]
    assert sc.ploidy("chr5", 10, "male", par) == 2
    assert sc.ploidy("chrX", 152_000_000, "male", par) == 1
    assert sc.ploidy("chrX", 1_000_000, "male", par) == 2        # PAR1
    assert sc.ploidy("chrX", 152_000_000, "female", par) == 2
    assert sc.ploidy("chrY", 10_000_000, "female", par) == 0
    fl = [("chr1", 0, 10, 30.0), ("chr2", 0, 10, 31.0), ("chrX", 100_000_000, 100_000_010, 15.0)]
    assert sc.infer_sex(fl, par)[0] == "male"
    fl[-1] = ("chrX", 100_000_000, 100_000_010, 30.0)
    assert sc.infer_sex(fl, par)[0] == "female"


def test_cn_segments_and_oscillation():
    p = {"min_baseline_depth": 5, "min_segment_bp": 1000, "min_cn_delta": 0.6}
    rows = [("chr1", 0, 20000, "flank_left", 30.0), ("chr1", 21000, 25000, "segment", 15.0),
            ("chr1", 25500, 30000, "segment", 30.0), ("chr1", 30500, 34000, "segment", 14.0),
            ("chr1", 34500, 39000, "segment", 31.0), ("chr1", 40000, 60000, "flank_right", 30.0)]
    segs = sc.cn_segments(rows, "female", [], p)
    assert [s[4] for s in segs] == [1, 2, 1, 2] and [s[6] for s in segs] == [True, False, True, False]
    assert sc.oscillation(segs) == (4, 2)
    f = sc.f_cn_oscillation(segs, W["ct_cn_oscillation"], "sex female (given)")
    assert f.points == 0.30
    # male chrX: flank depth 15 is one copy; a segment at 30 is a gain to 2 copies
    rows = [("chrX", 0, 20000, "flank_left", 15.0), ("chrX", 21000, 25000, "segment", 30.0),
            ("chrX", 26000, 46000, "flank_right", 15.0)]
    segs = sc.cn_segments([(c, s + 10_000_000, e + 10_000_000, k, d) for c, s, e, k, d in rows], "male", [], p)
    assert segs[0][5] == 1 and segs[0][4] == 2 and segs[0][6]
    assert sc.cn_segments([("chr1", 0, 100, "segment", 30.0)], "female", [], p) is None


def test_clustering_test_and_bh():
    eff = sc.EffectiveCoords(100_000_000, [(10_000_000, 20_000_000)])
    assert eff.length == 90_000_000 and eff.masked(15_000_000) and not eff.masked(5)
    assert eff.eff(25_000_000) == 15_000_000 and eff.eff(5) == 5
    locus = [50_000_000, 50_010_000, 50_020_000, 50_030_000]
    res = sc.clustering_test(locus, locus + [1_000, 80_000_000], eff, 200, 1)
    assert res["k"] == 4 and res["n_chrom"] == 6 and res["p_scan"] < 0.01
    again = sc.clustering_test(locus, locus + [1_000, 80_000_000], eff, 200, 1)
    assert again == res                                                     # seeded: deterministic
    q = sc.bh_qvalues([0.01, None, 0.04, 0.5])
    assert q[1] is None and abs(q[0] - 0.03) < 1e-9 and abs(q[2] - 0.06) < 1e-9 and q[3] == 0.5
    f = sc.f_clustering(res, 0.1, W["ct_clustering"])
    assert f.points == 0.10 and sc.f_clustering(res, 0.3, W["ct_clustering"]).points == 0.0


def test_hpo_match_with_ancestors(tmp_path):
    obo = tmp_path / "hp.obo"
    obo.write_text("[Term]\nid: HP:0000001\n\n[Term]\nid: HP:0000118\nis_a: HP:0000001\n\n"
                   "[Term]\nid: HP:0001250\nis_a: HP:0000118\n\n[Term]\nid: HP:0002373\nis_a: HP:0001250\n\n"
                   "[Term]\nid: HP:0001263\nis_a: HP:0000118\n")
    parents = sc.read_obo_parents(str(obo))
    # gene annotated with febrile seizure (child of seizure) matches a patient with seizure
    assert sc.hpo_match(["HP:0001250"], {"HP:0002373"}, parents) == {"HP:0001250"}
    # the shared root alone never matches
    assert sc.hpo_match(["HP:0001263"], {"HP:0002373"}, parents) == set()
    f = sc.f_hpo(["HP:0001250", "HP:0001263"], {"HP:0001250"}, W["hpo_match"])
    assert f.points == 0.45
    assert sc.f_hpo([], set(), W["hpo_match"]).points == 0.0


def test_inheritance_and_omim_and_tiers():
    w = W["inheritance"]
    assert sc.f_inheritance({}, False, w).points == 0.0
    assert sc.f_inheritance({"mother": False, "father": False}, False, w).points == 0.45
    assert sc.f_inheritance({"mother": True, "father": False}, False, w).points == -0.30
    assert sc.f_inheritance({"mother": False}, False, w).points == 0.0
    assert sc.f_omim({"G"}, False, W["omim"]).points == 0.0
    assert sc.f_omim({"G"}, True, W["omim"]).points == 0.15
    tiers = yaml.safe_load((REPO / "config/scoring.yaml").read_text())["tiers"]
    assert sc.tier_of(1.2, tiers) == "Tier1_high_priority" and sc.tier_of(0.99, tiers) == "Tier1_high_priority"
    assert sc.tier_of(0.6, tiers) == "Tier2" and sc.tier_of(0.1, tiers) == "Tier3"
    assert sc.tier_of(-0.4, tiers) == "Benign_likely"


# ─── catalogue matching uses the panel matcher ───────────────────────────────
def test_catalogue_point_sites_and_query(tmp_path):
    p = tmp_path / "cat.tsv.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("chrom1\tpos1\tchrom2\tpos2\tsvtype\tsvlen\tac\tan\taf\tsource\tversion\taf_AFR\n"
                 "chr1\t1000\tchr1\t5000\tDEL\t4000\t2\t100\t0.02\tX\t1\t0.05\n"
                 "chr1\t90000\tchr1\t95000\tINV\t5000\t.\t.\t.\tX\t1\t.\n")
    sites = cm.read_catalogue_sites(str(p), {"chr1": [1200]}, 500)
    assert len(sites) == 1 and sites[0]["af"] == 0.02 and sites[0]["group_af"] == {"AFR": 0.05}
    idx = cm.SiteIndex(cm.read_catalogue_sites(str(p)))
    j = {"svtype": "DEL", "chrom1": "chr1", "pos1": 1300, "chrom2": "chr1", "pos2": 4800, "svlen": 3500}
    assert len(idx.query(j, 500, 0.5)) == 1
    assert idx.query({**j, "svtype": "DUP"}, 500, 0.5) == []
    assert idx.query({**j, "pos1": 1600}, 500, 0.5) == []
    k, n, _ = cm.match_fraction([j, {**j, "pos1": 90000, "pos2": 95000, "svtype": "INV"}], idx, 500, 0.5)
    assert (k, n) == (2, 2)


# ─── scripts end to end on a tiny locus ──────────────────────────────────────
def _write(path, text):
    path.write_text(text)
    return str(path)


def test_match_catalogues_and_score_cli(tmp_path):
    t = tmp_path
    L = "chr1:100000-130000"
    junctions = _write(t / "j.tsv", "\t".join(cj.JUNCTION_COLUMNS) + "\n"
                       + f"chr1\t100000\tchr1\t110000\tINV\t10000\t{L}\t101\t1\t2\tr1,r2\n"
                       + f"chr1\t104000\tchr1\t120000\tINV\t16000\t{L}\t011\t1\t1\tr3\n"
                       + f"chr1\t125000\tchr1\t130000\tDEL\t5000\t{L}\t111\t1\t1\tr4\n"
                       + "chr1\t500000\tchr1\t501000\tDEL\t1000\t.\t111\t1\t1\tr5\n")
    master = ["Chr", "Start", "End", "Size", "SV_Types", "Cmplx_Class", "CxSV_Criteria", "N_Junctions",
              "N_Breakends", "Max_Interleaved", "N_Inter_Junctions", "Complexity_Tier", "Panel_Carriers",
              "Panel_Callable", "Panel_Individuals", "Panel_Status"]
    status = _write(t / "status.tsv", "\t".join(master) + "\n" + "\t".join(
        ["chr1", "100000", "130000", "30000", "DEL,INV", "DEL+INV", "C1_min_breakpoints,C4_CN_plus_orient",
         "3", "6", "2", "0", "extreme", "0", "8", "10", "private"]) + "\n")
    cat = t / "C" / "1"
    cat.mkdir(parents=True)
    table = cat / "C.1.junction_sites.tsv.gz"
    with gzip.open(table, "wt") as fh:
        fh.write("chrom1\tpos1\tchrom2\tpos2\tsvtype\tsvlen\tac\tan\taf\tsource\tversion\n"
                 "chr1\t125100\tchr1\t130050\tDEL\t4950\t5\t100\t0.05\tC\t1\n")
    import hashlib
    sha = hashlib.sha256(table.read_bytes()).hexdigest()
    manifest = _write(cat / "C.1.manifest.tsv", f"field\tvalue\nversion\t1\noutput_sha256\t{sha}\n")
    out = {k: str(t / f"{k}.tsv") for k in ("loci", "jhits", "check")}
    scripts = REPO / "workflow/scripts"
    run = lambda *args: subprocess.run([sys.executable, *map(str, args)], capture_output=True, text=True)
    r = run(scripts / "match_catalogues.py", "--junctions", junctions, "--loci", status,
            "--catalogue", f"C:1:{table}:{manifest}", "--out-loci", out["loci"],
            "--out-junctions", out["jhits"], "--out-check", out["check"])
    assert r.returncode == 0, r.stderr
    row = dict(zip(*[l.split("\t") for l in open(out["loci"]).read().splitlines()]))
    assert row["C_Junctions_Matched"] == "1/3" and row["C_Max_AF"] == "0.05" and row["C_Contains"] == "no"
    assert row["Catalogues_Containing"] == "none"
    # a tampered table is refused
    bad = _write(cat / "bad.manifest.tsv", "field\tvalue\nversion\t1\noutput_sha256\t00\n")
    r = run(scripts / "match_catalogues.py", "--junctions", junctions, "--loci", status,
            "--catalogue", f"C:1:{table}:{bad}", "--out-loci", out["loci"],
            "--out-junctions", out["jhits"], "--out-check", out["check"])
    assert r.returncode != 0 and "do not match their manifests" in r.stderr

    r = run(scripts / "match_catalogues.py", "--junctions", junctions, "--loci", status,
            "--catalogue", f"C:1:{table}:{manifest}", "--out-loci", out["loci"],
            "--out-junctions", out["jhits"], "--out-check", out["check"])
    vcf = _write(t / "cx.vcf", "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                 "chr1\t100000\tr1\tN\tN]chr1:110000]\t.\tPASS\tSVTYPE=INV;END=110000;SUPP_VEC=101\n"
                 "chr1\t104000\tr3\tN\t[chr1:120000[N\t.\tPASS\tSVTYPE=INV;END=120000;SUPP_VEC=011\n"
                 "chr1\t125000\tr4\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=130000;SUPP_VEC=111\n")
    depth = _write(t / "depth.tsv", "".join(f"chr1\t{s}\t{e}\t{L}\t{k}\t0\t{int(d * (e - s))}\n" for s, e, k, d in (
        (79000, 99000, "flank_left", 30), (100200, 103800, "segment", 15), (104200, 109800, "segment", 30),
        (110200, 119800, "segment", 30), (120200, 124800, "segment", 30), (125200, 129800, "segment", 15),
        (131000, 151000, "flank_right", 30))))
    jstatus = _write(t / "js.tsv", "chrom1\tpos1\tchrom2\tpos2\tsvtype\tsvlen\tpanel_carriers\n"
                     "chr1\t100000\tchr1\t110000\tINV\t10000\t0\nchr1\t500000\tchr1\t501000\tDEL\t1000\t4\n")
    genes = _write(t / "genes.bed", "chr1\t101000\t140000\tHIG\nchr1\t600000\t700000\tFAR\n")
    exons = _write(t / "exons.bed", "chr1\t105000\t105500\tHIG\n")
    reps = _write(t / "rep.bed", "chr1\t99990\t100010\tSimple_repeat\n")
    segd = _write(t / "sd.bed", "chr1\t900000\t910000\n")
    fai = _write(t / "ref.fa.fai", "chr1\t1000000\t6\t60\t61\n")
    dos = _write(t / "dos.tsv", "gene_symbol\tgene_id\thi_score\tts_score\tlast_evaluated\nHIG\t1\t3\t0\t2020\n")
    con = _write(t / "con.tsv", "gene_symbol\ttranscript\tloeuf\tpli\ttranscript_choice\nHIG\tT\t0.1\t1\tmane_select\n")
    cfg = yaml.safe_load((REPO / "config/scoring.yaml").read_text())
    params = {k: cfg[k] for k in ("weights", "tiers", "copy_number", "clustering", "foldback")}
    params.update({"scoring_resources": {"clingen_dosage": dos, "constraint": con},
                   "panel_match": {"tolerance": 500, "ins_len_ratio": 0.5, "min_junction_fraction": 0.75},
                   "collapse_distance": 300, "composite_fraction": 0.25, "panel_release_id": "TEST",
                   "patient": {"sex": "female", "hpo_terms": []}})
    params["clustering"]["n_simulations"] = 100
    pj = _write(t / "params.json", json.dumps(params))
    res = {k: str(t / f"score.{k}") for k in ("tsv", "features", "md")}
    r = run(scripts / "score_cxsv.py", "--status", status, "--catalogue-hits", out["loci"], "--catalogue-check",
            out["check"], "--junctions", junctions, "--junction-status", jstatus, "--vcf", vcf, "--depth", depth,
            "--params", pj, "--genes", genes, "--exons", exons, "--repeats", reps, "--segdups", segd, "--fai", fai,
            "--out-tsv", res["tsv"], "--out-features", res["features"], "--out-report", res["md"])
    assert r.returncode == 0, r.stderr
    lines = open(res["tsv"]).read().splitlines()
    row = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
    assert row["Locus_ID"] == L and row["Rank"] == "1"
    assert row["gene_dosage_points"] == "0.90" and row["Genes_Affected"] == "HIG:broken"
    assert row["panel_privacy_points"] == "0.24"                  # 0.30 x 8/10 callable
    assert row["ct_cn_oscillation_points"] == "0.10"              # two losses, no 3-segment alternation
    total = sum(float(row[f"{n}_points"]) for n in sc.FEATURE_ORDER)
    assert abs(total - float(row["Total_Score"])) < 0.011
    assert row["Tier"] == "Tier1_high_priority"
    feats = open(res["features"]).read().splitlines()
    assert len(feats) == 1 + len(sc.FEATURE_ORDER)
    md = open(res["md"]).read()
    assert "not a validated cutoff" in md and "| gene_dosage |" in md and L in md
