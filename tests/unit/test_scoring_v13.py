"""panel-v1.3 scoring changes (roadmap item 16; docs/PLAN_panel_v1.3.md):

  (e) tier ceiling: a public locus is never reported above Tier 2 (PI decision D1)
  (h) type-agnostic caller agreement: callers agree on a junction when their breakends
      match within the matching tolerance, whatever SV type each gave
  (c) one sex source, the copy-number track's call; ambiguous-sex flag in the report
  (f) depth-noise gate: CN sigma above the per-stratum threshold -> conservative mode
      for copy-number-derived loci
"""
import json
import subprocess
import sys

import pytest
import yaml

import cxsv_junctions as cj
import cxsv_scoring as sc
from conftest import REPO

TIERS = [{"name": "Tier1_high_priority", "min": 0.99}, {"name": "Tier2", "min": 0.50},
         {"name": "Tier3", "min": 0.00}, {"name": "Benign_likely", "min": -1000}]


def _scorer():
    sys.path.insert(0, str(REPO / "workflow" / "scripts"))
    import score_cxsv  # noqa: E402
    return score_cxsv


def _J(c1, p1, c2, p2, svtype, supp):
    return cj.Junction(c1, p1, c2, p2, svtype, p2 - p1, supp, 0, ["r"])


# ─── (e) tier ceiling ────────────────────────────────────────────────────────
def test_tier_ceiling_caps_only_above_the_ceiling():
    cap = {"public": "Tier2"}
    assert sc.tier_ceiling("Tier1_high_priority", "public", cap, TIERS) == ("Tier2", "Tier2")
    assert sc.tier_ceiling("Tier2", "public", cap, TIERS) == ("Tier2", None)          # at the ceiling
    assert sc.tier_ceiling("Tier3", "public", cap, TIERS) == ("Tier3", None)
    assert sc.tier_ceiling("Tier1_high_priority", "private", cap, TIERS) == ("Tier1_high_priority", None)
    assert sc.tier_ceiling("Tier1_high_priority", "shared_junctions", cap, TIERS)[1] is None
    assert sc.tier_ceiling("Tier1_high_priority", "public", {}, TIERS)[1] is None      # no ceiling configured
    with pytest.raises(ValueError):
        sc.tier_ceiling("Tier1_high_priority", "public", {"public": "Tier9"}, TIERS)
    assert sc.tier_label("Tier1_high_priority") == "Tier 1" and sc.tier_label("Tier2") == "Tier 2"


def test_brsk2_like_public_locus_capped_with_reason():
    """BOC_2057 chr11:1,465,980-1,509,727 (BRSK2): public, 65/1,582 carriers, 1.06 points."""
    scx = _scorer()
    row = {"Panel_Status": "public", "Panel_Carriers": "65", "Panel_Individuals": "1582"}
    loci = [{"row": row, "total": 1.06, "score_tier": sc.tier_of(1.06, TIERS)},
            {"row": {**row, "Panel_Status": "private", "Panel_Carriers": "0"}, "total": 1.46,
             "score_tier": "Tier1_high_priority"}]
    scx.apply_tier_ceiling(loci, {"tiers": TIERS, "tier_ceiling_by_status": {"public": "Tier2"}})
    assert loci[0]["score_tier"] == "Tier1_high_priority" and loci[0]["tier"] == "Tier2"
    assert loci[0]["ceiling_note"].startswith("capped at Tier 2: seen in 65 healthy panel individuals")
    assert loci[0]["total"] == 1.06                                          # points unchanged
    assert loci[1]["tier"] == "Tier1_high_priority" and loci[1]["ceiling"] is None
    # params written before item 16e (no key) get the D1 default; an explicit {} turns it off
    again = [dict(loci[0])]
    scx.apply_tier_ceiling(again, {"tiers": TIERS})
    assert again[0]["tier"] == "Tier2"
    scx.apply_tier_ceiling(again, {"tiers": TIERS, "tier_ceiling_by_status": {}})
    assert again[0]["tier"] == "Tier1_high_priority"


def test_config_sets_public_ceiling_and_drops_flank_sex():
    cfg = yaml.safe_load((REPO / "config" / "scoring.yaml").read_text())
    assert cfg["tier_ceiling_by_status"] == {"public": "Tier2"}
    assert {t["name"] for t in cfg["tiers"]} >= {"Tier2"}
    assert "sex_ratio_male" not in cfg["copy_number"] and "sex_ratio_female" not in cfg["copy_number"]
    assert "infer_sex" not in (REPO / "workflow" / "scripts" / "score_cxsv.py").read_text()


# ─── (h) type-agnostic caller agreement ──────────────────────────────────────
def test_nf1_like_del_and_inv_from_two_callers_agree():
    """BOC_4433 NF1 chr17:31,902,782-32,067,453: Sniffles2 DEL + SVIM INV at the same
    breakends, no cuteSV (100 kb limit) -> was 0/2 (-0.10), now 2/2 (+0.10)."""
    w = {"multi_fraction": 0.5, "multi": 0.10, "none": -0.10}
    sniffles_del = _J("chr17", 31902782, "chr17", 32067453, "DEL", "001")      # SUPP_VEC cuteSV,SVIM,Sniffles2
    svim_inv = _J("chr17", 31902790, "chr17", 32067441, "INV", "010")
    kept = [sniffles_del, svim_inv]
    old = sc.f_caller_agreement([j.supp_vec for j in kept], w)                # type-aware (before 16h)
    assert old.points == -0.10 and old.value.startswith("0/2")
    callers = sc.junction_callers(kept, kept, 500)
    assert callers == [{1, 2}, {1, 2}]
    new = sc.f_caller_agreement(callers, w, 500)
    assert new.points == 0.10 and new.value == "2/2 junctions with >= 2 callers (breakends within 500 bp, any SV type)"


def test_caller_agreement_needs_both_breakends_within_tolerance():
    a = _J("chr1", 1000, "chr1", 50000, "DEL", "001")
    far = _J("chr1", 1000, "chr1", 52000, "INV", "010")          # second end 2 kb away
    swapped = _J("chr1", 50100, "chr1", 1100, "BND", "100")        # BND written the other way round
    other_chrom = _J("chr2", 1000, "chr2", 50000, "DEL", "010")
    assert sc.junction_callers([a], [a, far], 500) == [{2}]
    assert sc.junction_callers([a], [a, swapped], 500) == [{0, 2}]
    assert sc.junction_callers([a], [a, other_chrom], 500) == [{2}]
    # a composite record (not kept) still lends its caller
    assert sc.f_caller_agreement(sc.junction_callers([a], [a, swapped], 500), {"multi_fraction": 0.5, "multi": 0.1,
                                                                             "none": -0.1}).points == 0.1


def test_library_qc_gate_uses_type_agnostic_count():
    scx = _scorer()

    class J:
        supp_vec = "001"
    loci = [{"L": "chr17:1-2", "tier": "Tier1_high_priority", "kept": [J(), J()], "n_multi": 2, "segs": None,
             "row": {}}]
    P = {"tiers": TIERS, "library_qc": {"conservative_mode": True, "conservative": {
        "min_multi_caller_fraction": 0.5, "min_cn_changed": 1, "demote_to": "Tier2"}}}
    scx.apply_library_qc_gate(loci, {"qc_flag": "FLAG"}, P)
    assert loci[0]["tier"] == "Tier1_high_priority" and loci[0]["qc_gate"].startswith("kept: 2/2")


# ─── (c) sex from the copy-number track ──────────────────────────────────────
def test_resolve_sex_uses_cn_call_and_flags_ambiguous_and_discordant():
    st = {"sex": "female", "sex_call": "female", "x_ratio": 1.0126, "y_ratio": 0.0}
    r = sc.resolve_sex(st, None)
    assert (r["sex"], r["sex_call"], r["ambiguous"], r["discordant"]) == ("female", "female", False, False)
    assert "chrX ratio 1.013" in r["note"] and "copy-number track" in r["note"]
    r = sc.resolve_sex(st, "male")                                   # sheet says male: reported, not used
    assert r["sex"] == "female" and r["discordant"] and r["given"] == "male"
    amb = sc.resolve_sex({"sex": "female", "sex_call": "ambiguous", "x_ratio": 0.70, "y_ratio": 0.01}, "female")
    assert amb["ambiguous"] and amb["sex"] == "female" and not amb["discordant"]
    assert sc.resolve_sex(None, "male")["sex"] == "male"              # no CN track: the sheet
    none = sc.resolve_sex(None, None)
    assert none["sex"] == "unknown" and none["sex_call"] == "NA"


def test_sample_qc_section_prints_sex_call_ratios_and_flag():
    scx = _scorer()
    P = {"tiers": TIERS, "tier_ceiling_by_status": {"public": "Tier2"}}
    SX = sc.resolve_sex({"sex": "male", "sex_call": "ambiguous", "x_ratio": 0.62, "y_ratio": 0.09}, None)
    text = "\n".join(scx.sample_qc_section(SX, sc.depth_noise({}, {"enabled": False}), [], P))
    assert "sex_call ambiguous, x_ratio 0.62, y_ratio 0.09" in text
    assert sc.AMBIGUOUS_SEX_FLAG in text
    assert "public at most Tier 2; 0 locus/loci capped" in text


# ─── (f) depth-noise gate ────────────────────────────────────────────────────
def test_depth_noise_states():
    cfg = {"enabled": True, "sigma_max": {"ONT_LRSC_R10": 0.25, "ONT_LRSC_R9": None}, "default_sigma_max": None}
    assert sc.depth_noise({"sigma": 0.32, "stratum": "ONT_LRSC_R10"}, cfg)["state"] == "high"
    assert sc.depth_noise({"sigma": 0.21, "stratum": "ONT_LRSC_R10"}, cfg)["state"] == "normal"
    assert sc.depth_noise({"sigma": 0.32, "stratum": "ONT_LRSC_R9"}, cfg)["state"] == "not_calibrated"
    assert sc.depth_noise({"sigma": 0.32, "stratum": "X"}, cfg)["state"] == "not_calibrated"
    assert sc.depth_noise({"stratum": "ONT_LRSC_R10"}, cfg)["state"] == "not_assessed"
    assert sc.depth_noise({"sigma": 0.32}, {"enabled": False})["state"] == "off"
    hi = sc.depth_noise({"sigma": 0.32, "stratum": "ONT_LRSC_R10"}, cfg)
    assert hi["text"].startswith("Depth noise high: copy-number-derived calls unreliable")


def test_depth_noise_gate_demotes_only_cn_derived_tier1_without_junction_support():
    scx = _scorer()

    def x(cls, crit, tier, n_multi, n):
        return {"row": {"Cmplx_Class": cls, "CxSV_Criteria": crit}, "tier": tier, "kept": [None] * n,
                "n_multi": n_multi}
    loci = [x("CN_linked_loss", "CN_linked;C4", "Tier1_high_priority", 0, 2),      # CN-derived, no support
            x("DUP-TRP_INV-DUP", "C2", "Tier1_high_priority", 2, 2),             # CN class, junction support
            x("DEL+INV", "C2,C4", "Tier1_high_priority", 0, 2),                  # not CN-derived
            x("CN_linked_gain", "CN_linked;", "Tier2", 0, 1)]                    # CN-derived, below Tier 1
    cfg = yaml.safe_load((REPO / "config" / "scoring.yaml").read_text())["depth_noise"]
    P = {"tiers": TIERS, "depth_noise": cfg}
    hi = {"state": "high", "text": "Depth noise high: copy-number-derived calls unreliable"}
    scx.apply_depth_noise_gate(loci, hi, P)
    assert [l["tier"] for l in loci] == ["Tier2", "Tier1_high_priority", "Tier1_high_priority", "Tier2"]
    assert loci[0]["dn_gate"].startswith("demoted") and loci[1]["dn_gate"].startswith("kept")
    assert loci[2]["dn_gate"] == "not_applied" and loci[3]["dn_gate"].startswith("flagged")
    text = "\n".join(scx.sample_qc_section(sc.resolve_sex(None, None), hi, loci, P))
    assert "**Depth noise high: copy-number-derived calls unreliable.**" in text and "1 demoted, 1 kept" in text
    for l in loci:                                                     # normal noise: nothing happens
        l["tier"] = "Tier1_high_priority"
    scx.apply_depth_noise_gate(loci, {"state": "normal"}, P)
    assert all(l["dn_gate"] == "not_applied" and l["tier"] == "Tier1_high_priority" for l in loci)


def test_depth_noise_config_covers_every_pon_stratum_and_is_calibrated():
    sco = yaml.safe_load((REPO / "config" / "scoring.yaml").read_text())["depth_noise"]
    cal = yaml.safe_load((REPO / "config" / "calling.yaml").read_text())["cn_params"]
    strata = set(cal["default_stratum"].values()) | set(cal["chemistry_stratum"].values())
    assert strata <= set(sco["sigma_max"]), strata - set(sco["sigma_max"])
    # panel-v1.3 design check D3 (docs/design_check_v13): median + 3.5 scaled MAD of the panel per stratum
    assert sco["sigma_max"] == {"ONT_LRSC_R10": 0.233, "ONT_LRSC_R9": 0.220, "ONT_1kGP_1019": 0.401,
                                "PacBio_HiFi": 0.198}
    assert sco["default_sigma_max"] is None        # every PoN stratum is listed
    assert sco["enabled"] is True and sco["conservative"]["demote_to"] == "Tier2"


# ─── end to end: score_cxsv.py ───────────────────────────────────────────────
def _score_fixture(t, status="public", carriers="65", stats=None):
    """One DEL+INV locus over a HI-3 gene: Sniffles2 DEL and SVIM INV at the same breakends."""
    L = "chr1:100000-200000"
    def w(name, text):
        (t / name).write_text(text)
        return str(t / name)
    junctions = w("j.tsv", "\t".join(cj.JUNCTION_COLUMNS) + "\n"
                  + f"chr1\t100000\tchr1\t200000\tDEL\t100000\t{L}\t001\t-1\t1\tS.DEL.1\n"
                  + f"chr1\t100010\tchr1\t199990\tINV\t99980\t{L}\t010\t0\t1\tV.INV.1\n")
    cols = ["Chr", "Start", "End", "Size", "SV_Types", "Cmplx_Class", "CxSV_Criteria", "N_Junctions", "N_Breakends",
            "Max_Interleaved", "N_Inter_Junctions", "Complexity_Tier", "Panel_Carriers", "Panel_Callable",
            "Panel_Individuals", "Panel_Status"]
    st = w("status.tsv", "\t".join(cols) + "\n" + "\t".join(
        ["chr1", "100000", "200000", "100000", "DEL,INV", "DEL+INV", "C2_nested_SV,C4_CN_plus_orient", "2", "2", "1",
         "0", "simple", carriers, "1500", "1582", status]) + "\n")
    cats = w("cats.tsv", "Locus_ID\tCatalogues_Containing\tCatalogue_Max_AF\n" + f"{L}\tnone\t.\n")
    check = w("check.tsv", "catalogue\tversion\tstatus\n")
    vcf = w("cx.vcf", "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t100000\tS.DEL.1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=200000;SUPP_VEC=001\n"
            "chr1\t100010\tV.INV.1\tN\t<INV>\t.\tPASS\tSVTYPE=INV;END=199990;SUPP_VEC=010\n")
    depth = w("depth.tsv", "".join(f"chr1\t{s}\t{e}\t{L}\t{k}\t0\t{int(d * (e - s))}\n" for s, e, k, d in (
        (79000, 99000, "flank_left", 30), (100200, 199800, "segment", 15), (201000, 221000, "flank_right", 30))))
    js = w("js.tsv", "chrom1\tpos1\tchrom2\tpos2\tsvtype\tsvlen\tpanel_carriers\n")
    genes = w("genes.bed", "chr1\t150000\t250000\tHIG\n")
    exons = w("exons.bed", "chr1\t150000\t150500\tHIG\n")
    reps = w("rep.bed", "chr1\t900000\t900010\tSimple_repeat\n")
    segd = w("sd.bed", "chr1\t900000\t910000\n")
    fai = w("ref.fa.fai", "chr1\t1000000\t6\t60\t61\n")
    dos = w("dos.tsv", "gene_symbol\tgene_id\thi_score\tts_score\tlast_evaluated\nHIG\t1\t3\t0\t2020\n")
    cfg = yaml.safe_load((REPO / "config/scoring.yaml").read_text())
    P = {k: cfg[k] for k in ("weights", "tiers", "copy_number", "clustering", "foldback", "tier_ceiling_by_status",
                             "depth_noise")}
    # every locus Tier 1 by points (order of the tiers kept), so the ceiling is what decides
    P["tiers"] = [dict(t, min={"Tier1_high_priority": -10.0, "Tier2": -20.0, "Tier3": -30.0}.get(t["name"], t["min"]))
                  for t in P["tiers"]]
    P.update({"scoring_resources": {"clingen_dosage": dos}, "collapse_distance": 300, "composite_fraction": 0.25,
              "panel_match": {"tolerance": 500, "ins_len_ratio": 0.5, "min_junction_fraction": 0.75},
              "panel_release_id": "TEST", "patient": {"sex": "female", "hpo_terms": []}})
    P["clustering"]["n_simulations"] = 50
    pj = w("params.json", json.dumps(P))
    args = ["--status", st, "--catalogue-hits", cats, "--catalogue-check", check, "--junctions", junctions,
            "--junction-status", js, "--vcf", vcf, "--depth", depth, "--params", pj, "--genes", genes, "--exons", exons,
            "--repeats", reps, "--segdups", segd, "--fai", fai, "--sample", "S1",
            "--out-tsv", str(t / "s.tsv"), "--out-features", str(t / "f.tsv"), "--out-report", str(t / "r.md")]
    if stats is not None:
        args += ["--cn-stats", w("stats.json", json.dumps(stats))]
    r = subprocess.run([sys.executable, str(REPO / "workflow/scripts/score_cxsv.py"), *args],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    lines = (t / "s.tsv").read_text().splitlines()
    return dict(zip(lines[0].split("\t"), lines[1].split("\t"))), (t / "r.md").read_text()


def test_score_cli_caps_public_and_agrees_across_sv_types(tmp_path):
    stats = {"stratum": "ONT_LRSC_R10", "sex": "male", "sex_call": "ambiguous", "x_ratio": 0.66, "y_ratio": 0.08,
             "sigma": 0.2}
    row, md = _score_fixture(tmp_path, stats=stats)
    assert row["caller_agreement_points"] == "0.10" and row["caller_agreement_value"].startswith("2/2")
    assert row["Score_Tier"] == "Tier1_high_priority" and row["Tier"] == "Tier2"
    assert row["Tier_Ceiling"].startswith("capped at Tier 2: seen in 65 healthy panel individuals")
    assert "capped at Tier 2: seen in 65 healthy panel individuals" in md and "Tier2 (capped)" in md
    assert "sex_call ambiguous, x_ratio 0.66, y_ratio 0.08" in md and sc.AMBIGUOUS_SEX_FLAG in md
    assert "sample-sheet sex: female" in md                     # reported, not used for ploidy
    assert "sex male (copy-number track; call ambiguous" in row["ct_cn_oscillation_value"]
    # config/scoring.yaml depth_noise (design check D3): ONT_LRSC_R10 sigma_max 0.233; fixture sigma 0.2
    assert "depth noise normal (CN sigma 0.2 <= 0.233, stratum ONT_LRSC_R10)" in md
    priv, _ = _score_fixture(tmp_path, status="private", carriers="0", stats=stats)
    assert priv["Tier"] == "Tier1_high_priority" and priv["Tier_Ceiling"] == "not_applied"
