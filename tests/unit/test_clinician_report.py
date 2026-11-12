"""Clinician report of the patient workflow (roadmap item 16l; scripts/clinician_report.py,
lib/cxsv_report_figure.py, rules/clinician_report.smk). Synthetic rows only."""
import json
import sys
import zipfile

import pytest
import yaml

import cxsv_scoring as sc
from conftest import REPO

sys.path.insert(0, str(REPO / "workflow" / "scripts"))
import clinician_report as cr  # noqa: E402

TIERS = [{"name": "Tier1_high_priority", "min": 0.99}, {"name": "Tier2", "min": 0.50},
         {"name": "Tier3", "min": 0.00}, {"name": "Benign_likely", "min": -1000}]


def row(locus, tier, total, status="private", carriers="0", cls="DEL+INV", genes="NF1:broken,RNF135:lost", **kw):
    c, se = locus.split(":")
    s, e = se.split("-")
    r = {"Locus_ID": locus, "Chr": c, "Start": s, "End": e, "Size": str(int(e) - int(s)), "Cmplx_Class": cls,
         "SV_Types": "DEL,INV", "N_Junctions": "2", "Panel_Status": status, "Panel_Carriers": carriers,
         "Panel_Callable": "1361", "Panel_Individuals": "1582", "Catalogues_Containing": "none",
         "Catalogue_Max_AF": ".", "Genes_Affected": genes, "Total_Score": f"{total:.2f}", "Tier": tier,
         "Library_QC_Gate": "not_applied", "Tier_Ceiling": "not_applied", "Depth_Noise_Gate": "not_applied",
         "CN_Derived": "no", "hpo_match_value": "no patient HPO terms", "omim_value": "off (no licence key)",
         "inheritance_value": "parents not sequenced",
         "gene_dosage_value": "NF1:broken(HI 3,TS 0,LOEUF 0.449); RNF135:lost(HI 0,TS 0,LOEUF 1.49)",
         "region_dosage_value": "none", "catalogue_frequency_value": "HPRC 0/2"}
    for f in sc.FEATURE_ORDER:
        r.setdefault(f + "_points", "0.00")
        r.setdefault(f + "_value", "NA")
    r.update({"gene_dosage_points": "0.90", "panel_privacy_points": "0.26", "exons_points": "0.15",
              "complexity_class_points": "0.10", "caller_agreement_points": "0.10",
              "caller_agreement_value": "2/2 junctions with >= 2 callers (breakends within 500 bp, any SV type)"})
    r.update(kw)
    return r


def data(rows, cn=None, sigma_max=None, libqc=None):
    T = yaml.safe_load((REPO / "config" / "clinician_report.yaml").read_text())["clinician_report"]
    cfg = yaml.safe_load((REPO / "config" / "scoring.yaml").read_text())
    dn = dict(cfg["depth_noise"])
    if sigma_max is not None:
        dn["sigma_max"] = {"ONT_LRSC_R10": sigma_max}
    P = {"tiers": TIERS, "panel_release_id": "cxsv-panel-TEST", "patient": {"hpo_terms": []}, "depth_noise": dn,
         "tier_ceiling_by_status": cfg["tier_ceiling_by_status"]}
    return {"sample": "S1", "T": T, "P": P, "scores": rows,
            "cn": cn if cn is not None else {"sex": "female", "sex_call": "female", "x_ratio": 1.0126, "y_ratio": 0.0,
                                             "sigma": 0.2066, "stratum": "ONT_LRSC_R10"},
            "cn_only": [], "libqc": libqc or {"qc_flag": "PASS", "qc_group": "ONT_LRSC_R10", "qc_group_source": "header",
                                              "median_read_len": "3385", "reads_sampled": "74808"},
            "coverage": 25.4, "fp": ("fdac9ff80401" + "0" * 52, "COMPATIBLE"), "bundle": "cxsv-bundle-TEST",
            "commit": "84d380e08cc4", "date": "2026-11-11"}


def text_of(blocks):
    out = []
    for b in blocks:
        for x in b[1:]:
            if isinstance(x, str):
                out.append(x)
            elif isinstance(x, list):
                out += [str(c) for r in x for c in (r if isinstance(r, list) else [r])]
    return "\n".join(out)


def test_one_tier1_report_has_every_section_and_no_comparator():
    rows = [row("chr17:31902782-32067453", "Tier1_high_priority", 1.66),
            row("chr11:1465980-1509727", "Tier2", 1.06, status="public", carriers="65", cls="Chromothripsis_like",
                genes="BRSK2:broken", Tier_Ceiling="capped at Tier 2: seen in 65 healthy panel individuals "
                                                  "(panel status public; Tier 1 by points)"),
            row("chr1:10-20000", "Tier3", 0.2), row("chr2:10-20000", "Benign_likely", -1.2)]
    fig = {"path": "x.png", "n_total": 1, "n_shown": 1, "width_px": 10, "height_px": 4}
    B = cr.sample_blocks(data(rows), fig)
    kinds = [b[0] for b in B]
    heads = [b[1] for b in B if b[0] == "h"]
    assert heads == ["1. Sample quality", "2. Findings", "3. What the tiers mean", "4. Limitations"]
    assert kinds[:3] == ["title", "kv", "box"] and "figure" in kinds
    box = next(b for b in B if b[0] == "box")
    assert box[3] == "alert" and box[2].startswith("1 Tier 1 finding requiring expert review")
    t = text_of(B)
    assert "cxsv-bundle-TEST" in t and "cxsv-panel-TEST" in t and "fdac9ff80401 (COMPATIBLE with panel)" in t
    assert "NF1 disrupted: ClinGen HI 3, TS 0, LOEUF 0.449" in t
    assert "gene dosage +0.90" in t and "caller agreement +0.10" in t       # every non-zero feature
    assert "sex_call female, x_ratio 1.013, y_ratio 0" in t
    # config/scoring.yaml depth_noise (design check D3): ONT_LRSC_R10 sigma_max 0.233; fixture sigma 0.2066
    assert "normal (CN sigma 0.207, threshold 0.233 for ONT_LRSC_R10)" in t
    assert "capped at Tier 2" in t and "Tier 3: 1 loci; likely benign: 1 loci (of 4 CxSV loci scored)" in t
    assert "\u2265 0.99" in t and "0.50 to < 0.99" in t                      # thresholds from the params
    assert "No phenotype (HPO) terms, OMIM annotation or parental samples were used" in t
    for word in ("LOFO", "per genome", "cohort", "healthy-panel leave"):   # no per-cohort comparator text
        assert word not in t
    assert "diagnosis" in t and "does not make a diagnosis" in t


def test_zero_tier1_prints_a_note_and_no_figure():
    B = cr.sample_blocks(data([row("chr1:10-20000", "Tier3", 0.2)]), None)
    assert "figure" not in [b[0] for b in B]
    box = next(b for b in B if b[0] == "box")
    assert box[2] == "No Tier 1 complex structural variant detected." and box[3] == "neutral"
    assert "Tier 1: none. (No Tier 1 locus, so no Tier 1 points graphic.)" in text_of(B)


def test_ambiguous_sex_and_noisy_depth_are_flagged():
    cn = {"sex": "male", "sex_call": "ambiguous", "x_ratio": 0.66, "y_ratio": 0.08, "sigma": 0.32,
          "stratum": "ONT_LRSC_R10"}
    rows = [row("chr17:0-2", "Tier2", 0.6, cls="CN_linked_loss", CN_Derived="yes")]
    B = cr.sample_blocks(data(rows, cn=cn, sigma_max=0.25), None)
    t = text_of(B)
    assert sc.AMBIGUOUS_SEX_FLAG in t and "AMBIGUOUS (ambiguous)" in t
    assert "Depth noise high: copy-number-derived calls unreliable" in t
    assert "malformed locus coordinates" in t
    box = next(b for b in B if b[0] == "box")
    assert box[3] == "caution" and "sex-chromosome check ambiguous" in box[2] and "noisy depth" in box[2]


def test_figure_caps_rows_and_skips_zero(tmp_path):
    pytest.importorskip("matplotlib")
    pytest.importorskip("PIL")
    import cxsv_report_figure as rf
    rows = [row(f"chr1:{i * 100000}-{i * 100000 + 50000}", "Tier1_high_priority", 1.0 + i / 100) for i in range(1, 21)]
    f = rf.tier1_points_figure(rows, TIERS, str(tmp_path / "t1.png"), max_rows=15)
    assert f["n_total"] == 20 and f["n_shown"] == 15 and (tmp_path / "t1.png").stat().st_size > 1000
    shown, _ = rf.select_rows(rows, 15)
    assert float(shown[0]["Total_Score"]) == 1.20                           # highest first
    assert rf.tier1_points_figure([], TIERS, str(tmp_path / "none.png")) is None
    cat = rf.category_points(rows[0])
    assert cat["gene/region dosage"] == pytest.approx(0.90) and cat["privacy"] == pytest.approx(0.26)
    assert set(cat) == set(rf.GORDER) and set(rf.FGROUP) == set(sc.FEATURE_ORDER)


def test_renders_docx_and_pdf(tmp_path):
    pytest.importorskip("docx")
    pytest.importorskip("reportlab")
    pytest.importorskip("matplotlib")
    import cxsv_report_figure as rf
    rows = [row("chr17:31902782-32067453", "Tier1_high_priority", 1.66), row("chr1:10-20000", "Tier3", 0.2)]
    fig = rf.tier1_points_figure(rows[:1], TIERS, str(tmp_path / "f.png"))
    B = cr.sample_blocks(data(rows), fig)
    cr.to_docx(B, str(tmp_path / "r.docx"), "Research report")
    cr.to_pdf(B, str(tmp_path / "r.pdf"), "Research report")
    assert (tmp_path / "r.pdf").read_bytes()[:5] == b"%PDF-"
    z = zipfile.ZipFile(tmp_path / "r.docx")
    assert any(n.startswith("word/media/") for n in z.namelist())
    assert "Sex-chromosome check" in z.read("word/document.xml").decode()


def test_rule_is_wired_into_the_patient_workflow():
    snake = (REPO / "workflow" / "patient" / "Snakefile").read_text()
    assert 'include: "../rules/clinician_report.smk"' in snake and '"config", "clinician_report.yaml"' in snake
    assert "clinician_report.docx" in snake and "clinician_report.pdf" in snake
    rule = (REPO / "workflow" / "rules" / "clinician_report.smk").read_text()
    assert "{{sample}}.clinician_report.docx" in rule and "{{sample}}.tier1_points.png" in rule
    T = yaml.safe_load((REPO / "config" / "clinician_report.yaml").read_text())["clinician_report"]
    assert set(T["tier_meaning"]) == {t["name"] for t in TIERS}
    assert T["footer"].startswith("Research report") and "accredited laboratory" in T["footer"]
    # the scripts the rule runs are not in the calling fingerprint
    import cxsv_fingerprint as fpm
    assert not any("report" in f for f in fpm.CALLING_CODE)
    json.dumps(T)                                                           # serialisable for the text rule
