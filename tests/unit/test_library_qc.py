"""Library QC (docs/LIBRARY_QC.md): call-level counts, coverage-adjusted flag,
conservative Tier-1 gate, config and the calling fingerprint's independence.
The read-level part needs pysam and is tested only where pysam is installed."""
import math
import sys

import pytest
import yaml

from conftest import REPO

sys.path.insert(0, str(REPO / "workflow" / "scripts"))
sys.path.insert(0, str(REPO / "workflow" / "lib"))
import library_qc as lq  # noqa: E402
import cxsv_scoring as sc  # noqa: E402
import cxsv_fingerprint as fpm  # noqa: E402

HDR = "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\n"


def rec(svtype, supp):
    return f"chr1\t100\tx\tN\t<{svtype}>\t.\tPASS\tSUPP={supp.count('1')};SUPP_VEC={supp};SVTYPE={svtype}\tGT\t./.\t./.\t./.\n"


def test_call_metrics_counts(tmp_path):
    p = tmp_path / "m.vcf"
    p.write_text(HDR + rec("TRA", "010") + rec("BND", "010") + rec("TRA", "100") + rec("TRA", "111")
                 + rec("DEL", "010") + rec("DEL", "011") + rec("INV", "001"))
    m = lq.call_metrics(str(p))
    assert m == {"merged_records": 7, "bnd_records": 4, "single_caller_bnd": 3, "svim_only_bnd": 2,
                 "single_caller_any": 5, "svim_only_any": 3}


def test_zscore_and_flag():
    cal = {"svim_only_bnd_per_gb": {"a": -1.0, "b": 0.5, "r0": 0.0, "mad": 0.1}}
    row = {"svim_only_bnd_per_gb": 10 ** (-1.0 + 0.5 * math.log10(100) + 0.7) - 1e-3, "coverage_est": 100.0}
    z = lq.zscore(row["svim_only_bnd_per_gb"], row["coverage_est"], cal["svim_only_bnd_per_gb"], 1e-3)
    assert z == pytest.approx(7.0, abs=1e-6)
    assert lq.flag(row, cal, 6, 1e-3)[0] == "FLAG"
    assert lq.flag(row, cal, 8, 1e-3)[0] == "PASS"
    assert lq.flag(row, {}, 6, 1e-3)[0] == "NA"
    # no coverage estimate (no window on the reference): not assessable
    assert lq.flag({"svim_only_bnd_per_gb": 5.0, "coverage_est": float("nan")}, cal, 6, 1e-3)[0] == "NA"


def test_conservative_gate():
    w = {"min_multi_caller_fraction": 0.5, "min_cn_changed": 1}
    assert sc.library_qc_gate(2, 3, 0, w)[0]          # 2/3 > 0.5
    assert not sc.library_qc_gate(1, 2, 0, w)[0]      # 1/2 is not > 0.5, no CN change
    assert sc.library_qc_gate(1, 2, 1, w)[0]          # depth support
    assert sc.library_qc_gate(0, 2, 3, w)[0]
    assert not sc.library_qc_gate(0, 0, 0, w)[0]


def test_config_default_on_and_groups():
    cfg = yaml.safe_load((REPO / "config" / "library_qc.yaml").read_text())["library_qc"]
    assert cfg["conservative_mode"] is True      # PI decision 2026-11-05
    assert set(cfg["groups"]) == {"ONT_1kGP_1019", "ONT_LRSC_R10", "ONT_LRSC_R9", "PacBio_HiFi"}
    for g in cfg["groups"].values():
        assert set(g) == {"svim_only_bnd_per_gb", "svim_only_any_per_gb"}
        for c in g.values():
            assert c["mad"] > 0
    assert set(cfg["default_group"]) == {"ONT", "PACBIO"}
    assert all(v in cfg["groups"] for v in cfg["default_group"].values())


def test_windows_file():
    rows = [l.split("\t") for l in (REPO / "config" / "library_qc_windows.hs1.bed").read_text().splitlines()
            if l and not l.startswith("#")]
    assert len(rows) == 20 and all(int(e) - int(s) == 1_000_000 for _, s, e, _ in rows)


def test_not_in_calling_fingerprint():
    assert "library_qc" not in fpm.CALLING_PARAM_KEYS
    assert not any("library_qc" in f for f in fpm.CALLING_CODE)


def test_read_metrics_with_pysam(tmp_path):
    pysam = pytest.importorskip("pysam")
    hdr = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 10000}, {"SN": "chr2", "LN": 10000}]}
    path = str(tmp_path / "t.bam")
    with pysam.AlignmentFile(path, "wb", header=hdr) as out:
        for i, (pos, mapq, sa, flag) in enumerate([
                (900, 60, None, 0),                        # starts before the window: not counted
                (1000, 60, None, 0),                       # plain
                (1100, 60, "chr1,5000,+,50M,60,0;", 0),    # SA intra-chromosomal
                (1200, 60, "chr2,500,+,50M,60,0;", 0),     # SA to another chromosome
                (1300, 10, "chr2,500,+,50M,60,0;", 0),     # MAPQ < 20: not counted
                (1400, 60, None, 2048)]):                  # supplementary: not counted
            a = pysam.AlignedSegment()
            a.query_name, a.reference_id, a.reference_start = f"r{i}", 0, pos
            a.mapping_quality, a.flag = mapq, flag
            a.cigarstring, a.query_sequence = "100M", "A" * 100
            if sa:
                a.set_tag("SA", sa)
            out.write(a)
    pysam.index(path)
    m = lq.read_metrics(path, None, [("chr1", 1000, 2000), ("chrZ", 0, 10)], 20)
    assert m["windows_used"] == 1 and m["reads_sampled"] == 3
    assert m["frac_sa"] == pytest.approx(2 / 3) and m["frac_sa_interchrom"] == pytest.approx(1 / 3)
    assert m["coverage_est"] == pytest.approx(300 / 1000)


class _J:
    def __init__(self, supp_vec):
        self.supp_vec = supp_vec


def _loci():
    def x(L, tier, vecs, segs):
        return {"L": L, "tier": tier, "kept": [_J(v) for v in vecs], "segs": segs, "total": 1.5, "rank": 1,
                "row": {"Cmplx_Class": "BND_cluster"}}
    return [x("chr1:100-200", "Tier1_high_priority", ["010", "100"], None),          # single-caller, no depth
            x("chr2:100-200", "Tier1_high_priority", ["111", "011", "010"], None),   # 2/3 multi-caller
            x("chr3:100-200", "Tier2", ["010"], None)]                               # not Tier 1


def _params(mode=True):
    return {"tiers": [{"name": "Tier1_high_priority", "min": 0.99}, {"name": "Tier2", "min": 0.5},
                      {"name": "Tier3", "min": 0.0}],
            "library_qc": {"conservative_mode": mode,
                           "conservative": {"min_multi_caller_fraction": 0.5, "min_cn_changed": 1,
                                            "demote_to": "Tier2"}}}


def _scorer():
    sys.path.insert(0, str(REPO / "workflow" / "scripts"))
    import score_cxsv  # noqa: E402
    return score_cxsv


def test_gate_and_report_flagged():
    scx = _scorer()
    loci, P = _loci(), _params(True)
    scx.apply_library_qc_gate(loci, {"qc_flag": "FLAG"}, P)
    assert [x["tier"] for x in loci] == ["Tier2", "Tier1_high_priority", "Tier2"]
    assert loci[0]["qc_gate"].startswith("demoted") and loci[1]["qc_gate"].startswith("kept")
    assert loci[2]["qc_gate"] == "not_applied"
    summary, section = scx.library_qc_report({"qc_flag": "FLAG", "qc_group": "ONT_LRSC_R10"}, P, loci)
    text = "\n".join(summary + section)
    assert "Conservative mode DEMOTED 1 Tier-1 locus/loci" in text
    assert "loci demoted by conservative mode" in text and "chr1:100-200" in text
    assert "| Tier1_high_priority | Tier2 |" in text and "chr2:100-200" not in "\n".join(section)


def test_gate_not_applied_when_pass_or_off():
    scx = _scorer()
    for qc, mode in (({"qc_flag": "PASS"}, True), ({"qc_flag": "NA"}, True), ({"qc_flag": "FLAG"}, False)):
        loci, P = _loci(), _params(mode)
        scx.apply_library_qc_gate(loci, qc, P)
        assert [x["tier"] for x in loci] == ["Tier1_high_priority", "Tier1_high_priority", "Tier2"]
        summary, section = scx.library_qc_report(qc, P, loci)
        assert section == []
        if qc["qc_flag"] == "PASS":
            assert "did not act here" in summary[-1]
        if not mode:
            assert "Conservative mode is OFF" in summary[-1]
