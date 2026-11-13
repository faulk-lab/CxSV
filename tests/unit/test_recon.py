"""Structural reconstruction of Tier 1 loci (lib/cxsv_recon.py, lib/cxsv_recon_figure.py,
clinician report section). Synthetic inputs only: hand-made junctions and depth, and inputs
derived from the positive-control truth design (tests/positive_control/events.yaml), never
reads or calls from a human genome."""
import sys
from collections import Counter

import pytest
import yaml

from conftest import REPO

sys.path.insert(0, str(REPO / "workflow" / "scripts"))
sys.path.insert(0, str(REPO / "tests" / "positive_control"))
import cxsv_recon as R          # noqa: E402
import recon_vs_truth as V      # noqa: E402

PARAMS = {"collapse_distance": 300, "composite_fraction": 0.25,
          "copy_number": {"min_baseline_depth": 5, "min_cn_delta": 0.6, "min_segment_bp": 1000, "margin_bp": 200}}
BASE = 30.0


class NoGenes:
    def overlapping(self, c, s, e):
        return []


def jrow(c1, p1, c2, p2, svtype, rid, locus, supp="001"):
    return {"chrom1": c1, "pos1": p1, "chrom2": c2, "pos2": p2, "svtype": svtype, "svlen": abs(p2 - p1),
            "locus_id": locus, "supp_vec": supp, "dosage": 1, "n_records": 1, "record_ids": rid, "loci": {locus}}


def rec(caller, rid, svtype, c1, p1, c2, p2, orients, precise=True, reads=20):
    return R.Record(caller, rid, svtype, c1, p1, c2, p2, list(orients), reads, precise)


def depth(chrom, a, b, cn_by_interval):
    """Locus-depth rows: 20 kb flanks at CN 2 and the given (start, end, CN) segments."""
    rows = [(chrom, a - 21000, a - 1000, "flank_left", BASE), (chrom, b + 1000, b + 21000, "flank_right", BASE)]
    rows += [(chrom, s + 200, e - 200, "segment", BASE * cn / 2) for s, e, cn in cn_by_interval]
    return rows


def recon(locus, rows, recs, depth_rows, cn_rows=()):
    lrow = {"Locus_ID": locus, "gene_dosage_value": ""}
    return R.reconstruct(lrow, rows, recs, {}, depth_rows, list(cn_rows), NoGenes(), PARAMS, "female")


# ─── record orientation ─────────────────────────────────────────────────────
def test_record_orientation_from_alt_and_symbolic_types():
    bnd = R.record_ends("chr17", 31902782, "N[chr17:32067453[", {"SVTYPE": "BND"})
    assert bnd == [("chr17", 31902782, "chr17", 32067453, [("+", "-")], True)]      # deletion-type join
    mate = R.record_ends("chr17", 32067453, "]chr17:31902782]N", {"SVTYPE": "BND"})
    assert mate[0][4] == [("-", "+")]                    # the same join seen from the mate (sides at pos, mate)
    assert R.record_ends("chr1", 10, "<DEL>", {"SVTYPE": "DEL", "END": "500"})[0][4] == [("+", "-")]
    assert R.record_ends("chr1", 10, "<DUP>", {"SVTYPE": "DUP", "END": "500"})[0][4] == [("-", "+")]
    inv = R.record_ends("chr1", 10, "<INV>", {"SVTYPE": "INV", "END": "500", "STRAND": "+-"})
    assert inv[0][4] == [("+", "+"), ("-", "-")] and inv[0][5] is False              # both junctions, imprecise
    cute = R.record_ends("chr1", 10, "<INV>", {"SVTYPE": "INV", "END": "500", "STRAND": "++", "_caller": "cuteSV"})
    assert cute[0][4] == [("+", "+")] and cute[0][5] is True
    assert R.record_ends("chr1", 10, "<INS>", {"SVTYPE": "INS"}) == []


# ─── synthetic structures ───────────────────────────────────────────────────
def test_simple_deletion_called_del_and_inv_is_a_deletion():
    """BOC_4433-like: Sniffles2 DEL + SVIM BND pair (typed INV in the merged set), CN 1."""
    L = "chr17:100000-265000"
    rows = [jrow("chr17", 100000, "chr17", 265000, "DEL", "Sniffles2.DEL.1", L),
            jrow("chr17", 100000, "chr17", 265001, "INV", "svim.BND.1,svim.BND.2", L, "010")]
    recs = [rec("Sniffles2", "Sniffles2.DEL.1", "DEL", "chr17", 100000, "chr17", 265000, [("+", "-")], reads=15),
            rec("SVIM", "svim.BND.1", "BND", "chr17", 100000, "chr17", 265001, [("+", "-")], reads=14),
            rec("SVIM", "svim.BND.2", "BND", "chr17", 265001, "chr17", 100000, [("-", "+")], reads=14)]
    res = recon(L, rows, recs, depth("chr17", 100000, 265000, [(100000, 265000, 1.09)]))
    assert res["status"] == "resolved" and res["paths"] == {"der": ["p+", "q+"]}
    assert res["summary"].startswith("deletion of A")
    j = res["junctions"][0]
    assert j["orientations"] == ["+- deletion-type"] and j["callers"] == {"Sniffles2": 15, "SVIM": 14}
    note = " ".join(res["notes"])
    assert "label this junction differently" in note and "supports a deletion" in note


def test_del_vs_symbolic_inv_is_resolved_by_copy_number():
    """A DEL and a symbolic INV at the same breakends: a loss means deletion, normal CN inversion."""
    L = "chr5:1000000-1100000"
    rows = [jrow("chr5", 1000000, "chr5", 1100000, "DEL", "Sniffles2.DEL.9", L),
            jrow("chr5", 1000010, "chr5", 1100005, "INV", "svim.INV.9", L, "010")]
    recs = [rec("Sniffles2", "Sniffles2.DEL.9", "DEL", "chr5", 1000000, "chr5", 1100000, [("+", "-")]),
            rec("SVIM", "svim.INV.9", "INV", "chr5", 1000010, "chr5", 1100005, [("+", "+"), ("-", "-")], precise=False)]
    loss = recon(L, rows, recs, depth("chr5", 1000000, 1100000, [(1000000, 1100000, 1.0)]))
    assert loss["paths"]["der"] == ["p+", "q+"] and loss["status"] == "resolved"
    assert "supports a deletion" in " ".join(loss["notes"])
    kept = recon(L, rows, recs, depth("chr5", 1000000, 1100000, [(1000000, 1100000, 2.0)]))
    assert kept["paths"]["der"] == ["p+", "A-", "q+"]
    assert "as an inversion would" in " ".join(kept["notes"])


def test_del_inv_conflict_without_depth_is_ambiguous():
    L = "chr5:1000000-1000600"
    rows = [jrow("chr5", 1000000, "chr5", 1000600, "DEL", "Sniffles2.DEL.9", L),
            jrow("chr5", 1000000, "chr5", 1000600, "INV", "svim.INV.9", L, "010")]
    recs = [rec("Sniffles2", "Sniffles2.DEL.9", "DEL", "chr5", 1000000, "chr5", 1000600, [("+", "-")]),
            rec("SVIM", "svim.INV.9", "INV", "chr5", 1000000, "chr5", 1000600, [("+", "+"), ("-", "-")], precise=False)]
    res = recon(L, rows, recs, [])                                   # 600 bp: no depth estimate
    assert res["status"] == "ambiguous"
    assert any("read depth cannot decide" in r for r in res["reasons"])


def test_simple_inversion():
    L = "chr12:500000-560000"
    rows = [jrow("chr12", 500000, "chr12", 560000, "INV", "Sniffles2.INV.1", L)]
    recs = [rec("Sniffles2", "Sniffles2.INV.1", "INV", "chr12", 500000, "chr12", 560000, [("+", "+"), ("-", "-")],
                precise=False)]
    res = recon(L, rows, recs, depth("chr12", 500000, 560000, [(500000, 560000, 2.0)]))
    assert res["status"] == "resolved" and res["paths"]["der"] == ["p+", "A-", "q+"]


def test_del_inv_with_a_missing_junction_is_ambiguous_with_a_gap():
    """delINVdel with only one of its two inversion-type junctions called: open end, gap marker."""
    L = "chr17:1000-30000"
    rows = [jrow("chr17", 1000, "chr17", 25000, "INV", "Sniffles2.BND.1", L)]      # p+ | B-  only
    recs = [rec("Sniffles2", "Sniffles2.BND.1", "BND", "chr17", 1000, "chr17", 25000, [("+", "+")])]
    rows.append(jrow("chr17", 5000, "chr17", 30000, "DEL", "Sniffles2.DEL.2", L))  # spans the second junction area
    recs.append(rec("Sniffles2", "Sniffles2.DEL.2", "DEL", "chr17", 5000, "chr17", 30000, [("+", "-")]))
    res = recon(L, rows, recs, depth("chr17", 1000, 30000, [(1000, 5000, 1.0), (5000, 25000, 2.0), (25000, 30000, 1.0)]))
    assert res["status"] == "ambiguous"
    assert any("could not be completed" in r or "not on the drawn derivative" in r for r in res["reasons"])


def test_translocation_is_drawn_as_two_derivatives():
    L = "chr2:199854292-199856296"
    rows = [jrow("chr2", 199854293, "chr7", 115869905, "BND", "Sniffles2.BND.1", L),
            jrow("chr2", 199856293, "chr7", 115868955, "BND", "Sniffles2.BND.2", L)]
    recs = [rec("Sniffles2", "Sniffles2.BND.1", "BND", "chr2", 199854293, "chr7", 115869905, [("+", "-")]),
            rec("Sniffles2", "Sniffles2.BND.2", "BND", "chr2", 199856293, "chr7", 115868955, [("-", "+")])]
    res = recon(L, rows, recs, [])
    assert res["inter"] and res["status"] == "resolved"
    assert res["paths"] == {"der(2)": ["p_chr2+", "q_chr7+"], "der(7)": ["p_chr7+", "q_chr2+"]}
    assert "joins between chromosomes chr2/chr7" in res["summary"]


def test_no_junction_gives_no_structure():
    res = recon("chr1:0-100000", [], [], [])
    assert res["status"] == "ambiguous" and not res["paths"]


# ─── truth-derived inputs (positive-control design) ─────────────────────────
def truth_inputs(ev):
    """Junctions (one precise BND record per novel join), locus depth and the locus id of a
    planted event, derived from its truth path only."""
    segs, paths = ev["segments"], ev["paths"]
    visits = Counter(V.parse(e)[0] for p in paths.values() for e in p)
    rows, recs, k = [], [], 0
    lo = {}
    for n, (c, a, b) in segs.items():
        if not V.is_flank(n):
            lo.setdefault(c, [a, b])
            lo[c] = [min(lo[c][0], a), max(lo[c][1], b)]
    locus = "truth:" + ev["id"]
    for p in paths.values():
        for x, y in zip(p, p[1:]):
            (a, sa), (b, sb) = V.parse(x), V.parse(y)
            A, B = segs[a], segs[b]
            if (sa == sb == "+" and A[0] == B[0] and A[2] == B[1]) or (sa == sb == "-" and A[0] == B[0] and A[1] == B[2]):
                continue
            e1 = (A[0], A[2], "+") if sa == "+" else (A[0], A[1], "-")
            e2 = (B[0], B[1], "-") if sb == "+" else (B[0], B[2], "+")
            if (R._natural(e2[0]), e2[1]) < (R._natural(e1[0]), e1[1]):
                e1, e2 = e2, e1
            o = (e1[2], e2[2])
            t = "BND" if e1[0] != e2[0] else {("+", "-"): "DEL", ("-", "+"): "DUP"}.get(o, "INV")
            k += 1
            rid = f"truth.BND.{k}"
            rows.append(jrow(e1[0], e1[1], e2[0], e2[1], t, rid, locus))
            recs.append(rec("Sniffles2", rid, "BND", e1[0], e1[1], e2[0], e2[1], [o]))
    drows = []
    for c, (a, b) in lo.items():
        drows += [(c, a - 21000, a - 1000, "flank_left", BASE), (c, b + 1000, b + 21000, "flank_right", BASE)]
    for n, (c, a, b) in segs.items():
        if not V.is_flank(n):
            drows.append((c, a + 200, b - 200, "segment", BASE * (visits.get(n, 0) + 1) / 2))
    return locus, rows, recs, drows


EVENTS = [e for e in yaml.safe_load((REPO / "tests" / "positive_control" / "events.yaml").read_text())["events"]
          if e.get("paths")]


@pytest.mark.parametrize("ev", EVENTS, ids=[e["id"] for e in EVENTS])
def test_truth_derived_reconstruction_reproduces_the_planted_layout(ev):
    locus, rows, recs, drows = truth_inputs(ev)
    res = R.reconstruct({"Locus_ID": locus}, rows, recs, {}, drows, [], NoGenes(), PARAMS, "female")
    cmp = V.compare(res, ev)
    assert cmp["outcome"] in ("exact", "in_alternatives"), (cmp["truth"], cmp["drawn"], cmp["alternatives"])
    if ev["id"] in ("E01", "E02", "E04"):
        assert cmp["outcome"] == "exact", (cmp["truth"], cmp["drawn"])


def test_merge_unbroken_truth_segments():
    E04 = next(e for e in EVENTS if e["id"] == "E04")
    segs, paths, _ = V.merge_unbroken(E04["segments"], E04["paths"])
    assert "A" not in segs and segs["p"][2] == E04["segments"]["A"][2]     # A never leaves the left flank
    assert [f"{n}{s}" for n, s in paths["der"]] == ["p+", "C-", "D+", "E+", "E+", "q+"]


# ─── figure and report ──────────────────────────────────────────────────────
def test_recon_figure_and_layout(tmp_path):
    pytest.importorskip("matplotlib")
    import cxsv_recon_figure as RF
    E01 = next(e for e in EVENTS if e["id"] == "E01")
    locus, rows, recs, drows = truth_inputs(E01)
    res = R.reconstruct({"Locus_ID": locus}, rows, recs, {}, drows, [], NoGenes(), PARAMS, "female")
    ev = R.layout(res, title="Tier 1 #1 test", subtitle="MECP2")
    assert set(ev["paths"]) == {"der"} and ev["junction_labels"]["der"].count(None) < len(ev["paths"]["der"])
    f = RF.recon_figure(ev, str(tmp_path / "r.png"))
    assert (tmp_path / "r.png").stat().st_size > 5000 and f["width_px"] > f["height_px"]
    tr = R.reconstruct({"Locus_ID": "x"}, [jrow("chr2", 1000, "chr7", 5000, "BND", "b1", "x")],
                       [rec("Sniffles2", "b1", "BND", "chr2", 1000, "chr7", 5000, [("+", "-")])], {}, [], [],
                       NoGenes(), PARAMS)
    RF.recon_figure(R.layout(tr, title="t"), str(tmp_path / "t.png"))      # inter-chromosomal, not to scale
    assert (tmp_path / "t.png").exists()


def test_report_section_and_page_rule():
    import clinician_report as cr
    assert cr.page_limit(0) == 2 and cr.page_limit(1) == 3 and cr.page_limit(2) == 3 and cr.page_limit(3) == 4
    assert cr.page_limit(3, per_page=3, base_pages=7) == 8
    L = "chr17:100000-265000"
    rows = [jrow("chr17", 100000, "chr17", 265000, "DEL", "Sniffles2.DEL.1", L)]
    recs = [rec("Sniffles2", "Sniffles2.DEL.1", "DEL", "chr17", 100000, "chr17", 265000, [("+", "-")], reads=15)]
    res = recon(L, rows, recs, depth("chr17", 100000, 265000, [(100000, 265000, 1.0)]))
    row = {"Locus_ID": L, "Chr": "chr17", "Start": "100000", "End": "265000"}
    B = cr.recon_blocks([{"rank": 1, "row": row, "res": res, "figure": {"path": "x.png", "width_px": 2, "height_px": 1}}], 1)
    assert B[0] == ("h2", "Structural reconstruction of the Tier 1 loci")
    fig = next(b for b in B if b[0] == "recon_figure")
    assert fig[2].startswith(cr.RECON_CAPTION)
    assert cr.RECON_CAPTION == "Inferred from split-read junctions and read depth in this sample; not orthogonally confirmed."
    kv = dict((r[0], r[1]) for r in next(b for b in B if b[0] == "kv2")[1])
    assert "Sniffles2 15 reads" in kv["Junctions (support)"] and kv["Inferred structure"].startswith("deletion of A")


def test_zero_tier1_has_no_reconstruction():
    import clinician_report as cr

    class A:
        junctions = "junctions.tsv"
    assert cr.make_reconstructions(A(), {"T": {}, "P": {}}, []) == []


def test_reconstruction_code_is_not_fingerprinted():
    import cxsv_fingerprint as fpm
    assert not any("recon" in f or "report" in f for f in fpm.CALLING_CODE)
