"""CxSV locus structure and classes on distinct junctions (workflow/lib/cxsv_classes.py).

Coordinates of several cases are those of the positive-control calls
(tests/lab/pos_control_regression.sh, item-3 run): E12 (simple inversion,
Sniffles2 reciprocal records), E02 (delINVdel with an SVIM composite), E07
(reciprocal translocation with junction deletions)."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "workflow", "lib"))
import cxsv_classes as cc  # noqa: E402
import cxsv_junctions as cj  # noqa: E402

ORDER = {f"chr{i}": i for i in range(1, 23)} | {"chrX": 23}


def J(c1, p1, c2, p2, t, supp="001", rid=None):
    return cj.Junction(c1, p1, c2, p2, t, abs(p2 - p1) if t != "BND" else 0, supp, 1,
                       [rid or f"{t}{p1}_{p2}"])


def analyse(js, chrom="chr1", **kw):
    return cc.analyse_locus(js, chrom, ORDER, kw or None)


def test_simple_inversion_reciprocal_records_is_one_junction_not_cxsv():
    # E12: Sniffles2.BND.81SC / 83SC, identical after SURVIVOR
    res = analyse([J("chr12", 74754175, "chr12", 74814325, "INV", "111", "a"),
                   J("chr12", 74754175, "chr12", 74814325, "INV", "001", "b")], "chr12")
    assert res["n_records"] == 2 and res["n_junctions"] == 1 and res["n_breakends"] == 2
    assert not res["is_cxsv"] and res["class"] == "NA" and res["tier"] == "simple"


def test_reciprocal_records_within_collapse_distance():
    res = analyse([J("chr1", 1000, "chr1", 9000, "INV", rid="a"), J("chr1", 1002, "chr1", 9001, "INV", rid="b")])
    assert res["n_junctions"] == 1 and not res["is_cxsv"]


def test_del_inv_del_with_svim_composite():
    # E02: tt 32243849-32265774, hh 32247849-32271849, SVIM composite 32246057-32269126
    js = [J("chr17", 32243849, "chr17", 32265774, "INV", "101", "s140"),
          J("chr17", 32265774, "chr17", 32243849, "INV", "001", "s142"),
          J("chr17", 32246057, "chr17", 32269126, "INV", "010", "svim19"),
          J("chr17", 32247849, "chr17", 32271849, "INV", "101", "s141"),
          J("chr17", 32247851, "chr17", 32271850, "INV", "001", "s143")]
    js = [cj.Junction(*cj.ordered_ends(j.chrom1, j.pos1, j.chrom2, j.pos2, ORDER), j.svtype, j.svlen,
                      j.supp_vec, 1, j.record_ids) for j in js]
    res = analyse(js, "chr17")
    assert res["n_junctions"] == 2 and res["n_composite"] == 1
    assert [x.record_ids for x in res["composite"]] == [["svim19"]]
    assert res["n_breakends"] == 4 and res["max_interleaved"] == 2
    assert res["criteria"] == ["C1_min_breakpoints"]
    assert res["class"] == "Complex_inversion" and res["tier"] == "intermediate"


def test_composite_needs_single_caller_and_midpoint():
    a, b = J("chr1", 1000, "chr1", 20000, "INV", "101"), J("chr1", 5000, "chr1", 25000, "INV", "101")
    mid1 = J("chr1", 3000, "chr1", 22500, "INV", "010")
    assert cc.composite_of(mid1, [a, b, mid1], 300, 0.25) is not None
    two_callers = J("chr1", 3000, "chr1", 22500, "INV", "110")
    assert cc.composite_of(two_callers, [a, b, two_callers], 300, 0.25) is None
    off = J("chr1", 1500, "chr1", 22500, "INV", "010")      # pos1 far from the midpoint
    assert cc.composite_of(off, [a, b, off], 300, 0.25) is None
    apart = [J("chr1", 1000, "chr1", 3000, "INV", "101"), J("chr1", 20000, "chr1", 26000, "INV", "101")]
    between = J("chr1", 10000, "chr1", 14000, "INV", "010")   # a and b do not overlap: not a composite
    assert cc.composite_of(between, apart + [between], 300, 0.25) is None
    other_type = J("chr1", 3000, "chr1", 22500, "DEL", "010")
    assert cc.composite_of(other_type, [a, b, other_type], 300, 0.25) is None


def test_two_nearby_deletions_are_not_a_cxsv():
    res = analyse([J("chr1", 1000, "chr1", 3000, "DEL", "111"), J("chr1", 20000, "chr1", 26000, "DEL", "111")])
    assert res["n_breakends"] == 4 and not res["is_cxsv"]


def test_three_deletions_are_c1():
    res = analyse([J("chr1", 1000, "chr1", 3000, "DEL"), J("chr1", 9000, "chr1", 9800, "DEL"),
                   J("chr1", 20000, "chr1", 26000, "DEL")])
    assert res["criteria"][0] == "C1_min_breakpoints" and res["class"] == "Multi_breakpoint_cluster"


def test_translocation_balanced_vs_unbalanced():
    # balanced: der(2) and der(7) junctions share their breakpoints -> one junction
    bal = analyse([J("chr2", 199854293, "chr7", 115869903, "BND", "110", "t1"),
                   J("chr2", 199854294, "chr7", 115869904, "BND", "010", "t2")], "chr2")
    assert bal["n_junctions"] == 1 and not bal["is_cxsv"]
    # E07: junction deletions of 2.0 kb (chr2) and 0.95 kb (chr7)
    unb = analyse([J("chr2", 199854293, "chr7", 115869903, "BND", "110", "t1"),
                   J("chr2", 199856295, "chr7", 115868955, "BND", "010", "t2")], "chr2")
    assert unb["n_junctions"] == 2 and unb["n_breakends"] == 4 and unb["n_inter"] == 2
    assert unb["criteria"] == ["C1_min_breakpoints"] and unb["class"] == "BND_cluster"


def test_interleaving_relation():
    a, b = J("chr1", 100, "chr1", 1000, "INV"), J("chr1", 500, "chr1", 1500, "INV")
    assert cc.interleaved(a, b, 50) and cc.interleaved(b, a, 50)
    nested = J("chr1", 200, "chr1", 800, "INV")
    assert not cc.interleaved(a, nested, 50)
    shared_end = J("chr1", 500, "chr1", 1020, "INV")           # ends within tolerance: not interleaved
    assert not cc.interleaved(a, shared_end, 50)
    disjoint = J("chr1", 2000, "chr1", 3000, "INV")
    assert not cc.interleaved(a, disjoint, 50)


def test_chromothripsis_needs_extreme_tier_and_interleaving():
    # 6 interleaved inversions (12 breakends): Chromothripsis_like
    js = [J("chr1", 10000 * i + 1000, "chr1", 10000 * i + 25000, "INV", "101") for i in range(6)]
    res = analyse(js)
    assert res["tier"] == "extreme" and res["max_interleaved"] == 6
    assert res["class"] == "Chromothripsis_like"
    # the same breakend count without interleaving (nested): Multi_breakpoint_cluster
    nested = [J("chr1", 1000 + 3000 * i, "chr1", 60000 - 3000 * i, "INV", "101") for i in range(3)]
    res = analyse(nested)
    assert res["tier"] == "extreme" and res["max_interleaved"] == 1
    assert res["class"] == "Multi_breakpoint_cluster"


def test_compact_c5_event_is_not_chromothripsis_like():
    # 3 junctions in 20 kb fire C5, but 4-5 breakends are the intermediate tier
    js = [J("chr1", 1000, "chr1", 9000, "INV"), J("chr1", 1000, "chr1", 6000, "INV"),
          J("chr1", 15000, "chr1", 21000, "DUP")]
    res = analyse(js)
    assert "C5_chromoanagenesis_cluster" in res["criteria"]
    assert res["tier"] == "intermediate" and res["class"] == "DUP+INV"


def test_class_precedence_type_labels():
    f = cc.complexity_class
    assert f(["DEL", "DUP", "INV"], 5, ["C1_min_breakpoints"], 1, 0, 0, 1) == "DEL+DUP+INV"
    assert f(["DEL", "INV"], 4, ["C4_CN_plus_orient"], 1, 0, 0, 1) == "DEL+INV"
    assert f(["INS", "DEL"], 3, ["C3_overlapping_types"], 0, 0, 0, 0) == "INS+DEL"
    assert f(["BND"], 4, ["C1_min_breakpoints"], 0, 2, 2, 0) == "BND_cluster"
    assert f(["INV"], 5, ["C1_min_breakpoints", "C2_nested_SV"], 1, 0, 0, 3) == "Nested_SV"
    assert f(["DEL", "BND"], 3, ["C4_CN_plus_orient"], 0, 1, 1, 0) == "CN_plus_orientation_change"
    assert f(["INV"], 4, ["C1_min_breakpoints"], 2, 0, 0, 2) == "Complex_inversion"
    assert f(["DEL"], 4, ["C1_min_breakpoints"], 0, 0, 0, 0) == "Complex_other"
    # the alternative interleaving path: >= 3 interleaved + >= 4 interchromosomal
    p = {"min_interleaved": 6}
    assert f(["INV", "BND"], 9, ["C1_min_breakpoints"], 3, 4, 4, 3, p) == "Chromothripsis_like"
    assert f(["INV", "BND"], 9, ["C1_min_breakpoints"], 3, 3, 3, 3, p) == "BND_cluster"
    for c in ("Overlapping_types",):
        assert c in cc.CLASSES


def test_breakends_merge_adjacent_segment_ends():
    # E09-like: the end of one segment and the start of the next (1 bp apart) are one break
    js = [J("chr1", 1000, "chr1", 50000, "INV"), J("chr1", 1001, "chr1", 30000, "INV")]
    assert len(cc.breakends(js, 300)) == 3
