"""Panel junction-site clustering and breakpoint-level matching (workflow/lib/cxsv_matching.py)."""
import cxsv_matching as cm

ORDER = {"chr1": 0, "chr2": 1, "chrX": 2}


def obs(ind, pos1, pos2, svtype="DEL", chrom="chr1", svlen=None, dosage=1):
    return {"svtype": svtype, "chrom1": chrom, "pos1": pos1, "chrom2": chrom, "pos2": pos2,
            "svlen": svlen if svlen is not None else abs(pos2 - pos1), "dosage": dosage,
            "individual": ind, "sample": ind}


def index_of(sites):
    for s in sites:
        s["carrier_set"] = {m["individual"] for m in s["members"]}
    return cm.SiteIndex(sites)


def test_clustering_is_type_aware_leader_based_and_counts_individuals():
    o = [obs("A", 1000, 5000), obs("B", 1040, 5010), obs("A", 1010, 4990),      # A twice: one carrier
         obs("C", 1450, 5400),                                                  # within 500 of A's leader
         obs("D", 1600, 5600),                                                  # 600 from the leader: new site
         obs("E", 1000, 5000, svtype="DUP")]                                    # same place, other type
    sites = cm.cluster_sites(o, 500, 0.5, ORDER)
    dels = [s for s in sites if s["svtype"] == "DEL"]
    assert [len({m["individual"] for m in s["members"]}) for s in dels] == [3, 1]
    assert (dels[0]["min_pos1"], dels[0]["max_pos1"]) == (1000, 1450)
    assert dels[0]["pos1"] == 1010                    # lower median of 1000, 1010, 1040, 1450
    assert [s["svtype"] for s in sites].count("DUP") == 1


def test_allele_counts_use_dosage_and_sex_ploidy():
    inds = {"M": {"sex": "male", "superpopulation": "EUR", "platforms": {"ONT"}},
            "F": {"sex": "female", "superpopulation": "AFR", "platforms": {"ONT"}},
            "G": {"sex": "female", "superpopulation": "AFR", "platforms": {"PACBIO"}}}
    groups = ["platform:ONT", "platform:PACBIO", "superpopulation:AFR", "superpopulation:EUR"]
    auto = cm.cluster_sites([obs("M", 100, 900, dosage=2), obs("F", 100, 900)], 500, 0.5, ORDER)[0]
    k = cm.site_counts(auto, inds, groups)
    assert (k["ac"], k["an"], k["carriers"]) == (3, 6, ["F", "M"])
    assert (k["groups"]["platform:PACBIO"]["ac"], k["groups"]["platform:PACBIO"]["an"]) == (0, 2)
    x = cm.cluster_sites([obs("M", 100, 900, chrom="chrX", dosage=2)], 500, 0.5, ORDER)[0]
    k = cm.site_counts(x, inds, groups)
    assert (k["ac"], k["an"]) == (1, 5)                # hemizygous male: dosage capped at 1; an 1+2+2


def test_matching_respects_type_tolerance_and_ins_length():
    sites = cm.cluster_sites([obs("A", 1000, 5000), obs("B", 3000, 3000, svtype="INS", svlen=300)], 500, 0.5, ORDER)
    idx = index_of(sites)
    q = lambda **kw: {"svtype": "DEL", "chrom1": "chr1", "pos1": 1000, "chrom2": "chr1", "pos2": 5000,
                      "svlen": 4000, **kw}
    assert cm.junction_carriers(q(pos1=1400), idx, 500, 0.5)[0] == {"A"}
    assert cm.junction_carriers(q(pos2=5600), idx, 500, 0.5)[0] == set()
    assert cm.junction_carriers(q(svtype="DUP"), idx, 500, 0.5)[0] == set()
    ins = q(svtype="INS", pos1=3100, pos2=3100)
    assert cm.junction_carriers(dict(ins, svlen=280), idx, 500, 0.5)[0] == {"B"}
    assert cm.junction_carriers(dict(ins, svlen=1000), idx, 500, 0.5)[0] == set()


def locus(*ends, svtype="INV"):
    return [{"svtype": svtype, "chrom1": "chr1", "pos1": a, "chrom2": "chr1", "pos2": b, "svlen": b - a}
            for a, b in ends]


def test_hotspot_decoy_stays_private_and_common_event_is_public():
    # panel: A and B carry a 3-junction rearrangement (INV) at 100-130 kb; C carries a
    # DIFFERENT rearrangement at the same place (decoy: other breakpoints)
    common = [(100_000, 120_000), (101_000, 125_000), (110_000, 130_000)]
    decoy = [(103_000, 118_000), (112_000, 127_000)]
    o = [obs(i, a, b, svtype="INV") for i in ("A", "B") for a, b in common]
    o += [obs("C", a, b, svtype="INV") for a, b in decoy]
    idx = index_of(cm.cluster_sites(o, 500, 0.5, ORDER))
    pub = cm.classify_locus(locus(*common), idx, 500, 0.5, 0.8, 1)
    assert (pub["status"], pub["carriers"], pub["n_seen"]) == ("public", ["A", "B"], 3)
    priv = cm.classify_locus(locus((100_200, 119_000), (106_000, 128_500)), idx, 500, 0.5, 0.8, 1)
    assert (priv["status"], priv["carriers"], priv["n_seen"]) == ("private", [], 0)
    # patient shares 1 of 3 junctions with A/B: seen, but A/B do not carry the rearrangement
    part = cm.classify_locus(locus(common[0], (140_000, 150_000), (145_000, 160_000)), idx, 500, 0.5, 0.8, 1)
    assert (part["status"], part["carriers"], round(part["best_fraction"], 3)) == ("shared_junctions", [], 0.333)


def test_fraction_threshold_and_min_carriers():
    js = [(10_000 * k, 10_000 * k + 5000) for k in range(1, 6)]           # 5 junctions
    o = [obs("A", a, b) for a, b in js[:4]] + [obs("B", a, b) for a, b in js]
    idx = index_of(cm.cluster_sites(o, 500, 0.5, ORDER))
    L = locus(*js, svtype="DEL")
    assert cm.classify_locus(L, idx, 500, 0.5, 0.8, 1)["carriers"] == ["A", "B"]    # 4/5 = 0.8 counts
    assert cm.classify_locus(L, idx, 500, 0.5, 1.0, 1)["carriers"] == ["B"]
    assert cm.classify_locus(L, idx, 500, 0.5, 1.0, 2)["status"] == "shared_junctions"
    assert cm.classify_locus([], idx, 500, 0.5, 0.8, 1)["status"] == "private"
