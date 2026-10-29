"""Junction model (workflow/lib/cxsv_junctions.py): extraction, ordering, collapsing, loci."""
import cxsv_junctions as cj

ORDER = {"chr1": 0, "chr2": 1, "chr7": 2, "chr22": 3, "chrX": 4}
FMT = "GT:DR"


def rec(chrom, pos, vid, alt, info, samples=("./.:0", "./.:0", "0/1:5")):
    return [chrom, str(pos), vid, "N", alt, "60", "PASS", info, FMT, *samples]


def test_record_types_and_schema_conventions():
    d = cj.record_to_junction(rec("chr1", 1000, "d", "<DEL>", "SVTYPE=DEL;END=3000;SVLEN=-2000;SUPP_VEC=111"), ORDER)
    assert (d.chrom1, d.pos1, d.chrom2, d.pos2, d.svtype, d.svlen) == ("chr1", 1000, "chr1", 3000, "DEL", 2000)
    i = cj.record_to_junction(rec("chr1", 5000, "i", "<INS>", "SVTYPE=INS;END=5000;SVLEN=320"), ORDER)
    assert (i.pos1, i.pos2, i.svlen) == (5000, 5000, 320)
    # SURVIVOR TRA (CHR2/END + bracket ALT) becomes BND; ends ordered by .fai contig order
    t = cj.record_to_junction(rec("chr7", 115868955, "t", "N[chr2:199856295[",
                                  "SVTYPE=TRA;CHR2=chr2;END=199856295;SVLEN=0"), ORDER)
    assert (t.chrom1, t.pos1, t.chrom2, t.pos2, t.svtype, t.svlen) == ("chr2", 199856295, "chr7", 115868955, "BND", 0)
    # an INV written as a reciprocal breakend (END < POS) is re-ordered
    v = cj.record_to_junction(rec("chr22", 51223041, "v", "[chr22:51213667[T",
                                  "SVTYPE=INV;CHR2=chr22;END=51213667;SVLEN=9374"), ORDER)
    assert (v.pos1, v.pos2, v.svtype) == (51213667, 51223041, "INV")
    assert cj.record_to_junction(rec("chr1", 1, "x", "<CNV>", "SVTYPE=FOO"), ORDER) is None


def test_dosage_prefers_sniffles_then_cutesv():
    assert cj.dosage_from_samples(FMT, ["0/1:3", "1/1:4", "1/1:9"]) == 2           # Sniffles2 column
    assert cj.dosage_from_samples(FMT, ["0/1:3", "1/1:4", "./.:0"]) == 1           # falls back to cuteSV
    assert cj.dosage_from_samples(FMT, ["./.:0", "./.:0", "./.:0"]) == 1           # merged call, no GT


def test_reciprocal_inversion_breakends_collapse_but_distinct_junctions_do_not():
    recs = [
        rec("chr22", 51213665, "Sniffles2.BND.9BS16", "[chr22:51223042[T",
            "SVTYPE=INV;CHR2=chr22;END=51223040;SVLEN=9376;SUPP_VEC=101"),
        rec("chr22", 51223041, "Sniffles2.BND.9ES16", "[chr22:51213667[T",
            "SVTYPE=INV;CHR2=chr22;END=51213667;SVLEN=9374;SUPP_VEC=001"),
        rec("chr22", 51223040, "Sniffles2.BND.9DS16", "G]chr22:51209790]",
            "SVTYPE=INV;CHR2=chr22;END=51209790;SVLEN=13250;SUPP_VEC=001"),
        rec("chr22", 51225790, "dup", "<DUP>", "SVTYPE=DUP;END=51230915;SVLEN=5108;SUPP_VEC=011"),
    ]
    js = cj.collapse([cj.record_to_junction(r, ORDER) for r in recs], 300, ORDER)
    got = [(j.svtype, j.pos1, j.pos2, j.supp_vec, len(j.record_ids)) for j in js]
    assert got == [("INV", 51209790, 51223040, "001", 1),
                   ("INV", 51213665, 51223040, "101", 2),       # leader = best-supported record
                   ("DUP", 51225790, 51230915, "011", 1)]


def test_collapse_is_leader_based_and_type_aware():
    def j(pos, svtype="DEL", supp="001"):
        return cj.Junction("chr1", pos, "chr1", pos + 5000, svtype, 5000, supp, 1, [f"{svtype}{pos}"])
    # 0 claims 250; 500 is > 300 from the leader, so it is not chained in
    js = cj.collapse([j(1000), j(1250), j(1500)], 300, ORDER)
    assert [x.pos1 for x in js] == [1000, 1500]
    # same coordinates, different types: never one junction
    js = cj.collapse([j(1000, "DEL"), j(1000, "DUP")], 300, ORDER)
    assert sorted(x.svtype for x in js) == ["DEL", "DUP"]
    # INS of very different lengths at one position stay separate
    a = cj.Junction("chr1", 100, "chr1", 100, "INS", 300, "001", 1, ["a"])
    b = cj.Junction("chr1", 110, "chr1", 110, "INS", 2000, "001", 1, ["b"])
    assert len(cj.collapse([a, b], 300, ORDER)) == 2


def test_cxsv_locus_membership_from_clusters(tmp_path):
    summary = tmp_path / "s.tsv"
    summary.write_text("ClusterID\tChr\tStart\tEnd\tSize\tN_Breakpoints\tN_SVtypes\tSV_Types\tHas_Inter\tIsCxSV\t"
                       "CxSV_Criteria\n"
                       "1\tchr2\t199854292\t199856296\t2004\t3\t1\tTRA\tYES\tYES\tC1\n"
                       "2\tchr7\t115868954\t115869904\t950\t3\t1\tTRA\tYES\tYES\tC1\n"
                       "3\tchr1\t100\t200\t100\t1\t1\tDEL\tNO\tNO\tNONE\n")
    clusters = tmp_path / "c.bed"
    clusters.write_text("chr2\t199854292\t199854294\tTRA\tt1\tinter\t1\n"
                        "chr7\t115869902\t115869904\tTRA\tt1_dest\tinter\t2\n"
                        "chr1\t100\t200\tDEL\td1\tintra\t3\n")
    m = cj.cxsv_loci_of_records(clusters, summary)
    assert m == {"t1": {"chr2:199854292-199856296", "chr7:115868954-115869904"}}


def test_sample_junctions_end_to_end(tmp_path):
    vcf = tmp_path / "m.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\tc\n"
                   + "\t".join(rec("chr1", 1000, "d1", "<DEL>", "SVTYPE=DEL;END=3000;SUPP_VEC=111")) + "\n"
                   + "\t".join(rec("chr1", 1001, "d2", "<DEL>", "SVTYPE=DEL;END=3001;SUPP_VEC=001")) + "\n")
    summ = tmp_path / "s.tsv"
    summ.write_text("ClusterID\tChr\tStart\tEnd\tIsCxSV\n7\tchr1\t1000\t3001\tYES\n")
    clu = tmp_path / "c.bed"
    clu.write_text("chr1\t1000\t3000\tDEL\td1\tintra\t7\nchr1\t1001\t3001\tDEL\td2\tintra\t7\n")
    js = cj.sample_junctions(vcf, dict(ORDER), 300, clu, summ)
    assert len(js) == 1 and js[0].record_ids == ["d1", "d2"] and js[0].locus_ids == {"chr1:1000-3001"}
    out = tmp_path / "j.tsv"
    cj.write_table(out, js)
    back = cj.read_table(out)
    assert back[0]["pos1"] == 1000 and back[0]["loci"] == {"chr1:1000-3001"} and back[0]["supp_vec"] == "111"


def test_orientation_follows_fai_order_not_numeric(tmp_path):
    fai = tmp_path / "hs1.fa.fai"           # hs1 .fai is ordered by length: ... chr7, chrX, chr9 ...
    fai.write_text("chr7\t160567428\t1\t60\t61\nchrX\t154259566\t2\t60\t61\nchr9\t150617247\t3\t60\t61\n")
    order = cj.contig_order(fai)
    b = cj.record_to_junction(rec("chr9", 1000, "b", "N[chrX:5000[", "SVTYPE=BND"), order)
    assert (b.chrom1, b.pos1, b.chrom2, b.pos2) == ("chrX", 5000, "chr9", 1000)
