"""Paired-breakend linkage of long records (cxsv_params.long_span_rule paired;
roadmap item 13 follow-up, docs/PAIRED_LINKAGE.md).

A DEL/DUP/INV longer than span_link_max takes part in discovery only when both of
its ends lie within cluster_distance of another record (not a duplicate of the same
junction, and a long partner must itself be paired); it then joins the two end
clusters into one cluster sized by its parts. Otherwise it is left out, as under
"exclude". bedtools is not available in the unit job: cluster() re-implements
`bedtools cluster -d` as in test_discovery_linkage.py.
"""
import csv
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(REPO / "workflow" / "lib"))
import link_paired_clusters as lpc  # noqa: E402
import vcf_to_bed as v2b            # noqa: E402

LINK = 500_000
HEADER = ("##fileformat=VCFv4.1\n##contig=<ID=chr17,length=84276897>\n"
          "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNULL\tSample\tS\n")
GT = "GT:PSV:LN:DR:ST:QV:TY:ID:RAL:AAL:CO"


def rec(chrom, pos, vid, svtype, end, supp="011", strands="+-"):
    svlen = end - pos if svtype != "DEL" else -(end - pos)
    col = f"0/1:NA:{abs(svlen)}:0,0:{strands}:20:{svtype}:{vid}:N:N:{chrom}_{pos}-{chrom}_{end}"
    empty = "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN"
    cols = [col if b == "1" else empty for b in supp]
    info = (f"SUPP={supp.count('1')};SUPP_VEC={supp};SVLEN={svlen};SVTYPE={svtype};"
            f"SVMETHOD=SURVIVOR1.0.7;CHR2={chrom};END={end};CIPOS=0,0;CIEND=0,0;STRANDS={strands}")
    return "\t".join([chrom, str(pos), vid, "N", f"<{svtype}>", "20", "PASS", info, GT] + cols)


def write(tmp, name, lines):
    p = tmp / name
    p.write_text(HEADER + "".join(l + "\n" for l in lines))
    return p


def cluster(bed_lines, d=50_000):
    rows = sorted((l.split("\t") for l in bed_lines), key=lambda f: (f[0], int(f[1]), int(f[2])))
    out, cid, chrom, max_end = [], 0, None, None
    for f in rows:
        s, e = int(f[1]), int(f[2])
        if f[0] != chrom or s - max_end > d:
            cid, chrom, max_end = cid + 1, f[0], e
        else:
            max_end = max(max_end, e)
        out.append("\t".join(f + [str(cid)]))
    return out


def run(tmp, vcf, rule):
    """vcf_to_bed (rule) -> cluster -> link_paired_clusters -> classify_cxsv; returns
    (long rows, summary rows)."""
    long_rows = []
    bed = list(v2b.process_vcf(str(vcf), 1, LINK, {}, rule, long_rows, 50_000, 300))
    longf = tmp / f"long.{rule}.tsv"
    longf.write_text("ID\tCHROM\tPOS\tEND\tSVTYPE\tLENGTH\tSUPP_VEC\tRULE\n"
                     + "".join("\t".join(map(str, r)) + "\n" for r in long_rows))
    linked, _, _ = lpc.link(cluster(bed), lpc.read_paired(longf))
    cl = tmp / f"clusters.{rule}.bed"
    cl.write_text("".join(l + "\n" for l in linked))
    fai = tmp / "ref.fai"
    fai.write_text("chr17\t84276897\t7\t60\t61\n")
    out = tmp / f"summary.{rule}.tsv"
    subprocess.run([sys.executable, str(SCRIPTS / "classify_cxsv.py"), "--bed", str(cl), "--vcf", str(vcf),
                    "--fai", str(fai), "--output", str(out), "--max-cluster-size", "500000"], check=True,
                   capture_output=True)
    with open(out, newline="") as fh:
        return long_rows, list(csv.DictReader(fh, delimiter="\t"))


# DUP-NML-DUP / inserted-duplication pattern: two long junctions whose ends pair up
# 30 kb apart at both sides (literature events NP18_P06, BAH19_BAB8923 have this shape)
PAIR = [rec("chr17", 10_000_000, "svim.DUP.a", "DUP", 14_500_000, "011", "-+"),
        rec("chr17", 10_030_000, "svim.DUP.b", "DUP", 14_470_000, "011", "-+")]
GIANT = rec("chr17", 16_716_382, "Sniffles2.DEL.giant", "DEL", 42_430_144, "001")
E02 = [rec("chr17", 32243849, "Sniffles2.DEL.e1", "DEL", 32247849, "111"),
       rec("chr17", 32247849, "Sniffles2.INV.e2", "INV", 32267849, "011", "++"),
       rec("chr17", 32267849, "Sniffles2.DEL.e3", "DEL", 32271849, "111")]
FILLER = [rec("chr17", p, f"cuteSV.DEL.f{p}", "DEL", p + 300, "100") for p in range(17_000_000, 42_000_000, 400_000)
          if abs(p - 32_250_000) > 150_000]


def test_paired_long_records_join_their_end_clusters(tmp_path):
    vcf = write(tmp_path, "m.vcf", PAIR)
    long_rows, summ = run(tmp_path, vcf, "paired")
    assert [(r[0], r[7]) for r in long_rows] == [("svim.DUP.a", "paired"), ("svim.DUP.b", "paired")]
    assert len(summ) == 1                                         # one cluster, not one per end
    loc = summ[0]
    assert loc["IsCxSV"] == "YES" and loc["N_Records"] == "2"
    assert (int(loc["Start"]), int(loc["End"])) == (9_999_999, 14_500_001)
    assert int(loc["Size"]) == (10_030_001 - 9_999_999) + (14_500_001 - 14_469_999)   # parts, not span
    # the same records under exclude take no part; under breakends each end is its own locus
    assert run(tmp_path, vcf, "exclude")[1] == []
    assert len(run(tmp_path, vcf, "breakends")[1]) == 2


def test_long_record_without_partners_is_handled_as_exclude(tmp_path):
    vcf = write(tmp_path, "m.vcf", [GIANT] + FILLER[:38] + E02 + FILLER[38:])
    long_rows, paired = run(tmp_path, vcf, "paired")
    assert [r[7] for r in long_rows] == ["exclude"]
    _, excl = run(tmp_path, vcf, "exclude")
    assert paired == excl                                         # identical discovery
    e02 = [r for r in paired if int(r["Start"]) <= 32250000 <= int(r["End"])]
    assert len(e02) == 1 and e02[0]["IsCxSV"] == "YES" and e02[0]["N_Records"] == "3"


def test_one_sided_partner_is_not_enough(tmp_path):
    near = [rec("chr17", 42_450_000, "svim.INV.n1", "INV", 42_455_000, "011", "++"),
            rec("chr17", 42_460_000, "svim.DEL.n2", "DEL", 42_462_000, "111")]
    vcf = write(tmp_path, "m.vcf", [GIANT] + near)
    long_rows, summ = run(tmp_path, vcf, "paired")
    assert [r[7] for r in long_rows] == ["exclude"]
    loc = [r for r in summ if int(r["Start"]) <= 42_450_000 <= int(r["End"])]
    assert len(loc) == 1 and loc[0]["N_Records"] == "2"           # the long DEL is not a junction here


def test_duplicate_of_the_same_junction_is_no_partner(tmp_path):
    # Sniffles2 BND pair rewritten as INV beside an SVIM INV at the same two ends
    dup = [rec("chr17", 20_000_000, "Sniffles2.INV.x", "INV", 60_000_000, "001", "++"),
           rec("chr17", 20_000_120, "svim.INV.x", "INV", 60_000_090, "010", "++")]
    long_rows, summ = run(tmp_path, write(tmp_path, "m.vcf", dup), "paired")
    assert [r[7] for r in long_rows] == ["exclude", "exclude"] and summ == []


def test_long_partner_must_itself_be_paired(tmp_path):
    # A pairs with B at its start and with a small DEL at its END; B has nothing at its END
    lines = [rec("chr17", 1_000_000, "svim.DUP.A", "DUP", 5_000_000, "011", "-+"),
             rec("chr17", 1_020_000, "svim.DEL.B", "DEL", 20_000_000, "011"),
             rec("chr17", 5_010_000, "svim.DEL.s", "DEL", 5_012_000, "111")]
    long_rows, _ = run(tmp_path, write(tmp_path, "m.vcf", lines), "paired")
    assert dict((r[0], r[7]) for r in long_rows) == {"svim.DUP.A": "exclude", "svim.DEL.B": "exclude"}
    # with a partner at B's END too, both are paired
    lines.append(rec("chr17", 20_030_000, "svim.INV.t", "INV", 20_040_000, "011", "++"))
    long_rows, summ = run(tmp_path, write(tmp_path, "m2.vcf", lines), "paired")
    assert dict((r[0], r[7]) for r in long_rows) == {"svim.DUP.A": "paired", "svim.DEL.B": "paired"}
    assert len(summ) == 1 and int(summ[0]["Size"]) < 200_000      # three parts joined, sized by parts


def test_link_is_a_no_op_without_paired_records(tmp_path):
    raw = ["chr17\t100\t200\tDEL\tr1\tintra\t1", "chr17\t300000\t300002\tINV\tr2_end\tintra\t2"]
    out, n, parts = lpc.link([l + "\n" for l in raw], set())
    assert out == raw and n == 0 and parts == 0


def test_extract_takes_only_paired_long_records(tmp_path):
    summ = tmp_path / "cxsv_summary.tsv"
    summ.write_text("ClusterID\tChr\tStart\tEnd\tSize\tN_Breakpoints\tN_SVtypes\tSV_Types\tHas_Inter\tIsCxSV\n"
                    "1\tchr17\t9999999\t14500001\t60004\t4\t1\tDUP\tNO\tYES\n")
    lines = PAIR + [rec("chr17", 10_010_000, "svim.DEL.u", "DEL", 30_000_000, "010")]   # long, unpaired
    vcf = write(tmp_path, "m.vcf", lines)
    longf = tmp_path / "long.tsv"
    longf.write_text("ID\tCHROM\tPOS\tEND\tSVTYPE\tLENGTH\tSUPP_VEC\tRULE\n"
                     "svim.DUP.a\tchr17\t10000000\t14500000\tDUP\t4500000\t011\tpaired\n"
                     "svim.DUP.b\tchr17\t10030000\t14470000\tDUP\t4440000\t011\tpaired\n"
                     "svim.DEL.u\tchr17\t10010000\t30000000\tDEL\t19990000\t010\texclude\n")
    out = tmp_path / "only.vcf"
    subprocess.run([sys.executable, str(SCRIPTS / "extract_cxsv_vcf.py"), "--vcf", str(vcf), "--summary", str(summ),
                    "--output", str(out), "--span-link-max", str(LINK), "--long-span-rule", "paired",
                    "--long-records", str(longf)], check=True, capture_output=True)
    got = {l.split("\t")[2] for l in out.read_text().splitlines() if l and not l.startswith("#")}
    assert got == {"svim.DUP.a", "svim.DUP.b"}


def test_paired_is_a_configurable_rule():
    assert "paired" in v2b.LONG_RULES
