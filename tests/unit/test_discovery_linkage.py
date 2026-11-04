"""Roadmap item 13: records longer than cxsv_params.span_link_max are clustered at their
two breakends (workflow/scripts/vcf_to_bed.py), so a multi-Mb record no longer chains
every SV it spans into one oversized cluster that classify_cxsv.py drops whole.

The cluster step is `bedtools cluster -d D` in the workflow; the unit job has no bedtools,
so `cluster()` below re-implements it (the same rule, verified against the released
clusters.bed of all 1,619 panel-v1.0 genomes in calibration/oversized_clusters).
Coordinates of the chaining case are those of the full-genome positive control: E02
(NF1 delINVdel, chr17:32.24-32.27 Mb) inside a 25.7 Mb Sniffles2-only DEL
chr17:16,716,382-42,430,144 that recurs in 645 panel genomes.
"""
import csv
import random
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))
sys.path.insert(0, str(REPO / "workflow" / "lib"))
import cxsv_junctions as cj  # noqa: E402
import vcf_to_bed as v2b      # noqa: E402

LINK = 500_000
HEADER = ("##fileformat=VCFv4.1\n##contig=<ID=chr17,length=84276897>\n"
          "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNULL\tSample\tS\n")
GT = "GT:PSV:LN:DR:ST:QV:TY:ID:RAL:AAL:CO"


def rec(chrom, pos, vid, svtype, end, supp="001", chr2=None, alt=None):
    svlen = end - pos if svtype != "DEL" else -(end - pos)
    col = f"0/1:NA:{abs(svlen)}:0,0:+-:20:{svtype}:{vid}:N:N:{chrom}_{pos}-{chrom}_{end}"
    empty = "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN"
    cols = [col if b == "1" else empty for b in supp]
    info = (f"SUPP={supp.count('1')};SUPP_VEC={supp};SVLEN={svlen};SVTYPE={svtype};"
            f"SVMETHOD=SURVIVOR1.0.7;CHR2={chr2 or chrom};END={end};CIPOS=0,0;CIEND=0,0;STRANDS=+-")
    return "\t".join([chrom, str(pos), vid, "N", alt or f"<{svtype}>", "20", "PASS", info, GT] + cols)


def write(tmp, name, lines):
    p = tmp / name
    p.write_text(HEADER + "".join(l + "\n" for l in lines))
    return p


def cluster(bed_lines, d=50_000):
    """bedtools cluster -d d on sorted BED rows: rows + cluster id (7th column)."""
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


def classify(tmp, vcf, link):
    bed = list(v2b.process_vcf(str(vcf), 1, link))
    cl = tmp / f"clusters.{link}.bed"
    cl.write_text("".join(l + "\n" for l in cluster(bed)))
    fai = tmp / "ref.fai"
    fai.write_text("chr17\t84276897\t7\t60\t61\n")
    out = tmp / f"summary.{link}.tsv"
    subprocess.run([sys.executable, str(SCRIPTS / "classify_cxsv.py"), "--bed", str(cl), "--vcf", str(vcf),
                    "--fai", str(fai), "--output", str(out), "--max-cluster-size", "500000"], check=True,
                   capture_output=True)
    with open(out, newline="") as fh:
        return bed, list(csv.DictReader(fh, delimiter="\t"))


# E02-like delINVdel (two deletions flanking an inversion), the 25.7 Mb artefact DEL,
# and unrelated small deletions spread across the artefact's interval
E02 = [rec("chr17", 32243849, "Sniffles2.DEL.e1", "DEL", 32247849, "111"),
       rec("chr17", 32247849, "Sniffles2.INV.e2", "INV", 32267849, "011"),
       rec("chr17", 32267849, "Sniffles2.DEL.e3", "DEL", 32271849, "111")]
GIANT = rec("chr17", 16716382, "Sniffles2.DEL.giant", "DEL", 42430144, "001")
FILLER = [rec("chr17", p, f"cuteSV.DEL.f{p}", "DEL", p + 300, "100") for p in range(17_000_000, 42_000_000, 400_000)
          if abs(p - 32_250_000) > 150_000]


def test_long_span_is_written_as_two_breakends(tmp_path):
    vcf = write(tmp_path, "m.vcf", [GIANT, rec("chr17", 1_000_000, "svim.DEL.s", "DEL", 1_100_000, "010")])
    rows = [l.split("\t") for l in v2b.process_vcf(str(vcf), 1, LINK)]
    giant = [r for r in rows if r[4].startswith("Sniffles2.DEL.giant")]
    assert [(r[1], r[2], r[4]) for r in giant] == [("16716381", "16716383", "Sniffles2.DEL.giant"),
                                                   ("42430143", "42430145", "Sniffles2.DEL.giant_end")]
    small = [r for r in rows if r[4] == "svim.DEL.s"]
    assert [(r[1], r[2]) for r in small] == [("1000000", "1100000")]      # 100 kb: still one interval
    # 0 restores the interval rule for every record
    rows0 = [l.split("\t") for l in v2b.process_vcf(str(vcf), 1, 0)]
    assert [(r[1], r[2]) for r in rows0 if r[4].startswith("Sniffles2.DEL.giant")] == [("16716382", "42430144")]


def test_canonical_id_strips_the_end_suffix():
    assert cj.canonical_id("Sniffles2.DEL.giant_end") == "Sniffles2.DEL.giant"
    assert cj.canonical_id("svim.BND.7_mate") == "svim.BND.7"
    assert cj.canonical_id("svim.DEL.3") == "svim.DEL.3"


def test_planted_event_chained_to_a_long_record_is_recovered(tmp_path):
    vcf = write(tmp_path, "m.vcf", [GIANT] + FILLER[:38] + E02 + FILLER[38:])
    # before item 13: one 25.7 Mb cluster, dropped whole
    _, before = classify(tmp_path, vcf, 0)
    hit = [r for r in before if int(r["Start"]) <= 32243849 <= int(r["End"])]
    assert len(hit) == 1 and hit[0]["CxSV_Criteria"] == "FILTERED_oversized_cluster"
    assert int(hit[0]["Size"]) > 25_000_000 and not any(r["IsCxSV"] == "YES" for r in before)
    # after: the planted event is its own CxSV locus, without the long record
    _, after = classify(tmp_path, vcf, LINK)
    hit = [r for r in after if int(r["Start"]) <= 32250000 <= int(r["End"])]
    assert len(hit) == 1 and hit[0]["IsCxSV"] == "YES"
    assert (int(hit[0]["Start"]), int(hit[0]["End"])) == (32243849, 32271849)
    assert hit[0]["N_Records"] == "3" and hit[0]["SV_Types"] == "DEL,INV"
    assert not any(r["CxSV_Criteria"] == "FILTERED_oversized_cluster" for r in after)
    # the long record's own breakend clusters hold one junction each: not CxSVs
    ends = [r for r in after if int(r["Start"]) in (16716381, 42430143)]
    assert len(ends) == 2 and all(r["IsCxSV"] == "NO" and r["N_Records"] == "1" for r in ends)


def test_long_record_joins_the_locus_at_its_breakend(tmp_path):
    # an inversion and a deletion 20-30 kb from the long DEL's END: the DEL is a
    # junction of that locus (like a translocation partner), not of what it spans
    near = [rec("chr17", 42450000, "svim.INV.n1", "INV", 42455000, "011"),
            rec("chr17", 42460000, "svim.DEL.n2", "DEL", 42462000, "111")]
    vcf = write(tmp_path, "m.vcf", [GIANT] + E02 + near)
    _, after = classify(tmp_path, vcf, LINK)
    loc = [r for r in after if int(r["Start"]) <= 42450000 <= int(r["End"])]
    assert len(loc) == 1 and loc[0]["N_Records"] == "3" and loc[0]["IsCxSV"] == "YES"
    assert int(loc[0]["Size"]) < 50_000                           # span of the local rows only
    e02 = [r for r in after if int(r["Start"]) <= 32250000 <= int(r["End"])]
    assert e02 and e02[0]["N_Records"] == "3"                     # E02 does not get the long DEL


def test_clusters_that_are_not_oversized_are_unchanged():
    # with span_link_max == max_cluster_size, a record that changes representation was
    # already in an oversized cluster, so every other cluster keeps its exact membership
    rnd = random.Random(13)
    lines = []
    for i in range(3000):
        s = rnd.randrange(1, 80_000_000)
        L = rnd.choice([60, 300, 2_000, 30_000, 200_000]) if rnd.random() > 0.004 else rnd.randrange(600_000, 9_000_000)
        lines.append(f"chr17\t{s}\t{s + L}\tDEL\tr{i}\tintra")

    def split(ls):
        out = []
        for l in ls:
            f = l.split("\t")
            s, e = int(f[1]), int(f[2])
            if e - s > LINK:
                out += [f"chr17\t{s - 1}\t{s + 1}\tDEL\t{f[4]}\tintra", f"chr17\t{e - 1}\t{e + 1}\tDEL\t{f[4]}_end\tintra"]
            else:
                out.append(l)
        return out

    def ok_clusters(ls):
        by = {}
        for l in cluster(ls):
            f = l.split("\t")
            by.setdefault(f[6], []).append((int(f[1]), int(f[2]), f[4]))
        return {frozenset(x[2] for x in v) for v in by.values()
                if max(x[1] for x in v) - min(x[0] for x in v) <= LINK}

    before, after = ok_clusters(lines), ok_clusters(split(lines))
    assert before and before <= after                      # nothing that was classified changes
    assert len(after) > len(before)                        # and the oversized ones are now classified


def test_extract_cxsv_vcf_uses_breakends_of_long_and_interchromosomal_records(tmp_path):
    summ = tmp_path / "cxsv_summary.tsv"
    summ.write_text("ClusterID\tChr\tStart\tEnd\tSize\tN_Breakpoints\tN_SVtypes\tSV_Types\tHas_Inter\tIsCxSV\n"
                    "1\tchr17\t32243849\t32271849\t28000\t4\t2\tDEL,INV\tNO\tYES\n"
                    "2\tchr17\t50000000\t50010000\t10000\t4\t2\tDEL,INV\tNO\tYES\n")
    lines = E02 + [GIANT,                                                     # spans locus 1, ends outside
                   rec("chr17", 32000000, "svim.DUP.edge", "DUP", 50005000, "010"),  # END breakend in locus 2
                   rec("chr17", 32100000, "svim.DEL.mid", "DEL", 32250000, "010"),   # 150 kb, overlaps: kept
                   rec("chr17", 30000000, "svim.TRA.x", "TRA", 40000000, "010", chr2="chr5")]  # END on chr5
    vcf = write(tmp_path, "m.vcf", lines)

    def ids(link):
        out = tmp_path / f"only.{link}.vcf"
        subprocess.run([sys.executable, str(SCRIPTS / "extract_cxsv_vcf.py"), "--vcf", str(vcf), "--summary",
                        str(summ), "--output", str(out), "--span-link-max", str(link)], check=True,
                       capture_output=True)
        return {l.split("\t")[2] for l in out.read_text().splitlines() if l and not l.startswith("#")}

    e02 = {"Sniffles2.DEL.e1", "Sniffles2.INV.e2", "Sniffles2.DEL.e3"}
    assert ids(LINK) == e02 | {"svim.DUP.edge", "svim.DEL.mid"}
    assert ids(0) == e02 | {"Sniffles2.DEL.giant", "svim.DUP.edge", "svim.DEL.mid", "svim.TRA.x"}


@pytest.mark.parametrize("key", ["span_link_max", "max_cluster_size"])
def test_calling_config_links_long_records_at_max_cluster_size(key):
    import yaml
    cfg = yaml.safe_load((REPO / "config" / "calling.yaml").read_text())["cxsv_params"]
    assert cfg[key] == 500_000 and cfg["span_link_max"] == cfg["max_cluster_size"]
