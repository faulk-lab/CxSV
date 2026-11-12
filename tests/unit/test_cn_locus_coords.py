"""Copy-number locus coordinates (panel-v1.3 item g; docs/PLAN_panel_v1.3.md).

panel-v1.2 took a CN locus's Chr:Start-End from its attached breakends alone and wrote one
locus per contig of those breakends. BOC_5883's 17p13.3 loss (segments chr17:10000-2944000
and 3392000-5037000 at 1/2), whose only junction was a translocation end at chr17:1, became
the Tier 2 locus chr17:0-2; the translocation's mate (Sniffles2 contig 'chrN') and the mates
of two other candidates became 2 bp loci on contigs without a copy-number change. These
tests pin the segment-anchored construction and the check that rejects malformed loci.
"""
import csv
import json
import subprocess
import sys

import pytest
import yaml

from conftest import REPO
import cxsv_cnlink as cl  # noqa: E402

LOSS17 = [dict(chrom="chr17", start=10_000, end=2_944_000, state=1, ploidy=2),
          dict(chrom="chr17", start=3_392_000, end=5_037_000, state=1, ploidy=2)]


def _cand(segs, i=290):
    return dict(id=i, segs=segs, records=[])


def test_check_locus_rejects_zero_length_and_out_of_segment():
    c = _cand(LOSS17)
    for s, e in ((5, 5), (7, 3), (-1, 100)):
        with pytest.raises(ValueError, match="zero or negative length"):
            cl.check_locus("chr17", s, e, c)
    with pytest.raises(ValueError, match="overlaps none"):
        cl.check_locus("chr17", 0, 2, c)                     # the v1.2 BOC_5883 locus: before the first segment
    with pytest.raises(ValueError, match="overlaps none"):
        cl.check_locus("chr17", 2_944_000, 3_392_000, c)     # the gap between the segments
    with pytest.raises(ValueError, match="overlaps none"):
        cl.check_locus("chrN", 62_454_005, 62_454_007, c)    # a mate contig without a segment
    with pytest.raises(ValueError, match="past the contig end"):
        cl.check_locus("chr17", 0, 90_000_000, c, length=84_276_897)
    assert cl.check_locus("chr17", 0, 5_037_000, c, length=84_276_897) == (0, 5_037_000)
    assert cl.check_locus("chr17", 2_943_999, 2_944_001, c) == (2_943_999, 2_944_001)


def test_locus_span_is_segments_plus_breakends():
    c = _cand(LOSS17)
    assert cl.locus_span(c, "chr17", [(0, 2)]) == (0, 5_037_000)
    assert cl.locus_span(c, "chr17", []) == (10_000, 5_037_000)
    gain = _cand([dict(chrom="chr1", start=238_516_000, end=238_599_000, state=3, ploidy=2),
                  dict(chrom="chr1", start=238_599_000, end=238_814_000, state=4, ploidy=2)])
    assert cl.locus_span(gain, "chr1", [(238_599_733, 238_599_735), (238_814_324, 238_814_326)]) == (238_516_000, 238_814_326)
    with pytest.raises(ValueError, match="no copy-number segment on chr5"):
        cl.locus_span(c, "chr5", [(1, 3)])
    merged = _cand(LOSS17 + [dict(chrom="chr5", start=1_000_000, end=1_200_000, state=3, ploidy=2)])
    assert cl.segment_contigs(merged) == ["chr17", "chr5"] and cl.segment_contigs(c) == ["chr17"]


TRA = ("chr17\t1\tSniffles2.BND.8FBS11\tN\tN[chr5:62454006[\t60\tPASS\tSUPP=1;SUPP_VEC=001;SVLEN=0;SVTYPE=TRA;"
       "SVMETHOD=SURVIVOR1.0.7;CHR2=chr5;END=62454006;CIPOS=0,0;CIEND=0,0;STRANDS=++\tGT\t./.\t./.\t0/1\n")
VCF_HEAD = ("##fileformat=VCFv4.2\n##contig=<ID=chr5,length=181265378>\n##contig=<ID=chr17,length=84276897>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSniffles2\tcuteSV\tSVIM\n")


def _run_linkage(tmp_path, summary_rows=()):
    import classify_cxsv
    cfg = yaml.safe_load(open(REPO / "config" / "calling.yaml"))
    (tmp_path / "params.json").write_text(json.dumps({k: cfg[k] for k in ("cn_params", "cxsv_params", "junction_params")}))
    (tmp_path / "ref.fa.fai").write_text("chr5\t181265378\t6\t60\t61\nchr17\t84276897\t200000000\t60\t61\n")
    (tmp_path / "merged_sv.vcf").write_text(VCF_HEAD + TRA)
    (tmp_path / "junction_only.vcf").write_text(VCF_HEAD)
    (tmp_path / "junction_clusters.bed").write_text("")
    with open(tmp_path / "junction_summary.tsv", "w") as fh:
        fh.write("\t".join(classify_cxsv.HEADER) + "\n")
        for r in summary_rows:
            fh.write("\t".join(str(r.get(k, "NA")) for k in classify_cxsv.HEADER) + "\n")
    cols = ["chrom", "start", "end", "state", "ploidy", "cn", "mean_ratio", "n_bins", "z", "nested",
            "start_mask_before", "start_mask_after", "end_mask_before", "end_mask_after"]
    with open(tmp_path / "segments.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for s in LOSS17:
            mb = 10_000 if s["start"] == 10_000 else 0                 # telomere masked before the first segment
            fh.write("\t".join(map(str, [s["chrom"], s["start"], s["end"], 1, 2, 1.0, 0.5, 100, -20, False,
                                         mb, 0, 0, 0])) + "\n")
    (tmp_path / "stats.json").write_text(json.dumps(dict(sex="female", sex_call="female")))
    out = {k: tmp_path / v for k, v in dict(clusters="clusters.bed", summary="cxsv_summary.tsv", vcf="cxsv_only.vcf",
                                             candidates="cn_candidates.tsv", loci="cn_loci.tsv", cn_only="cn_only.tsv").items()}
    cmd = [sys.executable, str(REPO / "workflow" / "scripts" / "cn_linkage.py"), "--vcf", tmp_path / "merged_sv.vcf",
           "--clusters", tmp_path / "junction_clusters.bed", "--summary", tmp_path / "junction_summary.tsv",
           "--cxsv-vcf", tmp_path / "junction_only.vcf", "--single", tmp_path / "none.tsv",
           "--segments", tmp_path / "segments.tsv", "--stats", tmp_path / "stats.json", "--fai", tmp_path / "ref.fa.fai",
           "--calling", tmp_path / "params.json", "--scripts", REPO / "workflow" / "scripts", "--work", tmp_path / "work",
           "--out-clusters", out["clusters"], "--out-summary", out["summary"], "--out-vcf", out["vcf"],
           "--out-candidates", out["candidates"], "--out-loci", out["loci"], "--out-cn-only", out["cn_only"]]
    r = subprocess.run([str(x) for x in cmd], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    rd = lambda p: list(csv.DictReader(open(p), delimiter="\t"))   # noqa: E731
    return r, rd(out["loci"]), rd(out["summary"]), open(out["clusters"]).read().splitlines()


def test_locus_span_clamped_to_contig_end():
    # HG03667 (panel-v1.2-r2): a gain reaching the last, partial 1 kb bin of chr2
    gain = _cand([dict(chrom="chr2", start=242_600_000, end=242_697_000, state=3, ploidy=2)])
    rows = [(242_559_540, 242_559_542)]
    assert cl.locus_span(gain, "chr2", rows) == (242_559_540, 242_697_000)
    s, e = cl.locus_span(gain, "chr2", rows, 242_696_752)
    assert (s, e) == (242_559_540, 242_696_752)
    assert cl.check_locus("chr2", s, e, gain, 242_696_752) == (s, e)
    with pytest.raises(ValueError):
        cl.check_locus("chr2", 242_559_540, 242_697_000, gain, 242_696_752)


def test_boc5883_loss_gets_segment_anchored_locus(tmp_path):
    r, loci, summ, bed = _run_linkage(tmp_path)
    assert [(L["locus_id"], L["junction_span"], L["written"]) for L in loci] == [("chr17:0-5037000", "chr17:0-2", "yes")]
    row = [x for x in summ if x["Chr"] == "chr17"]
    assert len(summ) == 1 and len(row) == 1
    assert (row[0]["Start"], row[0]["End"], row[0]["Size"]) == ("0", "5037000", "5037000")
    assert row[0]["Cmplx_Class"] == "CN_linked_loss" and row[0]["CxSV_Criteria"].startswith("CN_linked")
    for L in loci + [dict(locus_id=f"{x['Chr']}:{x['Start']}-{x['End']}") for x in summ]:   # every locus well-formed
        chrom, rng = L["locus_id"].split(":")
        s, e = map(int, rng.split("-"))
        cl.check_locus(chrom, s, e, _cand(LOSS17))
    assert all(l.split("\t")[0] == "chr17" for l in bed)            # no locus row on the mate contig
    assert "breakends on chr5 (no copy-number segment there)" in r.stderr
    assert "Sniffles2.BND.8FBS11" in open(tmp_path / "cxsv_only.vcf").read()   # the junction stays with the locus


def test_cn_locus_with_junction_track_span_is_not_written_again(tmp_path):
    jt = dict(ClusterID=7, Chr="chr17", Start=0, End=2, Size=2, IsCxSV="YES", Cmplx_Class="BND_cluster")
    r, loci, summ, bed = _run_linkage(tmp_path, [jt])
    assert len(loci) == 1 and loci[0]["written"].startswith("no:")
    assert [(x["Chr"], x["Start"], x["End"]) for x in summ] == [("chr17", "0", "2")]   # the junction-track row only
