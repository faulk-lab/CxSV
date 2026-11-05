"""The single-caller breakend rule of workflow/scripts/vcf_to_bed.py
(--single-caller-bnd-rule; command-line default keep = panel-v1.0 discovery;
pipeline default svim_lowq, config/calling.yaml, adopted 2026-11-05).

Seven LRSC R10 panel genomes carry 3,753-14,258 SVIM-only translocation records at
SVIM's quality floor (QUAL 5); clustered with real SVs they make BND_cluster loci
out of simple deletions (docs/R10_EXCESS.md). The rules leave such records out of
CxSV discovery only: they are not written to the clustering BED and are listed.
"""
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "workflow" / "scripts"))
import vcf_to_bed as v2b  # noqa: E402

HEADER = ("##fileformat=VCFv4.1\n##contig=<ID=chr2,length=242696752>\n##contig=<ID=chr7,length=159345973>\n"
          "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNULL\tSample\tS\n")
GT = "GT:PSV:LN:DR:ST:QV:TY:ID:RAL:AAL:CO"
E = "./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN"


def rec(chrom, pos, vid, svtype, supp, qual, chr2=None, end=None):
    end = end or pos + 1000
    info = (f"SUPP={supp.count('1')};SUPP_VEC={supp};SVLEN={end - pos if svtype != 'TRA' else 0};SVTYPE={svtype};"
            f"SVMETHOD=SURVIVOR1.0.7;CHR2={chr2 or chrom};END={end};CIPOS=0,0;CIEND=0,0;STRANDS=++")
    alt = f"N[{chr2}:{end}[" if svtype == "TRA" else f"<{svtype}>"
    col = f"0/1:NA:0:0,0:++:{qual}:{svtype}:{vid}:N:N:x"
    return "\t".join([chrom, str(pos), vid, "N", alt, str(qual), "PASS", info, GT] + [col if b == "1" else E for b in supp])


RECORDS = [
    rec("chr2", 1000, "svim.BND.1", "TRA", "010", 5, "chr7", 5000),     # SVIM-only, quality floor
    rec("chr2", 2000, "svim.BND.2", "TRA", "010", 13, "chr7", 6000),    # SVIM-only, QUAL 13 (as E07's second junction)
    rec("chr2", 3000, "Sniffles2.BND.3", "TRA", "001", 60, "chr7", 7000),  # Sniffles2-only
    rec("chr2", 4000, "svim.BND.4", "TRA", "110", 30, "chr7", 8000),    # two callers
    rec("chr2", 9000, "svim.DEL.5", "DEL", "010", 5),                    # SVIM-only DEL, quality floor
    rec("chr2", 20000, "cuteSV.DEL.6", "DEL", "111", 40),
]


def ids(tmp_path, rule):
    p = tmp_path / "m.vcf"
    p.write_text(HEADER + "".join(r + "\n" for r in RECORDS))
    rows = []
    bed = list(v2b.process_vcf(str(p), 1, 500_000, {}, "exclude", [], rule, 10, rows))
    return {l.split("\t")[4].replace("_dest", "") for l in bed}, {r[0] for r in rows}


ALL = {"svim.BND.1", "svim.BND.2", "Sniffles2.BND.3", "svim.BND.4", "svim.DEL.5", "cuteSV.DEL.6"}


@pytest.mark.parametrize("rule,left_out", [
    ("keep", set()),
    ("svim_lowq", {"svim.BND.1"}),
    ("svim_lowq_any", {"svim.BND.1", "svim.DEL.5"}),
    ("single_bnd", {"svim.BND.1", "svim.BND.2", "Sniffles2.BND.3"}),
])
def test_rules(tmp_path, rule, left_out):
    kept, listed = ids(tmp_path, rule)
    assert listed == left_out
    assert kept == ALL - left_out


def test_keep_is_main(tmp_path):
    """The default writes exactly what the rule-less call writes (main's discovery)."""
    p = tmp_path / "m.vcf"
    p.write_text(HEADER + "".join(r + "\n" for r in RECORDS))
    a = list(v2b.process_vcf(str(p), 1, 500_000, {}, "exclude", []))
    b = list(v2b.process_vcf(str(p), 1, 500_000, {}, "exclude", [], "keep", 10, []))
    assert a == b


def test_unknown_rule():
    with pytest.raises(ValueError):
        v2b.single_caller_excluded("bogus", "TRA", "010", 5, 10)


def test_calling_yaml_default_is_svim_lowq():
    """The pipeline default (config/calling.yaml) is svim_lowq at QUAL 10 (docs/R10_EXCESS.md)."""
    import yaml
    cp = yaml.safe_load((REPO / "config" / "calling.yaml").read_text())["cxsv_params"]
    assert cp["single_caller_bnd_rule"] == "svim_lowq"
    assert cp["svim_min_qual_discovery"] == 10
    assert cp["single_caller_bnd_rule"] in v2b.SINGLE_BND_RULES
