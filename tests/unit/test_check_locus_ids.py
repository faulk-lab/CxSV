"""Locus-id rule of bundle validation check 10 (tests/lab/check_locus_ids.py; PI decision 22,
2026-11-13): reference contig, 0 <= start < length, start < end <= length + 1."""
import sys

import pytest

from conftest import REPO

sys.path.insert(0, str(REPO / "tests" / "lab"))
import check_locus_ids as cli  # noqa: E402

L = {"chr1": 1000, "chrM": 16569, "chr7": 160567428}


@pytest.mark.parametrize("lid", [
    "chr1:0-2",                    # start 0
    "chr1:10-1000",                # end = length
    "chr1:10-1001",                # end = length + 1 (unclamped BND window at the last base)
    "chrM:2-16570",                # BOC_3537 (panel-v1.3 test run)
    "chr7:160563865-160567429",    # BOC_1882
    "chr1:999-1001",               # start on the last base, end = length + 1
])
def test_accepted(lid):
    assert cli.problem(lid, L) is None


@pytest.mark.parametrize("lid, why", [
    ("chr1:10-1002", "past the end"),                     # end = length + 2
    ("chr1:10-5000", "past the end"),
    ("chrN:62452837-62454007", "not in the reference"),   # Sniffles2 PAR artefact (item 18)
    ("chrUn_KI270742v1:1-20", "not in the reference"),
    ("chr1:500-500", "start >= end"),                     # zero length
    ("chr1:600-500", "start >= end"),
    ("chr1:1000-1001", "off the end"),                    # starts past the last base
    ("chr1:1001-1002", "off the end"),
    ("chr1:-5-10", "not CONTIG"),                         # negative start
    ("chr1-10-20", "not CONTIG"),
    ("", "not CONTIG"),
])
def test_rejected(lid, why):
    p = cli.problem(lid, L)
    assert p is not None and why in p


def test_main_exit_codes(tmp_path):
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("".join(f"{c}\t{n}\t0\t60\t61\n" for c, n in L.items()))
    t = tmp_path / "s.tsv"
    t.write_text("Rank\tLocus_ID\n1\tchr1:10-1001\n2\tchrM:7-16570\n")
    assert cli.main(["--fai", str(fai), str(t)]) == 0
    t.write_text("Rank\tLocus_ID\n1\tchr1:10-1002\n")
    assert cli.main(["--fai", str(fai), str(t)]) == 1
