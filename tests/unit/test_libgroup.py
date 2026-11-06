"""Library-QC calibration group detection (lib/cxsv_libgroup.py, docs/LIBRARY_QC.md 4b).
Headers are abridged copies of real ones from each panel platform and of the patient
stand-ins (2026-11-06)."""
import sys

from conftest import REPO

sys.path.insert(0, str(REPO / "workflow" / "lib"))
import cxsv_libgroup as lg  # noqa: E402

H_1019 = {"RG": [{"ID": "ERR15139516", "SM": "NA07019", "PL": "ONT", "PU": "ERR15139516",
                  "DS": "PRJEB89727_SRA_mirror"}],
          "PG": [{"ID": "minimap2", "PN": "minimap2", "VN": "2.24-r1122",
                  "CL": "minimap2 -t 24 -x map-ont -a --eqx -K 10g -R @RG\\tID:ERR15139516\\tPL:ONT "
                        "/home/minknow/refs/hs1/hs1.map-ont.k17.mmi -"}]}
H_LRSC_PMDV = {"PG": [{"ID": "minimap2", "PN": "minimap2", "VN": "2.24-r1122",
                       "CL": "minimap2 -t 16 -y -x map-ont -a --eqx -k 17 -K 10g /cromwell-executions/"
                             "cardEndToEndVcfMethyl/5ef812cb-6195-418d-bb39-b942468d3b5d/call-mm_align/"
                             "inputs/1438353101/T2T-CHM13v2.fasta -"}]}
H_LRSC_R9 = {"PG": [{"ID": "minimap2", "PN": "minimap2", "VN": "2.28-r1209",
                     "CL": "minimap2 -t 89 -y -a /n/dat/chm13/T2T-CHM13v2_lrhq.mmi "
                           "GM18507-ONT-WGS_R9-LSK110.dorado050.sup.5mCG.33.fastq.gz"},
                    {"ID": "longphase", "PN": "longphase", "CL": "longphase haplotag --snp-file=GM18507-ONT-WGS_"
                     "R9-LSK110.dorado050.sup.5mCG.33.clair3.phased.vcf.gz"}]}
H_PACBIO = {"RG": [{"ID": "16b2c9d5/26--26", "PL": "PACBIO",
                    "DS": "READTYPE=CCS;Ipd:CodecV1=ip;PulseWidth:CodecV1=pw;BINDINGKIT=103-496-900;"
                          "BASECALLERVERSION=5.0"}],
            "PG": [{"ID": "ccs", "PN": "ccs", "CL": "/opt/pacbio/bin/ccs --streamed"}]}
H_PC = {"RG": [{"ID": "3559", "SM": "pos_control_3559", "PL": "ONT",
                "DS": "dorado 0.8.3 dna_r10.4.1_e8.2_400bps_sup@v5.0.0 (5mC_5hmC,6mA); realigned with "
                      "minimap2 2.24-r1122 k17 to hs1"}]}


def hdr(h):
    return lg.classify_text(lg.header_text(h), "header")


def test_header_each_platform():
    assert hdr(H_PACBIO)[0] == lg.PACBIO
    assert hdr(H_PC)[0] == lg.R10 and "r10.4.1" in hdr(H_PC)[1]
    assert hdr(H_LRSC_R9)[0] == lg.R9 and "@PG CL" in hdr(H_LRSC_R9)[1]
    # no chemistry in the header: these two are decided by their reads
    assert hdr(H_1019) == (None, "header: no chemistry cue (PL ONT)")
    assert hdr(H_LRSC_PMDV) == (None, "header: no chemistry cue (no @RG)")


def test_header_no_false_hits_and_ambiguity():
    # chr9 / chr10 / -K 10g / hex paths must not read as R9 / R10
    h = {"PG": [{"ID": "x", "CL": "samtools view -o out.cram chr9:1-2 chr10 -K 10g /a/9e10ff/r10x"}]}
    assert hdr(h)[0] is None
    both = {"RG": [{"ID": "a", "PL": "ONT", "DS": "dna_r9.4.1_e8"}, {"ID": "b", "PL": "ONT", "DS": "dna_r10.4.1"}]}
    g, how = hdr(both)
    assert g is None and "ambiguous" in how
    assert hdr({"RG": [{"ID": "run1_dna_r10.4.1_e8.2_400bps_sup@v4.3.0", "PL": "ONT"}]})[0] == lg.R10
    assert hdr({"CO": ["library SQK-LSK109"]})[0] == lg.R9


def cues(n, q, name="0b1e7c2a-1111-2222-3333-444455556666", rg=None, rq=False):
    c = lg.ReadCues()
    for _ in range(n):
        c.add(name, q, rg, rq)
    return c.summary()


def test_reads_each_platform():
    rules = {}
    assert lg.classify_reads(cues(500, [33] * 50, "ERR15139516.17"), rules)[0] == lg.ONT1019
    assert lg.classify_reads(cues(500, list(range(30, 56))), rules)[0] == lg.R10        # mean 42.5
    assert lg.classify_reads(cues(500, list(range(12, 36))), rules)[0] == lg.R9         # mean 23.5
    assert lg.classify_reads(cues(500, list(range(20, 60)), "m84055_251224_081058_s4/1234/ccs"), rules)[0] == lg.PACBIO
    assert lg.classify_reads(cues(500, list(range(20, 60)), rq=True), rules)[0] == lg.PACBIO
    g, how = lg.classify_reads(cues(500, list(range(17, 42))), rules)                   # mean 29
    assert g is None and "between" in how
    assert lg.classify_reads(cues(50, [40] * 9 + [41]), rules)[0] is None               # too few reads
    # a dorado read-group tag naming the model decides before the qualities
    s = cues(500, list(range(12, 36)), rg="abc_dna_r10.4.1_e8.2_400bps_hac@v4.1.0")
    assert lg.classify_reads(s, rules)[0] == lg.R10


def test_resolve_order_and_evidence():
    dflt = {"ONT": lg.ONT1019, "PACBIO": lg.PACBIO}
    hdr_r10, rd_r9, none = (lg.R10, "header: x"), (lg.R9, "reads: y"), (None, "none")
    r = lg.resolve("ONT_LRSC_R9", hdr_r10, none, "ONT", dflt)
    assert r["qc_group"] == "ONT_LRSC_R9" and r["qc_group_source"] == "sample sheet qc_group"
    assert "auto-detection would give ONT_LRSC_R10" in r["qc_group_evidence"]
    r = lg.resolve(None, hdr_r10, rd_r9, "ONT", dflt)
    assert r["qc_group"] == lg.R10 and r["qc_group_source"] == "header" and "disagree" in r["qc_group_evidence"]
    r = lg.resolve(None, none, rd_r9, "ONT", dflt)
    assert (r["qc_group"], r["qc_group_source"]) == (lg.R9, "reads")
    r = lg.resolve(None, none, none, "ONT", dflt)
    assert (r["qc_group"], r["qc_group_source"]) == (lg.ONT1019, "platform default (ONT)")
    r = lg.resolve(None, (lg.PACBIO, "header: PL"), none, "ONT", dflt, groups={lg.PACBIO: {}})
    assert "sample sheet says platform ONT" in r["qc_group_evidence"]
    r = lg.resolve("MyGroup", none, none, "ONT", dflt, groups={lg.R9: {}})
    assert "group MyGroup has no calibration" in r["qc_group_evidence"]


def test_config_detection_block():
    import yaml
    cfg = yaml.safe_load((REPO / "config" / "library_qc.yaml").read_text())["library_qc"]
    d = cfg["group_detection"]
    assert set(d) == set(lg.DEFAULTS)
    assert d["r9_max_median_q"] < d["r10_min_median_q"]


def test_library_qc_py_end_to_end(tmp_path):
    """library_qc.py on a tiny BAM: detection from the header and from the reads."""
    import json
    import subprocess
    pysam = __import__("pytest").importorskip("pysam")
    ref = tmp_path / "ref.fa"
    ref.write_text(">chr1\n" + "ACGT" * 2500 + "\n")
    pysam.faidx(str(ref))
    win = tmp_path / "w.bed"
    win.write_text("chr1\t0\t10000\tw1\n")
    vcf = tmp_path / "m.vcf"
    vcf.write_text("##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    params = tmp_path / "p.json"
    params.write_text(json.dumps({"k": 6, "eps": 0.001, "default_group": {"ONT": lg.ONT1019},
                                  "group_detection": {"min_reads": 10}, "groups": {}}))
    script = REPO / "workflow" / "scripts" / "library_qc.py"

    def run(header, q, override=""):
        bam = tmp_path / "a.bam"
        with pysam.AlignmentFile(str(bam), "wb", header={"HD": {"VN": "1.6", "SO": "coordinate"},
                                                       "SQ": [{"SN": "chr1", "LN": 10000}], **header}) as out:
            for i in range(50):
                r = pysam.AlignedSegment()
                r.query_name, r.query_sequence = f"0b1e7c2a-1111-2222-3333-{i:012d}", "ACGT" * 25
                r.flag, r.reference_id, r.reference_start, r.mapping_quality = 0, 0, i * 100, 60
                r.cigartuples, r.query_qualities = [(0, 100)], pysam.qualitystring_to_array(chr(33 + q) * 100) \
                    if isinstance(q, int) else q
                out.write(r)
        pysam.index(str(bam))
        o = tmp_path / "o.tsv"
        subprocess.run([sys.executable, str(script), "--aln", str(bam), "--ref", str(ref), "--windows", str(win),
                        "--merged-vcf", str(vcf), "--out", str(o), "--thresholds", str(params),
                        "--group-override", override, "--platform", "ONT"], check=True, capture_output=True)
        lines = o.read_text().splitlines()
        return dict(zip(lines[0].split("\t"), lines[1].split("\t")))

    r = run({"RG": [{"ID": "3559", "PL": "ONT", "DS": "dna_r10.4.1_e8.2_400bps_sup@v5.0.0"}]}, 20)
    assert (r["qc_group"], r["qc_group_source"]) == (lg.R10, "header")
    r = run({}, 33)                                   # uniform qualities, no header cue
    assert (r["qc_group"], r["qc_group_source"]) == (lg.ONT1019, "reads")
    r = run({}, pysam.qualitystring_to_array("".join(chr(33 + 12 + i % 24) for i in range(100))))
    assert (r["qc_group"], r["qc_group_source"]) == (lg.R9, "reads")
    r = run({}, 33, override="ONT_LRSC_R10")
    assert (r["qc_group"], r["qc_group_source"]) == ("ONT_LRSC_R10", "sample sheet qc_group")
