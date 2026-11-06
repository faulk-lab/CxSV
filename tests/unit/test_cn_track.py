"""Roadmap item 15 Phase 2: the pipeline copy-number track (workflow/lib/cxsv_cn.py,
scripts cn_pon.py / cn_segment.py, the report's copy-number-only section and the
locus-depth boundaries). Synthetic depth only; numpy, no pandas."""
import gzip
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "workflow" / "lib"))
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))
import cxsv_cn as cn       # noqa: E402
import cxsv_cnlink as cnl  # noqa: E402

NB = 4000                  # 1 kb bins per contig (4 Mb)
CONTIGS = ("chr1", "chrX")
PAR = {"chrX": [(0, 200_000)]}


def _fine(rng, cov, sex, events=(), contigs=CONTIGS):
    out = {}
    for c in contigs:
        pl = 2 if (c != "chrX" or sex == "female") else 1
        lam = np.full(NB, cov * pl / 2.0)
        if c == "chrX" and sex == "male":
            lam[:200] = cov                       # PAR: two copies in a male
        for ec, s, e, k in events:
            if ec == c:
                lam[s:e] = cov * k / 2.0
        out[c] = rng.poisson(lam).astype(np.float32)
    return out


def _member(fine, par=PAR):
    coarse = cn.to_coarse(fine)
    gmed = cn.genome_median(coarse)
    sex, _ = cn.infer_sex(coarse, gmed)
    return cn.normalised(coarse, gmed, sex, par), sex


def _pon(rng, n=8, cov=30):
    return cn.Pon.build([_member(_fine(rng, cov, "male" if i % 2 else "female")) for i in range(n)], PAR,
                        meta={"stratum": "T"})


def _write_bins(path, fine):
    with gzip.open(path, "wt") as fh:
        for c, a in fine.items():
            for i, v in enumerate(a):
                fh.write(f"{c}\t{i * 1000}\t{i * 1000 + 1000}\t{v:.2f}\n")


def test_defaults_match_fingerprinted_config():
    import yaml
    c = yaml.safe_load(open(ROOT / "config" / "calling.yaml"))["cn_params"]
    assert c["segmentation"] == cn.DEFAULTS
    assert c["linkage"] == {k: v for k, v in cnl.CN_DEFAULTS.items() if k != "span_link_max"}
    assert c["depth"] == {"bin_bp": 1000, "min_mapq": 20}
    import csv
    with open(ROOT / "calibration" / "library_qc" / "flagged_v1.1.tsv", newline="") as fh:
        flagged = sorted(f"{r['cohort']}/{r['sample_id']}" for r in csv.DictReader(fh, delimiter="\t"))
    assert sorted(c["pon_exclude"]) == flagged and len(flagged) == 24
    par = cn.par_from_params(c["par"])
    assert par["chrX"] == [(0, 2394410), (153925834, 154259566)] and len(par["chrY"]) == 2


def test_infer_sex_and_par_ploidy():
    rng = np.random.default_rng(1)
    for sex in ("male", "female"):
        coarse = cn.to_coarse(_fine(rng, 30, sex))
        assert cn.infer_sex(coarse, cn.genome_median(coarse))[0] == sex
    assert cn.infer_sex(cn.to_coarse(_fine(rng, 30, "male", contigs=("chr1",))), 30.0)[0] == "unknown"
    p = cn.ploidy_track("chrX", 400, "male", PAR)
    assert p[0] == 2 and p[19] == 2 and p[20] == 1
    assert cn.ploidy_track("chrY", 10, "female", PAR).sum() == 0


def test_pon_pools_chrx_and_roundtrips(tmp_path):
    rng = np.random.default_rng(2)
    pon = _pon(rng)
    assert pon.meta["genomes"] == 8 and pon.meta["males"] == 4
    med, rel, n = pon.reference("chrX", 400)
    assert np.nanmedian(med) == pytest.approx(1.0, abs=0.05) and n.min() == 8   # both sexes pooled
    pon.save(tmp_path / "T.npz")
    assert (tmp_path / "T.npz").exists() and not (tmp_path / "T.npz.npz").exists()
    back = cn.Pon.load(tmp_path / "T.npz")
    assert back.meta == pon.meta and set(back.med) == set(pon.med)
    assert np.array_equal(back.n["chr1"], pon.n["chr1"])
    m2, _, n2 = back.reference("chr9", 5)              # absent contig: masked
    assert np.isnan(m2).all() and n2.sum() == 0


def test_segment_genome_dup_trp_male_and_null():
    rng = np.random.default_rng(3)
    pon = _pon(rng)
    assert cn.segment_genome(_fine(rng, 30, "female"), "female", PAR, pon)[0] == []
    ev = [("chrX", 1000, 1300, 2), ("chrX", 1300, 1500, 3), ("chrX", 1500, 1800, 2)]
    segs, st = cn.segment_genome(_fine(rng, 30, "male", ev), "male", PAR, pon)
    assert [s["state"] for s in segs] == [2, 3, 2] and all(s["ploidy"] == 1 for s in segs)
    assert st["masked_fraction"] < 0.05
    inv = dict(id="r", chrom="chrX", pos=1_300_000, chr2="chrX", end=1_500_000, svtype="INV",
               svlen=200_000, supp="", qual=None, inverted=True, line="")
    c = cnl.candidates(segs, [inv])
    assert len(c) == 1 and c[0]["cls"] == "DUP-TRP/INV-DUP" and c[0]["decision"] == "locus"


def test_empty_pon_masks_everything():
    """A stratum without members (e.g. the synthetic chr20 test) gives no segments."""
    rng = np.random.default_rng(4)
    empty = cn.Pon.build([], PAR)
    fine = _fine(rng, 30, "female", [("chr1", 1000, 1500, 3)])
    segs, st = cn.segment_genome(fine, "female", PAR, empty)
    assert segs == [] and (np.isnan(st["masked_fraction"]) or st["masked_fraction"] == 1.0)


def test_pon_and_segment_scripts(tmp_path):
    rng = np.random.default_rng(5)
    rows = ["key\tindividual_id\tstratum\tdepth"]
    for i in range(6):
        p = tmp_path / f"n{i}.bed.gz"
        _write_bins(p, _fine(rng, 30, "male" if i % 2 else "female"))
        rows.append(f"c/n{i}\tI{i}\tT\t{p}")
    (tmp_path / "members.tsv").write_text("\n".join(rows) + "\n")
    (tmp_path / "flag.tsv").write_text("cohort\tsample_id\tgroup\nc\tn5\tT\n")
    params = {"cn_params": {"segmentation": dict(cn.DEFAULTS), "par": {"chrX": [[0, 200000]]},
                            "pon_exclude": ["c/n4"]}}
    (tmp_path / "p.json").write_text(json.dumps(params))
    S = ROOT / "workflow" / "scripts"
    subprocess.run([sys.executable, S / "cn_pon.py", "--members", tmp_path / "members.tsv", "--exclude", tmp_path / "flag.tsv", "--params", tmp_path / "p.json",
                    "--out-dir", tmp_path / "pon", "--threads", "1"], check=True)
    pon = cn.Pon.load(tmp_path / "pon" / "T.npz")
    assert pon.meta["genomes"] == 4                     # n4 (params) and n5 (--exclude) left out
    mem = (tmp_path / "pon" / "pon_members.tsv").read_text()
    assert "c/n5\tI5\tT\tmale" in mem and "no: library-QC flag" in mem
    case = tmp_path / "case.bed.gz"
    _write_bins(case, _fine(rng, 30, "female", [("chr1", 2000, 2300, 3)]))
    subprocess.run([sys.executable, S / "cn_segment.py", "--depth", case, "--pon", tmp_path / "pon" / "T.npz",
                    "--params", tmp_path / "p.json", "--stratum", "T",
                    "--out", tmp_path / "seg.tsv", "--stats", tmp_path / "stats.json"], check=True)
    lines = (tmp_path / "seg.tsv").read_text().splitlines()
    assert len(lines) == 2 and lines[1].startswith("chr1\t") and lines[1].split("\t")[3] == "3"
    st = json.loads((tmp_path / "stats.json").read_text())
    assert st["stratum"] == "T" and st["sex"] == "female" and st["stratum_source"] == "configured"
    # patient mode: PoN directory, stratum from the library-QC group, fallback to the default
    seg = [sys.executable, S / "cn_segment.py", "--depth", case, "--pon", tmp_path / "pon", "--params",
           tmp_path / "p.json", "--stratum", "auto", "--out", tmp_path / "s2.tsv", "--stats", tmp_path / "s2.json"]
    (tmp_path / "lq.tsv").write_text("sample\tqc_group\tqc_group_source\nx\tT\theader\n")
    subprocess.run(seg + ["--library-qc", tmp_path / "lq.tsv", "--default-stratum", "U"], check=True)
    st = json.loads((tmp_path / "s2.json").read_text())
    assert st["stratum"] == "T" and st["stratum_source"] == "library_qc (header)"
    (tmp_path / "lq.tsv").write_text("sample\tqc_group\tqc_group_source\nx\tR10\treads\n")
    (tmp_path / "pon" / "T.npz").rename(tmp_path / "pon" / "D.npz")
    subprocess.run(seg + ["--library-qc", tmp_path / "lq.tsv", "--default-stratum", "D"], check=True)
    st = json.loads((tmp_path / "s2.json").read_text())
    assert st["stratum"] == "D" and st["stratum_source"].startswith("platform default (no PoN")
    r = subprocess.run(seg[:9] + ["Q"] + seg[10:], capture_output=True, text=True)
    assert r.returncode != 0 and "available: ['D']" in r.stderr


def test_report_lists_cn_only_untiered(tmp_path):
    import score_cxsv
    (tmp_path / "st.json").write_text(json.dumps(dict(stratum="ONT_LRSC_R10", stratum_source="library_qc (header)",
                                                      pon_genomes=240, sex="male", x_ratio=0.51,
                                                      masked_fraction=0.05, segments=12)))
    line = score_cxsv.cn_track_line(str(tmp_path / "st.json"))[0]
    assert "ONT_LRSC_R10 (library_qc (header); 240 genomes)" in line and "sex male" in line
    assert score_cxsv.cn_track_line(str(tmp_path / "none.json")) == []
    p = tmp_path / "cn_only.tsv"
    assert score_cxsv.cn_only_report(str(p))[1] == []           # absent file: count only
    p.write_text("candidate\tdecision\tcls\tspan\tsegments\trecord_types\n"
                 "cn1\tcn_only\tCN_linked_gain_loss\t400000\tchr1:1-2:3/2;chr1:2-3:1/2\t\n")
    summ, sec = score_cxsv.cn_only_report(str(p))
    assert "1" in summ[0] and any("not tiered" in l for l in sec)
    assert any(l.startswith("| cn1 | CN_linked_gain_loss | 400000 |") for l in sec)


def test_locus_depth_adds_cn_boundaries(tmp_path):
    """A CN-linked locus whose only junction is an INV still gets its DUP/TRP segments measured."""
    (tmp_path / "ref.fa.fai").write_text("chrX\t5000000\t6\t60\t61\n")
    (tmp_path / "loci.tsv").write_text("Chr\tStart\tEnd\nchrX\t1000000\t1800000\n")
    import cxsv_junctions as cj
    jc = cj.JUNCTION_COLUMNS
    row = dict(chrom1="chrX", pos1=1300000, chrom2="chrX", pos2=1500000, svtype="INV", svlen=200000,
               locus_id="chrX:1000000-1800000", supp_vec="111", dosage=0, n_records=1, record_ids="r1")
    (tmp_path / "j.tsv").write_text("\t".join(jc) + "\n" + "\t".join(str(row[c]) for c in jc) + "\n")
    (tmp_path / "cn_loci.tsv").write_text(
        "locus_id\tcandidate\tcn_class\tclassifier_class\tn_records\tn_long\tn_single_filtered\tsegments\twritten\n"
        "chrX:1000000-1800000\tc1\tDUP-TRP/INV-DUP\tx\t1\t0\t0\t"
        "chrX:1000000-1300000:3/2;chrX:1300000-1500000:4/2;chrX:1500000-1800000:3/2\tyes\n")
    S = ROOT / "workflow" / "scripts" / "locus_segments.py"
    base = [sys.executable, S, "--junctions", tmp_path / "j.tsv", "--loci", tmp_path / "loci.tsv",
            "--fai", tmp_path / "ref.fa.fai"]
    base += ["--flank", "50000", "--gap", "1000", "--margin", "1000"]
    subprocess.run(base + ["--out", tmp_path / "a.bed"], check=True)
    subprocess.run(base + ["--cn-loci", tmp_path / "cn_loci.tsv", "--out", tmp_path / "b.bed"], check=True)
    seg = lambda p: [l.split("\t") for l in p.read_text().splitlines() if "\tsegment\t" in l]
    assert len(seg(tmp_path / "a.bed")) == 1                    # INV ends only
    b = seg(tmp_path / "b.bed")
    assert len(b) == 3 and int(b[0][1]) == 1_001_000 and int(b[-1][2]) == 1_799_000
