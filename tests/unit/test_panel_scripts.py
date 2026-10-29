"""Panel aggregation, junction sites and patient labelling scripts (end to end on tiny inputs)."""
import csv
import subprocess
import sys

from conftest import REPO

SCRIPTS = REPO / "workflow" / "scripts"
HEADER = ["Chr", "Start", "End", "Size", "SV_Types", "Cmplx_Class", "CxSV_Criteria", "Has_Inter_Chromosomal"]


def master_table(path, loci):
    with open(path, "w") as fh:
        fh.write("\t".join(HEADER) + "\n")
        for chrom, s, e in loci:
            fh.write(f"{chrom}\t{s}\t{e}\t{e - s}\tDEL,INV\tDEL+INV\tC1\tNO\n")
    return path


def read(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def build_panel(tmp_path, samples):
    man = tmp_path / "manifest.tsv"
    with open(man, "w") as fh:
        fh.write("cohort\tsample_id\tindividual_id\tplatform\tcoverage_x\tmaster_table\n")
        for cohort, sid, ind, plat, loci in samples:
            t = master_table(tmp_path / f"{cohort}_{sid}.tsv", loci)
            fh.write(f"{cohort}\t{sid}\t{ind}\t{plat}\t30\t{t}\n")
    out = {k: tmp_path / f"{k}.tsv" for k in ("calls", "loci", "sharing", "summary")}
    subprocess.run([sys.executable, SCRIPTS / "build_panel_loci.py", "--manifest", man,
                    *[x for k, v in out.items() for x in (f"--{k}", v)]], check=True)
    return out


def test_carriers_are_distinct_individuals(tmp_path):
    out = build_panel(tmp_path, [
        ("cohA", "HG001", "HG001", "ONT", [("chr1", 1000, 5000)]),
        ("cohB", "HG001", "HG001", "PACBIO", [("chr1", 1200, 5200)]),   # same person, 2nd cohort
        ("cohA", "HG002", "HG002", "ONT", [("chr1", 1100, 4000), ("chr2", 10, 20)]),
    ])
    loci = read(out["loci"])
    chr1 = [l for l in loci if l["Chr"] == "chr1"][0]
    assert chr1["N_Individuals"] == "2" and chr1["N_Calls"] == "3"
    assert chr1["Individuals_By_Platform"] == "ONT:2;PACBIO:1"
    total = [r for r in read(out["summary"]) if r["level"] == "total"][0]
    assert (total["n_samples"], total["n_individuals"], total["n_cxsv_calls"]) == ("3", "2", "4")


JCOLS = ["chrom1", "pos1", "chrom2", "pos2", "svtype", "svlen", "locus_id", "supp_vec", "dosage",
         "n_records", "record_ids"]


def junction_table(path, rows):
    """rows: (pos1, pos2, svtype, locus_id) on chr1."""
    with open(path, "w") as fh:
        fh.write("\t".join(JCOLS) + "\n")
        for k, (a, b, svtype, locus) in enumerate(rows):
            fh.write(f"chr1\t{a}\tchr1\t{b}\t{svtype}\t{b - a}\t{locus}\t111\t1\t1\tr{k}\n")
    return path


def panel_release(tmp_path, samples):
    """samples: (sample, individual, platform, loci [(chr, s, e)], junction rows)."""
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("chr1\t1000000\t6\t60\t61\nchr2\t1000000\t1016674\t60\t61\n")
    man = tmp_path / "manifest.tsv"
    with open(man, "w") as fh:
        fh.write("cohort\tsample_id\tindividual_id\tplatform\tcoverage_x\tmaster_table\tjunctions\n")
        for sid, ind, plat, loci, jrows in samples:
            mt = master_table(tmp_path / f"{sid}_mt.tsv", loci)
            jt = junction_table(tmp_path / f"{sid}_j.tsv", jrows)
            fh.write(f"c\t{sid}\t{ind}\t{plat}\t30\t{mt}\t{jt}\n")
    reg = tmp_path / "individuals.tsv"
    reg.write_text("individual_id\tsex\tsuperpopulation\n" + "".join(
        f"{ind}\tfemale\t{'EUR' if k % 2 else 'AFR'}\n" for k, (_, ind, *_r) in enumerate(samples)))
    out = {k: tmp_path / f"{k}.tsv" for k in ("calls", "loci", "sharing", "summary", "sites", "carriers", "inds")}
    subprocess.run([sys.executable, SCRIPTS / "build_panel_loci.py", "--manifest", man, "--calls", out["calls"],
                    "--loci", out["loci"], "--sharing", out["sharing"], "--summary", out["summary"]], check=True)
    subprocess.run([sys.executable, SCRIPTS / "build_panel_sites.py", "--manifest", man, "--individuals", reg,
                    "--fai", fai, "--tolerance", "500", "--version", "test1", "--sites", out["sites"],
                    "--carriers", out["carriers"], "--panel-individuals", out["inds"]], check=True)
    return out


def test_sites_follow_the_shared_schema(tmp_path):
    out = panel_release(tmp_path, [
        ("A", "A", "ONT", [], [(1000, 9000, "DEL", ".")]),
        ("B", "B", "PACBIO", [], [(1100, 9050, "DEL", "."), (20000, 20000, "INS", ".")]),
    ])
    rows = read(out["sites"])
    assert list(rows[0].keys())[:11] == ["chrom1", "pos1", "chrom2", "pos2", "svtype", "svlen", "ac", "an", "af",
                                         "source", "version"]
    assert {"af_ONT", "af_PACBIO", "af_AFR", "af_EUR"} <= set(rows[0])
    d = [r for r in rows if r["svtype"] == "DEL"][0]
    assert (d["ac"], d["an"], d["af"], d["source"], d["version"]) == ("2", "4", "0.500000", "cxsv_panel", "test1")
    assert (d["af_ONT"], d["af_PACBIO"]) == ("0.500000", "0.500000")
    car = [r for r in read(out["carriers"]) if r["svtype"] == "DEL"][0]
    assert (car["n_carriers"], car["carriers"]) == ("2", "A,B")


def classify(tmp_path, out, loci, jrows, **opt):
    mt = master_table(tmp_path / "patient_mt.tsv", loci)
    jt = junction_table(tmp_path / "patient_j.tsv", jrows)
    res = {k: tmp_path / f"p_{k}" for k in ("all", "private", "bed", "junctions")}
    args = [sys.executable, SCRIPTS / "classify_patient_cxsv.py", "--master-table", mt, "--junctions", jt,
            "--panel-carriers", out["carriers"], "--panel-sites", out["sites"],
            "--panel-individuals", out["inds"], "--panel-calls", out["calls"],
            "--out-all", res["all"], "--out-private", res["private"], "--out-bed", res["bed"],
            "--out-junctions", res["junctions"]]
    for k, v in opt.items():
        args += [f"--{k.replace('_', '-')}", str(v)]
    subprocess.run(args, check=True)
    return read(res["all"]), res


def test_patient_labels_by_junctions_with_hotspot_flag(tmp_path):
    common = [(10_000, 20_000, "INV"), (12_000, 25_000, "DEL"), (21_000, 30_000, "INV")]
    decoy = [(41_000, 52_000, "INV"), (46_000, 58_000, "DUP")]            # another event at the 40-60 kb hotspot
    out = panel_release(tmp_path, [
        ("A", "A", "ONT", [("chr1", 10_000, 30_000)], [(*j, "chr1:10000-30000") for j in common]),
        ("B", "B", "ONT", [("chr1", 10_000, 30_000), ("chr1", 41_000, 58_000)],
         [(*j, "chr1:10000-30000") for j in common] + [(*j, "chr1:41000-58000") for j in decoy]),
        ("C", "C", "PACBIO", [], [(70_000, 72_000, "DEL", ".")]),        # common simple SV, not a CxSV
    ])
    patient_loci = [("chr1", 10_000, 30_000), ("chr1", 40_000, 60_000), ("chr1", 69_000, 80_000)]
    jrows = [(*j, "chr1:10000-30000") for j in common]
    jrows += [(40_000, 50_000, "INV", "chr1:40000-60000"), (44_000, 60_000, "INV", "chr1:40000-60000")]
    jrows += [(69_900, 72_100, "DEL", "chr1:69000-80000"), (75_000, 80_000, "INV", "chr1:69000-80000")]
    jrows += [(90_000, 90_000, "INS", ".")]
    rows, res = classify(tmp_path, out, patient_loci, jrows)
    got = [(r["Panel_Status"], r["Panel_Carriers"], r["Panel_Junctions"], r["Panel_Junctions_Seen"],
            r["Panel_Region_Carriers"], r["Panel_Hotspot"]) for r in rows]
    assert got == [("public", "2", "3", "3", "2", "no"),
                   ("private", "0", "2", "0", "1", "yes"),               # decoy masked it under region overlap
                   ("shared_junctions", "0", "2", "1", "0", "no")]
    assert rows[0]["Panel_Carrier_Freq"] == f"{2 / 3:.6f}"
    assert rows[0]["Panel_Carriers_By_Platform"] == "ONT:2"
    assert res["bed"].read_text() == "chr1\t40000\t60000\nchr1\t69000\t80000\n"
    js = read(res["junctions"])
    assert len(js) == 8 and js[-1]["panel_status"] == "private" and js[0]["panel_carriers"] == "2"
    rows2, _ = classify(tmp_path, out, patient_loci, jrows, min_carriers=3)
    assert rows2[0]["Panel_Status"] == "shared_junctions"


def test_site_conventions_shared_with_catalogues(tmp_path):
    """BND svlen '.', non-BND sites < 50 bp dropped (sv_params.min_sv_len), unsigned svlen."""
    out = panel_release(tmp_path, [
        ("A", "A", "ONT", [], [(1000, 1030, "DEL", "."), (5000, 9000, "BND", "."), (20000, 26000, "INV", ".")]),
    ])
    rows = read(out["sites"])
    assert [(r["svtype"], r["svlen"]) for r in rows] == [("BND", "."), ("INV", "6000")]
    assert [r["svlen"] for r in read(out["carriers"])] == [".", "6000"]
