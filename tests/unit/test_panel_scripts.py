"""Panel aggregation and patient labelling scripts."""
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


def classify(tmp_path, out, patient_loci, window=1000, min_carriers=1):
    t = master_table(tmp_path / "patient.tsv", patient_loci)
    res = {k: tmp_path / f"p_{k}" for k in ("all", "private", "bed")}
    subprocess.run([sys.executable, SCRIPTS / "classify_patient_cxsv.py", "--master-table", t,
                    "--panel-calls", out["calls"], "--panel-summary", out["summary"],
                    "--window", str(window), "--min-carriers", str(min_carriers),
                    "--out-all", res["all"], "--out-private", res["private"], "--out-bed", res["bed"]],
                   check=True)
    return read(res["all"]), res


def test_patient_labels_window_and_min_carriers(tmp_path):
    out = build_panel(tmp_path, [
        ("c", "A", "A", "ONT", [("chr1", 10_000, 20_000)]),
        ("c", "B", "B", "ONT", [("chr1", 10_500, 19_000)]),
        ("c", "C", "C", "PACBIO", [("chr3", 100, 900)]),
    ])
    rows, res = classify(tmp_path, out, [
        ("chr1", 12_000, 15_000),      # inside both A and B
        ("chr1", 20_800, 22_000),      # 800 bp past A's end: within the 1 kb window
        ("chr1", 25_000, 26_000),      # nothing nearby
        ("chr3", 950, 960),            # 50 bp past C
    ])
    got = [(r["Panel_Carriers"], r["Panel_Status"]) for r in rows]
    assert got == [("2", "public"), ("1", "public"), ("0", "private"), ("1", "public")]
    assert rows[0]["Panel_Carrier_Freq"] == f"{2 / 3:.6f}"
    assert rows[0]["Panel_Carriers_By_Platform"] == "ONT:2"
    assert res["bed"].read_text() == "chr1\t25000\t26000\n"
    rows2, _ = classify(tmp_path, out, [("chr1", 12_000, 15_000), ("chr1", 20_800, 22_000)], min_carriers=2)
    assert [r["Panel_Status"] for r in rows2] == ["public", "private"]
