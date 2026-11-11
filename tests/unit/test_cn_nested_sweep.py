"""nested_min_z is sweepable without code edits (panel-v1.3 item j, decision D5).

PC L01's nested triplication reads z 5.73 against the female chrX reference, below
nested_min_z 6.0; the v1.3 design check sweeps 5.0-6.0. The default stays 6.0. These tests
pin the two override routes: cn_segment.py --set (direct runs, the design check) and a
Snakemake --configfile overlay (pipeline runs; Snakemake merges config recursively), and
that the threshold governs the nested call. Synthetic depth only.
"""
import gzip
import json
import subprocess
import sys
import textwrap

import numpy as np
import pytest
import yaml

from conftest import REPO
import cxsv_cn as cn  # noqa: E402

NB = 3000
PAR = {"chrX": [(0, 200_000)]}
S = REPO / "workflow" / "scripts"


def _fine(rng, cov, events=()):
    out = {}
    for c in ("chr1", "chr2", "chrX"):
        lam = np.full(NB, float(cov))                 # female: two copies everywhere
        for ec, s, e, k in events:
            if ec == c:
                lam[s:e] = cov * k / 2.0
        out[c] = rng.poisson(lam).astype(np.float32)
    return out


def _write(path, fine):
    with gzip.open(path, "wt") as fh:
        for c, a in fine.items():
            fh.writelines(f"{c}\t{i * 1000}\t{i * 1000 + 1000}\t{v:.2f}\n" for i, v in enumerate(a))


# a DUP (3 copies, 800 kb) carrying a compressed TRP (+0.6 copy, 90 kb), as PC L01
EVENT = [("chr1", 1000, 1400, 3), ("chr1", 1400, 1490, 3.6), ("chr1", 1490, 1800, 3)]


def _nested_calls(segs):
    return [s for s in segs if s.get("nested")]


def test_default_is_six_and_config_matches():
    assert cn.DEFAULTS["nested_min_z"] == 6.0
    cfg = yaml.safe_load(open(REPO / "config" / "calling.yaml"))
    assert cfg["cn_params"]["segmentation"]["nested_min_z"] == 6.0


def test_nested_call_follows_the_threshold():
    rng = np.random.default_rng(21)
    members = []
    for _ in range(8):
        coarse = cn.to_coarse(_fine(rng, 30))
        gmed = cn.genome_median(coarse)
        members.append((cn.normalised(coarse, gmed, "female", PAR), "female"))
    pon = cn.Pon.build(members, PAR)
    fine = _fine(rng, 30, EVENT)
    n = {}
    for z in (0.5, 5.0, 5.5, 6.0, 1e6):
        segs, _ = cn.segment_genome(fine, "female", PAR, pon, {**cn.DEFAULTS, "nested_min_z": z})
        n[z] = len(_nested_calls(segs))
    vals = [n[z] for z in sorted(n)]
    assert vals == sorted(vals, reverse=True)          # a lower threshold never calls fewer nested levels
    assert n[0.5] >= 1 and n[1e6] == 0                 # the threshold decides the nested call


def _pon_dir(tmp_path, rng):
    rows = ["key\tindividual_id\tstratum\tdepth"]
    for i in range(6):
        p = tmp_path / f"n{i}.bed.gz"
        _write(p, _fine(rng, 30))
        rows.append(f"c/n{i}\tI{i}\tT\t{p}")
    (tmp_path / "members.tsv").write_text("\n".join(rows) + "\n")
    params = {"cn_params": {"segmentation": dict(cn.DEFAULTS), "par": {"chrX": [[0, 200000]]}, "linkage": {}}}
    (tmp_path / "p.json").write_text(json.dumps(params))
    subprocess.run([sys.executable, S / "cn_pon.py", "--members", tmp_path / "members.tsv", "--params", tmp_path / "p.json",
                    "--out-dir", tmp_path / "pon", "--threads", "1"], check=True, capture_output=True)


def test_cn_segment_set_override(tmp_path):
    rng = np.random.default_rng(22)
    _pon_dir(tmp_path, rng)
    case = tmp_path / "case.bed.gz"
    _write(case, _fine(rng, 30, EVENT))
    base = [sys.executable, S / "cn_segment.py", "--depth", case, "--pon", tmp_path / "pon" / "T.npz",
            "--params", tmp_path / "p.json", "--stratum", "T"]
    out = {}
    for z in ("6.0", "5.5", "5.0", "0.5", "1000000"):
        r = subprocess.run(base + ["--set", f"segmentation.nested_min_z={z}", "--out", tmp_path / f"s{z}.tsv",
                                   "--stats", tmp_path / f"s{z}.json"], capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
        st = json.loads((tmp_path / f"s{z}.json").read_text())
        assert st["nested_min_z"] == float(z) and st["overrides"] == {"segmentation.nested_min_z": float(z)}
        out[float(z)] = sum(l.split("\t")[9] == "True" for l in (tmp_path / f"s{z}.tsv").read_text().splitlines()[1:])
    assert out[0.5] >= out[5.0] >= out[5.5] >= out[6.0] >= out[1e6] == 0 and out[0.5] >= 1
    # no --set: the params file decides and stats.json carries no overrides
    r = subprocess.run(base + ["--out", tmp_path / "d.tsv", "--stats", tmp_path / "d.json"], capture_output=True, text=True)
    st = json.loads((tmp_path / "d.json").read_text())
    assert r.returncode == 0 and st["nested_min_z"] == 6.0 and "overrides" not in st
    for bad in ("segmentation.nested_min_zz=5", "nosuch.nested_min_z=5", "nested_min_z5"):
        r = subprocess.run(base + ["--set", bad, "--out", tmp_path / "x.tsv", "--stats", tmp_path / "x.json"],
                           capture_output=True, text=True)
        assert r.returncode != 0 and "--set" in r.stderr


def test_snakemake_configfile_overlay_overrides_only_nested_min_z(tmp_path):
    pytest.importorskip("snakemake")
    out = tmp_path / "cn_params.json"
    (tmp_path / "Snakefile").write_text(textwrap.dedent(f"""
        import json
        configfile: "{REPO / 'config' / 'calling.yaml'}"
        json.dump(config["cn_params"], open("{out}", "w"))
        rule all:
            input: []
    """))
    (tmp_path / "nested_5.5.yaml").write_text("cn_params:\n  segmentation:\n    nested_min_z: 5.5\n")
    r = subprocess.run([sys.executable, "-m", "snakemake", "-s", tmp_path / "Snakefile", "-d", tmp_path, "-n",
                        "--cores", "1", "--configfile", tmp_path / "nested_5.5.yaml"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-2000:]
    got = json.loads(out.read_text())
    ref = yaml.safe_load(open(REPO / "config" / "calling.yaml"))["cn_params"]
    assert got["segmentation"]["nested_min_z"] == 5.5
    ref["segmentation"]["nested_min_z"] = 5.5
    assert got == ref                                   # every other cn_params value unchanged
