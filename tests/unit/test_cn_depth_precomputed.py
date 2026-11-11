"""cn_depth.sh precomputed mode (panel-v1.3 item a; docs/PLAN_panel_v1.3.md).

The full-panel depth pass (cn_depth_v1) is shared by every rebuild. In v1.2 the rule
hard-linked its files, so the build's outputs shared the source inodes: a rebuild that
touched its outputs changed 1,619 cn_depth_v1 mtimes and once stopped the run. The rule
now copies them and checks the copy against the recorded sha256. These tests pin that the
source is only read (mtime, inode and link count unchanged) and that a mismatching copy
or source fails the step.
"""
import gzip
import hashlib
import os
import shutil
import subprocess

import pytest

from conftest import REPO

SCRIPT = REPO / "workflow" / "scripts" / "cn_depth.sh"
FILES = ("depth.regions.bed.gz", "depth.regions.bed.gz.csi", "depth.mosdepth.summary.txt")
NEED = ("bash", "zcat", "sha256sum", "awk", "cp")

pytestmark = pytest.mark.skipif(any(shutil.which(t) is None for t in NEED),
                                reason="needs bash, zcat, sha256sum, awk, cp")


def _precomputed(tmp_path, bins=5, binsize=1000):
    ref = tmp_path / "ref.fa"
    ref.write_text(">chrT\n" + "A" * (bins * binsize) + "\n")
    (tmp_path / "ref.fa.fai").write_text(f"chrT\t{bins * binsize}\t6\t{bins * binsize}\t{bins * binsize + 1}\n")
    pre = tmp_path / "cn_depth_v1" / "cohort" / "S1"
    pre.mkdir(parents=True)
    body = "".join(f"chrT\t{i * binsize}\t{(i + 1) * binsize}\t30.0\n" for i in range(bins)).encode()
    with gzip.open(pre / FILES[0], "wb") as fh:
        fh.write(body)
    (pre / FILES[1]).write_bytes(b"csi-index-placeholder")
    (pre / FILES[2]).write_text("chrom\tlength\tbases\tmean\tmin\tmax\ntotal\t5000\t150000\t30.0\t0\t30\n")
    sha = hashlib.sha256((pre / FILES[0]).read_bytes()).hexdigest()
    (pre / "depth.sha256").write_text(sha + "\n")
    (pre / "depth.ok").write_text("")
    old = 1_600_000_000_000_000_000                  # fixed source mtimes (ns), far in the past
    for f in FILES:
        os.utime(pre / f, ns=(old, old))
    return ref, pre, sha


def _run(ref, out, pre):
    return subprocess.run(["bash", str(SCRIPT), "unused.cram", str(ref), str(out), "1000", "20", "1", str(pre)],
                          capture_output=True, text=True)


def test_precomputed_depth_is_copied_and_source_untouched(tmp_path):
    ref, pre, sha = _precomputed(tmp_path)
    before = {f: os.stat(pre / f) for f in FILES}
    out = tmp_path / "build" / "S1" / "cn"
    r = _run(ref, out, pre)
    assert r.returncode == 0, r.stderr + r.stdout
    assert "copy verified" in r.stdout
    for f in FILES:
        src, dst = os.stat(pre / f), os.stat(out / f)
        assert src.st_mtime_ns == before[f].st_mtime_ns, f          # source mtime unchanged
        assert src.st_nlink == before[f].st_nlink == 1, f           # no hard link was added
        assert src.st_ino == before[f].st_ino, f
        assert dst.st_ino != src.st_ino, f                          # the output is its own file
        assert (out / f).read_bytes() == (pre / f).read_bytes(), f
    assert (out / "depth.sha256").read_text().strip() == sha
    # a re-run (outputs present) replaces the copies and still leaves the source alone
    os.utime(out / FILES[0], None)
    r = _run(ref, out, pre)
    assert r.returncode == 0, r.stderr
    for f in FILES:
        assert os.stat(pre / f).st_mtime_ns == before[f].st_mtime_ns
        assert os.stat(pre / f).st_nlink == 1


def test_precomputed_depth_with_wrong_sha_fails(tmp_path):
    ref, pre, sha = _precomputed(tmp_path)
    (pre / "depth.sha256").write_text("0" * 64 + "\n")
    r = _run(ref, tmp_path / "out", pre)
    assert r.returncode != 0 and "differs from depth.sha256" in r.stderr


def test_precomputed_depth_with_wrong_bin_count_fails(tmp_path):
    ref, pre, sha = _precomputed(tmp_path)
    (tmp_path / "ref.fa.fai").write_text("chrT\t9000\t6\t9000\t9001\n")
    r = _run(ref, tmp_path / "out", pre)
    assert r.returncode != 0 and "bin count" in r.stderr


def test_script_no_longer_hard_links():
    body = SCRIPT.read_text()
    code = "\n".join(l for l in body.splitlines() if not l.lstrip().startswith("#"))
    assert " ln " not in code and "ln \"$PRE" not in code
    assert "cp --reflink=auto" in code
