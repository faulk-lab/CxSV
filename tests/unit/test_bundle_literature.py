"""Literature set selection for the offline bundle (scripts/lib/literature_select.sh).

Until panel-v1.2, make_bundle.sh took `ls -d <resources>/literature/lit-v* | sort -V | tail -1`.
On the lab NAS that matched the build log lit-v1.0.build.log (a file), which sorts after
the lit-v1.0 directory, so every bundle since item 11 was built without the literature set
and without a message. These tests rebuild that layout and pin the fixed behaviour:
only directories named lit-v* with a MANIFEST.tsv are candidates, a missing set stops the
build unless --no-literature is given, and the validation checks that the set was shipped.
Needs only bash and coreutils.
"""
import hashlib
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
LIB = REPO / "scripts/lib/literature_select.sh"
MAKE_BUNDLE = REPO / "scripts/make_bundle.sh"
VALIDATE = REPO / "tests/lab/validate_bundle.sh"


def _select(resources, literature="", no_literature="0"):
    return subprocess.run(["bash", "-c", f'source "{LIB}"; cxsv_select_literature "$1" "$2" "$3"', "_",
                           str(resources), str(literature), no_literature],
                          capture_output=True, text=True)


def _verify(d):
    return subprocess.run(["bash", "-c", f'source "{LIB}"; cxsv_verify_literature "$1"', "_", str(d)],
                          capture_output=True, text=True)


def _make_set(d: Path, files=None):
    d.mkdir(parents=True)
    files = files or {"literature_events.tsv": "event_id\nE1\n", "literature_junctions.tsv": "j\n1\n"}
    rows = ["file\tbytes\tsha256"]
    for name, text in files.items():
        (d / name).write_text(text)
        rows.append(f"{name}\t{len(text.encode())}\t{hashlib.sha256(text.encode()).hexdigest()}")
    (d / "MANIFEST.tsv").write_text("\n".join(rows) + "\n")
    return d


def _lab_layout(tmp_path: Path) -> Path:
    """The NAS layout of 2026-11: the set, its build log, a superseded build, downloads."""
    lit = tmp_path / "resources/literature"
    _make_set(lit / "lit-v1.0")
    (lit / "lit-v1.0.build.log").write_text("built\n")
    _make_set(lit / "_superseded_lit-v1.0_build-8949cef")
    (lit / "_downloads").mkdir()
    return tmp_path / "resources"


def test_scripts_parse():
    for p in (LIB, MAKE_BUNDLE, VALIDATE):
        r = subprocess.run(["bash", "-n", str(p)], capture_output=True, text=True)
        assert r.returncode == 0, f"{p}: {r.stderr}"


def test_old_selection_picked_the_build_log(tmp_path):
    """Documents the bug: the pre-fix expression returns the log file, not the set."""
    res = _lab_layout(tmp_path)
    r = subprocess.run(["bash", "-c", f'ls -d "{res}"/literature/lit-v* | sort -V | tail -1'],
                       capture_output=True, text=True)
    assert Path(r.stdout.strip()).name == "lit-v1.0.build.log"


def test_directory_is_chosen_over_the_build_log(tmp_path):
    res = _lab_layout(tmp_path)
    r = _select(res)
    assert r.returncode == 0, r.stderr
    assert Path(r.stdout.strip()) == res / "literature/lit-v1.0"


def test_newest_version_wins_and_needs_a_manifest(tmp_path):
    res = _lab_layout(tmp_path)
    lit = res / "literature"
    _make_set(lit / "lit-v1.2")
    _make_set(lit / "lit-v1.10")
    (lit / "lit-v2.0").mkdir()                      # no MANIFEST.tsv: not a candidate
    (lit / "lit-v3.0.build.log").write_text("x\n")  # a file: not a candidate
    r = _select(res)
    assert r.returncode == 0, r.stderr
    assert Path(r.stdout.strip()).name == "lit-v1.10"


def test_no_set_stops_the_build(tmp_path):
    lit = tmp_path / "resources/literature"
    lit.mkdir(parents=True)
    (lit / "lit-v1.0.build.log").write_text("built\n")
    r = _select(tmp_path / "resources")
    assert r.returncode != 0
    assert "--no-literature" in r.stderr and r.stdout == ""
    r = _select(tmp_path / "missing")
    assert r.returncode != 0


def test_no_literature_flag(tmp_path):
    r = _select(tmp_path / "missing", "", "1")
    assert r.returncode == 0 and r.stdout == ""
    r = _select(tmp_path / "missing", str(tmp_path), "1")
    assert r.returncode != 0 and "mutually exclusive" in r.stderr


def test_explicit_literature_must_have_a_manifest(tmp_path):
    d = tmp_path / "lit-v9"
    d.mkdir()
    assert _select(tmp_path, d).returncode != 0
    assert _select(tmp_path, tmp_path / "absent").returncode != 0
    _make_set(tmp_path / "ok/lit-v1.0")
    r = _select(tmp_path, tmp_path / "ok/lit-v1.0")
    assert r.returncode == 0 and Path(r.stdout.strip()) == tmp_path / "ok/lit-v1.0"


def test_verify_catches_a_changed_file(tmp_path):
    d = _make_set(tmp_path / "lit-v1.0")
    assert _verify(d).returncode == 0
    (d / "literature_events.tsv").write_text("event_id\nE2\n")
    assert _verify(d).returncode != 0
    (d / "MANIFEST.tsv").write_text("file\tbytes\tsha256\n")
    assert _verify(d).returncode != 0


def test_make_bundle_uses_the_helper_and_stops_early():
    s = MAKE_BUNDLE.read_text()
    assert 'source "$REPO/scripts/lib/literature_select.sh"' in s
    assert "--no-literature) NO_LITERATURE=1" in s
    assert "ls -d \"$RESOURCES\"/literature/lit-v*" not in s
    # chosen before the commit is resolved, i.e. before any copying
    assert s.index("cxsv_select_literature") < s.index("git -C \"$REPO\" archive")
    assert "cxsv_verify_literature \"$LITERATURE\"" in s


def test_validation_has_a_literature_check():
    s = VALIDATE.read_text()
    assert 'hdr "9. literature set shipped and verified"' in s
    assert "cxsv_verify_literature \"$LIT_DIR\"" in s
    assert "cxsv_literature_hits.tsv" in s
    # check 6 runs with the same literature set as the bundle
    assert "literature_set:\\n  root:" in s


def test_bundle_literature_code_is_not_fingerprinted():
    import sys
    sys.path.insert(0, str(REPO / "workflow/lib"))
    import cxsv_fingerprint
    for f in cxsv_fingerprint.CALLING_CODE:
        assert "literature" not in f and "make_bundle" not in f
