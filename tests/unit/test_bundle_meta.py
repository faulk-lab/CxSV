"""bundle.json tag and release-manifest fields (workflow/lib/cxsv_bundle_meta.py, item 16n)."""
import hashlib
import json
import subprocess

import pytest

import cxsv_bundle_meta as bm
from conftest import REPO


def git(repo, *args):
    return subprocess.run(["git", "-C", str(repo), *args], check=True, capture_output=True, text=True).stdout.strip()


@pytest.fixture
def repo(tmp_path):
    r = tmp_path / "repo"
    r.mkdir()
    git(r, "init", "-q")
    git(r, "config", "user.email", "ci@example.org")
    git(r, "config", "user.name", "ci")
    (r / "a.txt").write_text("1\n")
    git(r, "add", "a.txt")
    git(r, "commit", "-q", "-m", "one")
    return r


def test_tags_at_commit_peels_annotated_tags(repo):
    head = git(repo, "rev-parse", "HEAD")
    assert bm.tags_at_commit(str(repo), head) == []
    git(repo, "tag", "-a", "panel-v9.9", "-m", "release")
    git(repo, "tag", "light")
    tags = bm.tags_at_commit(str(repo), head)
    assert [t["name"] for t in tags] == ["light", "panel-v9.9"]
    ann = [t for t in tags if t["name"] == "panel-v9.9"][0]
    assert ann["type"] == "tag" and ann["object"] == git(repo, "rev-parse", "panel-v9.9")
    assert [t for t in tags if t["name"] == "light"][0]["object"] == head


def test_source_tag_prefers_the_release_tag():
    tags = [{"name": "light", "object": "c1", "type": "commit"}, {"name": "panel-v1.2", "object": "t1", "type": "tag"}]
    f, w = bm.source_tag({"git_tag": "panel-v1.2", "git_tag_object": "t1"}, tags)
    assert (f["git_tag"], f["git_tag_object"], f["git_tag_type"]) == ("panel-v1.2", "t1", "tag")
    assert f["git_tags_at_commit"] == ["light", "panel-v1.2"] and w == []


def test_source_tag_untagged_and_ambiguous_and_moved():
    f, w = bm.source_tag({}, [])
    assert f["git_tag"] is None and w == []
    f, w = bm.source_tag({"git_tag": "panel-v1.2"}, [])
    assert f["git_tag"] is None and "bundle commit is not" in w[0]
    f, w = bm.source_tag({}, [{"name": "a", "object": "1", "type": "tag"}, {"name": "b", "object": "2", "type": "tag"}])
    assert f["git_tag"] is None and "several tags" in w[0]
    f, w = bm.source_tag({}, [{"name": "only", "object": "1", "type": "tag"}])
    assert f["git_tag"] == "only" and w == []
    f, w = bm.source_tag({"git_tag": "v", "git_tag_object": "old"}, [{"name": "v", "object": "new", "type": "tag"}])
    assert f["git_tag"] == "v" and "tag moved" in w[0]


def test_parse_tag_lines_matches_the_make_bundle_format():
    text = "panel-v1.2\t101023e4\ttag\nb\t84d380e0\tcommit\n"
    assert [t["name"] for t in bm.parse_tag_lines(text)] == ["b", "panel-v1.2"]


def test_release_fields(tmp_path):
    (tmp_path / "release_manifest.tsv").write_text("file\tbytes\tsha256\n")
    rel = {"release_id": "r", "source": {"commit": "c", "git_tag": "panel-v1.2", "git_tag_object": "t1"},
           "finalised_utc": "2026-11-09T17:04:02Z"}
    (tmp_path / "release.json").write_text(json.dumps(rel))
    f = bm.release_fields(str(tmp_path))
    assert f["git_tag"] == "panel-v1.2" and f["finalised_utc"].startswith("2026-11-09")
    assert f["release_manifest_sha256"] == hashlib.sha256(b"file\tbytes\tsha256\n").hexdigest()
    assert f["release_json_sha256"] == hashlib.sha256((tmp_path / "release.json").read_bytes()).hexdigest()
    (tmp_path / "release.json").write_text(json.dumps({"source": {"commit": "c"}}))
    f = bm.release_fields(str(tmp_path))
    assert f["git_tag"] is None and f["finalised_utc"] is None


def test_make_bundle_writes_the_fields():
    sh = (REPO / "scripts/make_bundle.sh").read_text()
    assert "for-each-ref --points-at \"$COMMIT\" --format='%(refname:short)%09%(objectname)%09%(objecttype)'" in sh
    assert 'import cxsv_bundle_meta as bm' in sh
    assert 'meta["source"].update(tag_fields)' in sh and "bm.release_fields(" in sh
