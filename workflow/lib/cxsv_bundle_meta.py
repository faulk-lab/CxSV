"""Provenance fields of an offline bundle's bundle.json (scripts/make_bundle.sh).

Until panel-v1.2 bundle.json recorded the commit but not its git tag, and named
the panel release without the sha256 of its release_manifest.tsv; both were
only inside the embedded panel_release/<id>/release.json (roadmap item 16n).
These helpers add, without changing any other field:

    source.git_tag, source.git_tag_object, source.git_tag_type, source.git_tags_at_commit
        the tag that points at the bundle commit (the release's tag when several do;
        null when the commit is untagged, e.g. an untagged test bundle)
    panel_release.git_tag, panel_release.git_tag_object, panel_release.finalised_utc
        the release's own tag fields (release.json source.*, set by finalising)
    panel_release.release_manifest_sha256, panel_release.release_json_sha256
        sha256 of the shipped release_manifest.tsv / release.json

Plain Python, unit-tested in tests/unit/test_bundle_meta.py.
"""
from __future__ import annotations

import hashlib
import json
import os
import subprocess


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def tags_at_commit(repo: str, commit: str) -> list[dict]:
    """Tags pointing at `commit` (annotated tags peeled): name, object, type."""
    out = subprocess.run(["git", "-C", repo, "for-each-ref", "--points-at", commit,
                          "--format=%(refname:short)\t%(objectname)\t%(objecttype)", "refs/tags"],
                         check=True, capture_output=True, text=True).stdout
    return parse_tag_lines(out)


def parse_tag_lines(text: str) -> list[dict]:
    """The tab-separated output of tags_at_commit's git command, already written to a file."""
    tags = []
    for line in text.splitlines():
        if line.strip():
            name, obj, typ = line.split("\t")
            tags.append({"name": name, "object": obj, "type": typ})
    return sorted(tags, key=lambda t: t["name"])


def source_tag(release_source: dict, tags: list[dict]) -> tuple[dict, list[str]]:
    """source.* tag fields for bundle.json plus warnings for the build log."""
    warnings = []
    by_name = {t["name"]: t for t in tags}
    rel_tag = release_source.get("git_tag")
    chosen = None
    if rel_tag and rel_tag in by_name:
        chosen = by_name[rel_tag]
    elif len(tags) == 1:
        chosen = tags[0]
    elif len(tags) > 1:
        warnings.append(f"bundle commit has several tags {sorted(by_name)} and none is the release's "
                        f"({rel_tag}); source.git_tag left null")
    if rel_tag and rel_tag not in by_name:
        warnings.append(f"the release is tagged {rel_tag} but the bundle commit is not "
                        f"(tags at the commit: {sorted(by_name) or 'none'})")
    if chosen and rel_tag and rel_tag == chosen["name"] and release_source.get("git_tag_object") \
            and release_source["git_tag_object"] != chosen["object"]:
        warnings.append(f"tag {rel_tag} now points to object {chosen['object']}, the release recorded "
                        f"{release_source['git_tag_object']} (tag moved?)")
    fields = {"git_tag": chosen["name"] if chosen else None,
              "git_tag_object": chosen["object"] if chosen else None,
              "git_tag_type": chosen["type"] if chosen else None,
              "git_tags_at_commit": [t["name"] for t in tags]}
    return fields, warnings


def release_fields(release_dir: str) -> dict:
    """panel_release.* tag and checksum fields from a (shipped) release directory."""
    rel = json.load(open(os.path.join(release_dir, "release.json")))
    src = rel.get("source", {})
    return {"git_tag": src.get("git_tag"), "git_tag_object": src.get("git_tag_object"),
            "finalised_utc": rel.get("finalised_utc"),
            "release_manifest_sha256": sha256_file(os.path.join(release_dir, "release_manifest.tsv")),
            "release_json_sha256": sha256_file(os.path.join(release_dir, "release.json"))}
