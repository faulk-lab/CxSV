"""Panel releases: identity, provenance and integrity (roadmap item 2).

A panel release is a directory the patient workflow reads and nothing else:

    release_id   <name>-<YYYYMMDD>-<first 12 hex digits of the calling fingerprint>
                 e.g. cxsv-panel-20261029-3f2a9c1d0b4e. Same fingerprint = the
                 same calling; the date and name tell apart rebuilds with more
                 samples.
    release.json metadata: id, source commit, fingerprint, counts by platform and
                 superpopulation, callability definition and summary, tool
                 versions, container runtime (docs/RELEASE_AND_BUNDLE.md)
    release_manifest.tsv   file, bytes, sha256 of EVERY other file in the
                 directory (relative paths, sorted). The patient workflow
                 refuses a release that does not verify: a missing, changed or
                 unlisted file (`verify`).

The source commit is read with git when the code is a git checkout, otherwise
from workflow/COMMIT, which `git archive` fills in (export-subst in
.gitattributes) — the form the code takes inside the offline bundle.
"""
from __future__ import annotations

import hashlib
import os
import re
import subprocess

MANIFEST = "release_manifest.tsv"
METADATA = "release.json"
MANIFEST_COLUMNS = ("file", "bytes", "sha256")
_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


def release_id(name: str, date: str, fingerprint_sha256: str) -> str:
    if not _NAME_RE.match(name or ""):
        raise ValueError(f"panel release name {name!r}: use letters, digits, '.', '_' or '-'")
    if not re.fullmatch(r"\d{8}", str(date)):
        raise ValueError(f"panel release date {date!r}: use YYYYMMDD")
    if not re.fullmatch(r"[0-9a-f]{64}", fingerprint_sha256 or ""):
        raise ValueError("calling fingerprint sha256 must be 64 hex digits")
    return f"{name}-{date}-{fingerprint_sha256[:12]}"


def sha256_file(path, block=1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(block), b""):
            h.update(chunk)
    return h.hexdigest()


def source_info(project_dir) -> dict:
    """Commit the code was taken from: git checkout, or a git-archive export."""
    try:
        commit = subprocess.run(["git", "-C", str(project_dir), "rev-parse", "HEAD"], capture_output=True,
                                text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(project_dir), "status", "--porcelain", "--untracked-files=no"],
                               capture_output=True, text=True, check=True).stdout.strip()
        return {"commit": commit, "dirty": bool(dirty), "method": "git"}
    except (OSError, subprocess.CalledProcessError):
        pass
    path = os.path.join(str(project_dir), "workflow", "COMMIT")
    if os.path.exists(path):
        text = open(path).read().strip()
        if re.fullmatch(r"[0-9a-f]{40}", text):
            return {"commit": text, "dirty": False, "method": "git-archive (workflow/COMMIT)"}
    return {"commit": "unknown", "dirty": None, "method": "unknown"}


def list_files(release_dir) -> list[str]:
    """Every regular file below release_dir except the manifest, as sorted relative paths."""
    out = []
    for root, dirs, files in os.walk(release_dir):
        dirs.sort()
        for f in files:
            rel = os.path.relpath(os.path.join(root, f), release_dir)
            if rel != MANIFEST:
                out.append(rel.replace(os.sep, "/"))
    return sorted(out)


def write_manifest(release_dir) -> str:
    """Checksum every file of the release; returns the manifest path."""
    path = os.path.join(release_dir, MANIFEST)
    rows = []
    for rel in list_files(release_dir):          # listed before the temp file exists
        full = os.path.join(release_dir, rel)
        rows.append(f"{rel}\t{os.path.getsize(full)}\t{sha256_file(full)}\n")
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        fh.write("\t".join(MANIFEST_COLUMNS) + "\n")
        fh.writelines(rows)
    os.replace(tmp, path)
    return path


def read_manifest(release_dir) -> list[dict]:
    with open(os.path.join(release_dir, MANIFEST)) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if tuple(header) != MANIFEST_COLUMNS:
            raise ValueError(f"{MANIFEST}: unexpected header {header}")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]


def verify(release_dir) -> dict:
    """Check a release against its manifest.

    Returns {"ok", "problems", "n_files", "bytes", "manifest_sha256"}. Problems:
    no manifest, a listed file missing / of another size / another sha256, a
    file present but not listed, or no release.json among the listed files."""
    problems = []
    mpath = os.path.join(release_dir, MANIFEST)
    if not os.path.isfile(mpath):
        return {"ok": False, "problems": [f"no {MANIFEST} in {release_dir}"], "n_files": 0, "bytes": 0,
                "manifest_sha256": None}
    rows = read_manifest(release_dir)
    listed = {r["file"] for r in rows}
    total = 0
    for r in rows:
        full = os.path.join(release_dir, r["file"])
        if not os.path.isfile(full):
            problems.append(f"missing: {r['file']}")
            continue
        size = os.path.getsize(full)
        total += size
        if str(size) != r["bytes"]:
            problems.append(f"size differs: {r['file']} ({size} bytes, manifest {r['bytes']})")
        elif sha256_file(full) != r["sha256"]:
            problems.append(f"sha256 differs: {r['file']}")
    for rel in list_files(release_dir):
        if rel not in listed:
            problems.append(f"not in manifest: {rel}")
    if METADATA not in listed:
        problems.append(f"{METADATA} is not listed in the manifest")
    return {"ok": not problems, "problems": problems, "n_files": len(rows), "bytes": total,
            "manifest_sha256": sha256_file(mpath)}
