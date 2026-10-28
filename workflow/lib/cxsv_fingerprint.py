"""Calling fingerprint: what must be identical between a panel and a patient run.

A patient CxSV is "private" only relative to a panel built the same way. If the
two sides differ in caller versions, calling/discovery parameters, reference
sequence, or the code that turns calls into CxSV loci, some patient calls look
private purely because of the difference. The fingerprint records exactly those
inputs; the panel release carries it and the patient workflow refuses a panel
whose fingerprint differs from its own.

Components
----------
params     calling and discovery parameters (CALLING_PARAM_KEYS), canonical JSON
reference  every contig's name, length and MD5 (from the htslib MD5 cache manifest)
tools      versions of the tools whose output defines the calls
code       SHA-256 of the files that define per-sample calling and CxSV loci,
           with comment and blank lines dropped so a comment edit does not
           invalidate a panel

Deliberately NOT included: thread counts, file paths, annotation resources
(they label loci but do not decide which loci exist) and report formatting.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

CALLING_PARAM_KEYS = ("sv_params", "filter_params", "survivor_params", "cxsv_params")

# Relative to the workflow/ directory. Everything that decides which CxSV loci a
# sample has; changing any of these needs a new panel release.
CALLING_CODE = (
    "rules/calling.smk",
    "rules/per_sample.smk",
    "scripts/vcf_to_bed.py",
    "scripts/classify_cxsv.py",
    "scripts/extract_cxsv_vcf.py",
    "scripts/annotate_cxsv.py",
)

TOOLS = ("sniffles", "cuteSV", "svim", "SURVIVOR", "bcftools", "bedtools", "samtools")


def calling_params(config: dict) -> dict:
    missing = [k for k in CALLING_PARAM_KEYS if k not in config]
    if missing:
        raise KeyError(f"config lacks calling parameter section(s): {missing}")
    return {k: config[k] for k in CALLING_PARAM_KEYS}


def _normalised_code(path: Path) -> bytes:
    keep = []
    for line in path.read_text().splitlines():
        s = line.strip()
        if s and not s.startswith("#"):
            keep.append(line.rstrip())
    return ("\n".join(keep) + "\n").encode()


def code_digest(workflow_dir: str | Path) -> dict:
    wd = Path(workflow_dir)
    per_file = {rel: hashlib.sha256(_normalised_code(wd / rel)).hexdigest() for rel in CALLING_CODE}
    combined = hashlib.sha256("".join(f"{k}:{v}\n" for k, v in sorted(per_file.items())).encode())
    return {"sha256": combined.hexdigest(), "files": per_file}


def reference_contigs(md5_manifest: str | Path) -> list[dict]:
    rows = []
    with open(md5_manifest) as fh:
        next(fh)
        for line in fh:
            sn, ln, m5 = line.rstrip("\n").split("\t")
            rows.append({"name": sn, "length": int(ln), "md5": m5})
    return rows


def read_tool_versions(path: str | Path) -> dict:
    out = {}
    with open(path) as fh:
        for line in fh:
            if line.strip():
                tool, version = line.rstrip("\n").split("\t", 1)
                out[tool] = version
    return out


def build(config: dict, workflow_dir, md5_manifest, tool_versions_tsv) -> dict:
    fp = {
        "params": calling_params(config),
        "reference": reference_contigs(md5_manifest),
        "tools": read_tool_versions(tool_versions_tsv),
        "code": code_digest(workflow_dir),
    }
    fp["sha256"] = digest(fp)
    return fp


def digest(fp: dict) -> str:
    body = {k: fp[k] for k in ("params", "reference", "tools", "code")}
    return hashlib.sha256(json.dumps(body, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def compare(ours: dict, theirs: dict) -> list[str]:
    """Human-readable differences; empty list means compatible."""
    diffs = []
    if ours.get("params") != theirs.get("params"):
        for section in CALLING_PARAM_KEYS:
            a, b = ours["params"].get(section, {}), theirs.get("params", {}).get(section, {})
            for key in sorted(set(a) | set(b)):
                if a.get(key) != b.get(key):
                    diffs.append(f"params.{section}.{key}: this run={a.get(key)!r} panel={b.get(key)!r}")
    ra = {c["name"]: (c["length"], c["md5"]) for c in ours.get("reference", [])}
    rb = {c["name"]: (c["length"], c["md5"]) for c in theirs.get("reference", [])}
    for name in sorted(set(ra) | set(rb)):
        if ra.get(name) != rb.get(name):
            diffs.append(f"reference.{name}: this run={ra.get(name)} panel={rb.get(name)}")
    ta, tb = ours.get("tools", {}), theirs.get("tools", {})
    for tool in sorted(set(ta) | set(tb)):
        if ta.get(tool) != tb.get(tool):
            diffs.append(f"tools.{tool}: this run={ta.get(tool)!r} panel={tb.get(tool)!r}")
    fa = ours.get("code", {}).get("files", {})
    fb = theirs.get("code", {}).get("files", {})
    for rel in sorted(set(fa) | set(fb)):
        if fa.get(rel) != fb.get(rel):
            diffs.append(f"code.{rel}: differs from the code the panel was built with")
    return diffs
