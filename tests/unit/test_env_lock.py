"""The three descriptions of the pipeline environment must agree (roadmap item 2).

envs/cxsv-pipeline.yaml (pins), envs/cxsv-pipeline.lock.yaml (solved versions)
and envs/cxsv-pipeline.explicit.txt (package files the container image is built
from). They drifted once: the lock kept the pre-`mamba install minimap2=2.24`
solve while the env had been re-solved, so a container built from the lock would
not have matched the env that built the panel.
"""
import re

import yaml

from conftest import REPO

ENVS = REPO / "envs"


def explicit_versions():
    out = {}
    for line in (ENVS / "cxsv-pipeline.explicit.txt").read_text().splitlines():
        if line.startswith("https://"):
            base = line.split("#")[0].rsplit("/", 1)[1]
            base = re.sub(r"\.(tar\.bz2|conda)$", "", base)
            name, version, _build = base.rsplit("-", 2)
            out[name] = version
    return out


def lock_versions():
    deps = yaml.safe_load((ENVS / "cxsv-pipeline.lock.yaml").read_text())["dependencies"]
    return dict(d.split("=", 1) for d in deps if isinstance(d, str))


def test_lock_and_explicit_list_the_same_packages_and_versions():
    exp, lock = explicit_versions(), lock_versions()
    assert exp, "explicit list is empty"
    assert exp == lock, {k: (exp.get(k), lock.get(k)) for k in set(exp) | set(lock) if exp.get(k) != lock.get(k)}


def test_yaml_pins_hold_in_the_lock():
    lock = lock_versions()
    deps = yaml.safe_load((ENVS / "cxsv-pipeline.yaml").read_text())["dependencies"]
    for d in deps:
        if not isinstance(d, str):          # the pip: section (test below)
            continue
        m = re.fullmatch(r"([\w.-]+)=([\w.]+)", d)
        if m:
            name, pin = m.groups()
            v = lock.get(name)
            assert v is not None and (v == pin or v.startswith(pin + ".")), f"{name}: pinned {pin}, lock {v}"


def test_image_is_built_from_the_explicit_list():
    df = (REPO / "container" / "Dockerfile").read_text()
    assert "envs/cxsv-pipeline.explicit.txt" in df and re.search(r"^FROM \S+@sha256:[0-9a-f]{64}", df, re.M)
    assert "envs/cxsv-pipeline.pip.txt" in df and "--require-hashes" in df and "--no-deps" in df


def test_dockerignore_admits_every_copy_source():
    """container/Dockerfile.dockerignore excludes everything (`*`) and re-admits the files the
    image needs. panel-v1.3: the new `COPY envs/cxsv-pipeline.pip.txt` was not re-admitted, so
    `make_bundle.sh image` failed ("not found") on main 1c9c80d. Every COPY source must be."""
    df = (REPO / "container" / "Dockerfile").read_text()
    ign = (REPO / "container" / "Dockerfile.dockerignore").read_text().splitlines()
    rules = [l.strip() for l in ign if l.strip() and not l.lstrip().startswith("#")]
    assert rules[0] == "*"
    admitted = {r[1:] for r in rules if r.startswith("!")}
    sources = []
    for m in re.finditer(r"^COPY\s+(?!--from)(.+)$", df, re.M):
        parts = [p for p in m.group(1).split() if not p.startswith("--")]
        sources += parts[:-1]
    assert sources, "no COPY lines found"
    assert set(sources) <= admitted, set(sources) - admitted


def pip_lock():
    """{name: (version, sha256)} of envs/cxsv-pipeline.pip.txt."""
    text = (ENVS / "cxsv-pipeline.pip.txt").read_text().replace("\\\n", " ")
    out = {}
    for line in text.splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            m = re.fullmatch(r"([\w.-]+)==([\w.]+)\s+--hash=sha256:([0-9a-f]{64})", line)
            assert m, line
            out[m.group(1).lower()] = (m.group(2), m.group(3))
    return out


def test_pip_section_matches_the_hashed_pip_lock():
    """The clinician-report wheels (roadmap item 16l): yaml pip: pins == pip.txt pins, every
    pip package hash-pinned, none of them also a conda package of the lock."""
    deps = yaml.safe_load((ENVS / "cxsv-pipeline.yaml").read_text())["dependencies"]
    pip_sec = [d["pip"] for d in deps if isinstance(d, dict) and "pip" in d]
    assert len(pip_sec) == 1 and "pip" in deps
    pins = dict(p.split("==", 1) for p in pip_sec[0])
    lock = pip_lock()
    assert {k.lower(): v for k, v in pins.items()} == {k: v for k, (v, _) in lock.items()}
    assert {"python-docx", "reportlab", "lxml"} <= set(lock)
    conda = lock_versions()
    assert not set(lock) & {k.lower() for k in conda}, set(lock) & set(conda)
    assert "pip" in conda
