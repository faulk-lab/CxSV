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
        m = re.fullmatch(r"([\w.-]+)=([\w.]+)", d)
        if m:
            name, pin = m.groups()
            v = lock.get(name)
            assert v is not None and (v == pin or v.startswith(pin + ".")), f"{name}: pinned {pin}, lock {v}"


def test_image_is_built_from_the_explicit_list():
    df = (REPO / "container" / "Dockerfile").read_text()
    assert "envs/cxsv-pipeline.explicit.txt" in df and re.search(r"^FROM \S+@sha256:[0-9a-f]{64}", df, re.M)
