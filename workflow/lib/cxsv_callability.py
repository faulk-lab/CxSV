"""Panel callability: could a panel sample have shown a junction at all?

Without callability, "absent from the panel" cannot be told apart from "the
panel could not have seen it there" (no reads, or only multi-mapping reads).

Per panel sample, the panel builder runs mosdepth once (rules/panel_callability.smk):

    mosdepth -n -x --mapq <sv_params.min_mapq> --quantize 0:<min_depth>: -f <ref> ...

and keeps the CALLABLE intervals (depth >= min_depth counting reads of MAPQ >=
min_mapq, the reads the SV callers use) as a BED: `callable.bed.gz`. These masks
are shipped in the panel release (callability/), so callability can be asked at
any coordinate, including at a patient's private junctions, which by
definition have no panel site.

Definitions (half-open, 0-based intervals; VCF-style 1-based positions)
------------------------------------------------------------------------
breakpoint boundary  b = pos (the junction lies between bases pos and pos+1,
                     1-based: after the last reference base kept, e.g. a DEL's
                     POS, and after a DEL's END)
breakpoint callable  a window of `flank` bp on at least ONE side of the
                     boundary, [b - flank, b) or [b, b + flank), lies entirely
                     inside one CALLABLE interval. One side is enough: reads
                     anchoring a junction come from either flank, and a
                     homozygous deletion carrier has no reads inside the
                     deleted segment.
junction callable    both breakpoints callable in the same sample
individual callable  at least one of its panel samples is callable
locus callable       callable at every junction of the locus

A carrier normally is callable; it can fall below min_depth at a breakpoint
(mosaic coverage dips), so carrier counts and callable counts are reported
side by side, never one capped by the other.
"""
from __future__ import annotations

import csv
import gzip
from collections import defaultdict

import numpy as np

MASK_SUFFIX = ".callable.bed.gz"
SAMPLES_TSV = "samples.tsv"
SAMPLES_COLUMNS = ("mask", "cohort", "sample_id", "individual_id", "platform", "callable_bp")


def mask_name(cohort: str, sample: str) -> str:
    """File name of one panel sample's mask inside release/callability/."""
    return f"{cohort}__{sample}{MASK_SUFFIX}"


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


class Mask:
    """CALLABLE intervals of one sample: per contig sorted, non-overlapping arrays."""

    def __init__(self, intervals: dict[str, list[tuple[int, int]]]):
        self.starts, self.ends = {}, {}
        for chrom, iv in intervals.items():
            iv = sorted(iv)
            merged = []
            for s, e in iv:
                if merged and s <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], e)
                else:
                    merged.append([s, e])
            arr = np.array(merged, dtype=np.int64).reshape(-1, 2)
            self.starts[chrom], self.ends[chrom] = arr[:, 0], arr[:, 1]

    @classmethod
    def from_bed(cls, path):
        iv = defaultdict(list)
        with _open(path) as fh:
            for line in fh:
                if not line.strip() or line.startswith(("#", "track")):
                    continue
                f = line.split("\t")
                iv[f[0]].append((int(f[1]), int(f[2])))
        return cls(iv)

    def callable_bp(self) -> int:
        return int(sum(int((self.ends[c] - self.starts[c]).sum()) for c in self.starts))

    def contains(self, chrom: str, ws, we):
        """Vectorised: is each window [ws, we) entirely inside one interval?"""
        ws, we = np.asarray(ws, dtype=np.int64), np.asarray(we, dtype=np.int64)
        if chrom not in self.starts or len(self.starts[chrom]) == 0:
            return np.zeros(ws.shape, dtype=bool)
        st, en = self.starts[chrom], self.ends[chrom]
        i = np.searchsorted(st, ws, side="right") - 1
        ok = i >= 0
        out = np.zeros(ws.shape, dtype=bool)
        out[ok] = en[i[ok]] >= we[ok]
        return out & (ws >= 0)

    def breakpoints_callable(self, chrom: str, pos, flank: int):
        """Vectorised breakpoint rule (module docstring) for 1-based positions."""
        b = np.asarray(pos, dtype=np.int64)
        return self.contains(chrom, b - flank, b) | self.contains(chrom, b, b + flank)


def junctions_callable(mask: Mask, junctions, flank: int) -> np.ndarray:
    """Boolean per junction (dicts with chrom1 pos1 chrom2 pos2): both ends callable."""
    n = len(junctions)
    end1 = np.zeros(n, dtype=bool)
    end2 = np.zeros(n, dtype=bool)
    by1, by2 = defaultdict(list), defaultdict(list)
    for k, j in enumerate(junctions):
        by1[j["chrom1"]].append(k)
        by2[j["chrom2"]].append(k)
    for chrom, idx in by1.items():
        end1[idx] = mask.breakpoints_callable(chrom, [int(junctions[k]["pos1"]) for k in idx], flank)
    for chrom, idx in by2.items():
        end2[idx] = mask.breakpoints_callable(chrom, [int(junctions[k]["pos2"]) for k in idx], flank)
    return end1 & end2


def read_samples(callability_dir) -> list[dict]:
    with open(f"{callability_dir}/{SAMPLES_TSV}", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def callable_individuals(callability_dir, junctions, flank: int):
    """For each junction, the set of panel individuals callable at both ends.

    Loads one sample mask at a time (memory: one mask), so it scales to large
    panels: the per-sample work is a vectorised search over all junctions."""
    out = [set() for _ in junctions]
    if not junctions:
        return out
    for s in read_samples(callability_dir):
        mask = Mask.from_bed(f"{callability_dir}/{s['mask']}")
        hit = junctions_callable(mask, junctions, flank)
        for k in np.flatnonzero(hit):
            out[k].add(s["individual_id"])
    return out


def callable_counts(callability_dir, junctions, flank: int, individuals: dict):
    """Per junction: n callable samples, n callable individuals, and individuals by platform.

    individuals: {individual_id: {"platforms": "ONT,PACBIO", ...}} (panel_individuals.tsv)."""
    n_samples = np.zeros(len(junctions), dtype=np.int64)
    inds = [set() for _ in junctions]
    by_platform = [defaultdict(set) for _ in junctions]
    for s in read_samples(callability_dir):
        mask = Mask.from_bed(f"{callability_dir}/{s['mask']}")
        hit = junctions_callable(mask, junctions, flank)
        n_samples += hit
        for k in np.flatnonzero(hit):
            inds[k].add(s["individual_id"])
            by_platform[k][s["platform"]].add(s["individual_id"])
    return n_samples, inds, by_platform
