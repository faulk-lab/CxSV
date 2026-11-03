"""Shared helpers for the graded-privacy calibration (roadmap item 10).

The calibration scores every panel genome's own CxSV loci with the patient
scorer against the panel with that genome's individual (leave-one-individual-
out, LOIO) or its whole family (leave-one-family-out, LOFO) removed. Nothing
here changes the panel release: the left-out individuals are subtracted from
the carrier sets of the released panel_junction_carriers.tsv at query time
(`ExcludingIndex`), and from the per-locus callability counts.
"""
from __future__ import annotations

import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO, "workflow", "lib"))

import cxsv_matching as cm  # noqa: E402

MODES = ("LOIO", "LOFO")


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, rows, cols):
    tmp = path + ".tmp"
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    os.replace(tmp, path)


def load_families(path):
    """{individual_id: row} from families.tsv (families.py)."""
    return {r["individual_id"]: r for r in read_tsv(path)}


def excluded(individual, families, mode):
    """Panel individuals removed when scoring `individual` in `mode`."""
    if mode == "LOIO":
        return {individual}
    if mode == "LOFO":
        fam = families[individual]["panel_family_members"]
        return set(fam.split(",")) | {individual}
    if mode == "NONE":
        return set()
    raise ValueError(mode)


class ExcludingIndex:
    """A cxsv_matching.SiteIndex whose sites lose the excluded carriers.

    query() returns the same matched sites as the wrapped index (the coordinate
    envelopes of the released panel are kept) with carrier_set minus the
    excluded individuals, so cxsv_matching.junction_carriers and
    classify_locus run unchanged. A site left with no carriers still matches
    but contributes no carrier, which is what removing it would do to the
    status (validated against true rebuilds in validate_rebuild.py)."""

    def __init__(self, index, exclude):
        self.index = index
        self.exclude = frozenset(exclude)

    def query(self, j, tolerance, ins_len_ratio):
        hits = self.index.query(j, tolerance, ins_len_ratio)
        if not self.exclude:
            return hits
        out = []
        for s in hits:
            if s["carrier_set"] & self.exclude:
                s = dict(s)
                s["carrier_set"] = s["carrier_set"] - self.exclude
            out.append(s)
        return out


def locus_id(row):
    return f"{row['Chr']}:{row['Start']}-{row['End']}"
