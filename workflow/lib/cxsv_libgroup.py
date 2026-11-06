"""Library-QC calibration group of a patient alignment (docs/LIBRARY_QC.md section 4b).

The library-QC flag compares a genome with the panel genomes of the same platform
group (config/library_qc.yaml `groups`). The group is chosen, in this order:

  1. sample sheet       the `qc_group` column (or patients.<sample>.qc_group): always wins
  2. header             alignment header: @RG PL / DS / ID, @PG CL / DS, @CO, read-group
                        tags of the sampled reads (dorado writes the basecall model there)
                          PacBio      PL:PACBIO, DS READTYPE=CCS, @PG ccs / jasmine / pbmm2
                          ONT R10     r10.4 / R10 / e8.2 / LSK114 / NBD114 / RBK114 in the text
                          ONT R9      r9.4 / R9 / LSK109 / LSK110 / NBD104 / RBK004 in the text
                        (both R9 and R10 cues = ambiguous, the header is not used)
  3. reads              the reads the QC samples anyway (primary, MAPQ >= 20, in the QC
                        windows; the first `max_reads`):
                          PacBio      read names <movie>/<zmw>/ccs or an `rq` tag
                          1kGP-1019   base qualities uniform within each read (SRA-lite
                                      style; the 1kGP ONT-1019 panel cohort has them)
                          ONT R10     median per-read mean base quality >= r10_min_median_q (30)
                          ONT R9      median per-read mean base quality <= r9_max_median_q (28)
                        (panel, 3,000 reads per genome: R9 20.6-25.6, R10 31.8-45.9, PacBio
                        38.2-39.7, 1kGP-1019 uniform Q 33; docs/library_qc/group_detection_panel.tsv)
  4. platform default   library_qc.default_group[<sample-sheet platform>]

Header and reads are both evaluated and reported; when both decide and disagree, the
header is used and the disagreement is written to the evidence. Pure functions: the
reads are summarised by library_qc.py, which owns pysam.
"""
from __future__ import annotations

import re
import statistics

PACBIO, R10, R9, ONT1019 = "PacBio_HiFi", "ONT_LRSC_R10", "ONT_LRSC_R9", "ONT_1kGP_1019"

DEFAULTS = {"max_reads": 20000, "min_reads": 200, "r10_min_median_q": 30.0, "r9_max_median_q": 28.0,
            "uniform_q_min_fraction": 0.9, "pacbio_min_fraction": 0.5}

_R10 = re.compile(r"(?<![a-z0-9])(r10(\.4(\.1)?)?|r1041)(?![0-9a-z])|e8\.2|lsk114|nbd114|rbk114", re.I)
_R9 = re.compile(r"(?<![a-z0-9])(r9(\.4(\.1)?)?|r941)(?![0-9a-z])|lsk109|lsk110|nbd104|rbk004", re.I)
_PB_PG = re.compile(r"^(ccs|jasmine|pbmm2|lima|pbtrim)", re.I)
_CCS_NAME = re.compile(r"^m\d+[A-Za-z]?_\d+_\d+(_s\d+)?/\d+/(ccs|\d+_\d+)")


def header_text(header: dict) -> list[tuple[str, str]]:
    """(where, text) pairs from a pysam-style header dict that may name the chemistry."""
    out = []
    for rg in header.get("RG", []) or []:
        for k in ("PL", "DS", "ID", "PU", "PM"):
            if rg.get(k):
                out.append((f"@RG {k}", str(rg[k])))
    for pg in header.get("PG", []) or []:
        for k in ("ID", "PN", "CL", "DS"):
            if pg.get(k):
                out.append((f"@PG {k}", str(pg[k])))
    for co in header.get("CO", []) or []:
        out.append(("@CO", str(co)))
    return out


def clean(t: str, n: int = 60) -> str:
    """Header text made safe for one TSV cell."""
    return re.sub(r"\s+", " ", str(t))[:n]


def classify_text(pairs: list[tuple[str, str]], what: str) -> tuple[str | None, str]:
    """Group from header-like text, with the evidence; None when it does not decide."""
    pairs = [(w, clean(t, 400)) for w, t in pairs]
    for where, t in pairs:
        if (where == "@RG PL" and t.upper().startswith("PACBIO")) or "READTYPE=CCS" in t.upper() \
                or (where in ("@PG ID", "@PG PN") and _PB_PG.match(t)):
            return PACBIO, f"{what}: {where} '{clean(t, 40)}'"
    r10 = [(w, m.group(0)) for w, t in pairs for m in [_R10.search(t)] if m]
    r9 = [(w, m.group(0)) for w, t in pairs for m in [_R9.search(t)] if m]
    if r10 and r9:
        return None, f"{what}: ambiguous (R10 cue {r10[0][1]!r} in {r10[0][0]}, R9 cue {r9[0][1]!r} in {r9[0][0]})"
    if r10:
        return R10, f"{what}: {r10[0][0]} names {r10[0][1]!r}"
    if r9:
        return R9, f"{what}: {r9[0][0]} names {r9[0][1]!r}"
    pl = sorted({t.upper() for w, t in pairs if w == "@RG PL"})
    return None, f"{what}: no chemistry cue" + (f" (PL {','.join(pl)})" if pl else " (no @RG)")


class ReadCues:
    """Running summary of the sampled reads (fed by library_qc.py)."""

    def __init__(self, max_reads: int = DEFAULTS["max_reads"]):
        self.max_reads, self.n, self.mean_q, self.uniform, self.no_qual = max_reads, 0, [], 0, 0
        self.pacbio, self.rg_tags = 0, {}

    def add(self, name: str, quals, rg=None, has_rq: bool = False) -> None:
        """One sampled read: its name, base qualities (array or None), RG tag, rq tag present."""
        if self.n >= self.max_reads:
            return
        self.n += 1
        if has_rq or _CCS_NAME.match(name or ""):
            self.pacbio += 1
        if rg is not None and len(self.rg_tags) < 50:
            self.rg_tags[str(rg)] = self.rg_tags.get(str(rg), 0) + 1
        if quals is None or len(quals) == 0:
            self.no_qual += 1
            return
        q = list(quals)
        self.mean_q.append(sum(q) / len(q))
        self.uniform += min(q) == max(q)

    def summary(self) -> dict:
        nq = len(self.mean_q)
        return {"reads": self.n, "median_mean_q": statistics.median(self.mean_q) if nq else float("nan"),
                "uniform_q_fraction": self.uniform / nq if nq else float("nan"),
                "pacbio_fraction": self.pacbio / self.n if self.n else float("nan"),
                "read_groups": sorted(self.rg_tags)}


def classify_reads(s: dict, rules: dict) -> tuple[str | None, str]:
    """Group from the read summary (ReadCues.summary), with the evidence; None when undecided."""
    r = {**DEFAULTS, **(rules or {})}
    if not s.get("reads") or s["reads"] < r["min_reads"]:
        return None, f"reads: too few sampled ({s.get('reads', 0)} < {r['min_reads']})"
    rg_group, rg_how = classify_text([("read RG tag", t) for t in s.get("read_groups", [])], "reads")
    if rg_group:
        return rg_group, rg_how
    if s["pacbio_fraction"] >= r["pacbio_min_fraction"]:
        return PACBIO, f"reads: {s['pacbio_fraction']:.0%} PacBio CCS read names or rq tags"
    mq, uq = s["median_mean_q"], s["uniform_q_fraction"]
    if mq != mq:
        return None, "reads: no base qualities"
    if uq >= r["uniform_q_min_fraction"]:
        return ONT1019, f"reads: base qualities uniform within {uq:.0%} of reads (SRA-lite style, as the 1kGP ONT-1019 cohort)"
    if mq >= r["r10_min_median_q"]:
        return R10, f"reads: median per-read mean base quality {mq:.1f} >= {r['r10_min_median_q']:g} (R10-like)"
    if mq <= r["r9_max_median_q"]:
        return R9, f"reads: median per-read mean base quality {mq:.1f} <= {r['r9_max_median_q']:g} (R9-like)"
    return None, f"reads: median per-read mean base quality {mq:.1f} between the R9 and R10 bounds"


def resolve(override: str | None, header_result, read_result, platform: str | None,
            default_group: dict | None, groups: dict | None = None) -> dict:
    """The group to use, how it was chosen, and the evidence of every step."""
    hg, hhow = header_result
    rg, rhow = read_result
    ev = [hhow, rhow]
    if hg and rg and hg != rg:
        ev.append(f"header ({hg}) and reads ({rg}) disagree; header used")
    plat = (platform or "").upper()
    det = hg or rg
    if det and plat and ((det == PACBIO) != (plat.startswith("PACBIO"))):
        ev.append(f"detected {det} but the sample sheet says platform {plat}")
    if override:
        g, src = override, "sample sheet qc_group"
        if det and det != override:
            ev.append(f"auto-detection would give {det}")
    elif hg:
        g, src = hg, "header"
    elif rg:
        g, src = rg, "reads"
    else:
        g, src = (default_group or {}).get(plat, "NA"), f"platform default ({plat or 'no platform'})"
    if groups is not None and g not in groups:
        ev.append(f"group {g} has no calibration")
    return {"qc_group": g, "qc_group_source": src, "qc_group_detected": det or "NA",
            "qc_group_evidence": "; ".join(ev)}
