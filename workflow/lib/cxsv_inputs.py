"""Sample discovery shared by the panel and patient workflows.

A *cohort* (panel) or *run* (patient) is one named set of alignments, given by
one spec dict from config. Two input modes:

Directory mode (default)
    Every ``*.bam`` / ``*.cram`` in ``<bamfiles_root>/<name>/``. Sample ID is the
    filename without its extension. A missing index is built next to the file by
    ``rules/indexing.smk``. ``individual_id`` defaults to the sample ID and
    ``platform`` to the spec's ``platform`` (default ``unknown``).

Sample-sheet mode   ``samples: <sheet.tsv>``
    Tab-separated with a header. Required: ``sample_id``, ``alignment``.
    Optional: ``index``, ``individual_id``, ``platform``, ``coverage_x``; other
    columns are ignored. Files may live anywhere (e.g. a read-only datasets
    tree); their indexes must already exist, because the pipeline never writes
    into a sheet's data directories. ``min_coverage: <x>`` keeps rows with
    ``coverage_x >= x``.

The functions here are plain Python (no Snakemake objects) so they can be
unit-tested; ``rules/common.smk`` calls them at DAG-build time.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

INDEX_EXT = {".bam": ".bai", ".cram": ".crai"}


@dataclass(frozen=True)
class SampleRecord:
    cohort: str
    sample_id: str
    alignment: str
    index: str
    individual_id: str
    platform: str
    coverage_x: str = ""

    @property
    def is_cram(self) -> bool:
        return self.alignment.endswith(".cram")


class SampleSheetError(ValueError):
    """A sample sheet or directory cannot be turned into a valid sample list."""


def default_index(alignment: str) -> str:
    return alignment + INDEX_EXT[Path(alignment).suffix]


def samples_from_dir(cohort: str, data_dir: str | Path, platform: str = "unknown") -> list[SampleRecord]:
    found: dict[str, SampleRecord] = {}
    for ext in INDEX_EXT:
        for p in sorted(Path(data_dir).glob(f"*{ext}")):
            sid = p.name[: -len(ext)]
            if sid in found:
                raise SampleSheetError(
                    f"{cohort}: sample '{sid}' appears as both {found[sid].alignment} and {p} "
                    f"in {data_dir}/ -- keep one alignment per sample")
            found[sid] = SampleRecord(cohort, sid, str(p), default_index(str(p)), sid, platform)
    return [found[k] for k in sorted(found)]


def _read_sheet(path: str | Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader((ln.replace("\r", "") for ln in fh), delimiter="\t"))


def samples_from_sheet(cohort: str, sheet: str | Path, min_coverage=None,
                       platform: str = "unknown", check_files: bool = True) -> tuple[list[SampleRecord], int]:
    """Return (records, n_dropped_below_min_coverage). Raises SampleSheetError
    listing every problem found, so a bad sheet is fixed in one pass."""
    rows = _read_sheet(sheet)
    if not rows:
        raise SampleSheetError(f"{cohort}: sample sheet {sheet} has no rows")
    for col in ("sample_id", "alignment"):
        if col not in rows[0]:
            raise SampleSheetError(f"{cohort}: sample sheet {sheet} lacks required column '{col}'")
    if min_coverage is not None and "coverage_x" not in rows[0]:
        raise SampleSheetError(f"{cohort}: min_coverage is set but {sheet} has no coverage_x column")

    out: dict[str, SampleRecord] = {}
    problems: list[str] = []
    dropped = 0
    for i, row in enumerate(rows, start=2):
        sid, aln = (row.get("sample_id") or "").strip(), (row.get("alignment") or "").strip()
        cov = (row.get("coverage_x") or "").strip()
        if min_coverage is not None:
            try:
                if float(cov) < float(min_coverage):
                    dropped += 1
                    continue
            except ValueError:
                problems.append(f"line {i}: {sid}: coverage_x '{cov}' is not a number")
                continue
        if not sid:
            problems.append(f"line {i}: empty sample_id")
            continue
        if Path(aln).suffix not in INDEX_EXT:
            problems.append(f"line {i}: {sid}: alignment must end in .bam or .cram: {aln}")
            continue
        idx = (row.get("index") or "").strip() or default_index(aln)
        if sid in out:
            problems.append(f"line {i}: duplicate sample_id '{sid}'")
        if check_files and not Path(aln).exists():
            problems.append(f"line {i}: {sid}: alignment not found: {aln}")
        if check_files and not Path(idx).exists():
            problems.append(f"line {i}: {sid}: index not found: {idx}")
        out[sid] = SampleRecord(
            cohort, sid, aln, idx,
            (row.get("individual_id") or "").strip() or sid,
            (row.get("platform") or "").strip() or platform,
            cov)
    if problems:
        shown = "\n  ".join(problems[:20])
        more = f"\n  ... and {len(problems) - 20} more" if len(problems) > 20 else ""
        raise SampleSheetError(f"{cohort}: sample sheet {sheet} has {len(problems)} problem(s):\n  {shown}{more}")
    return [out[k] for k in sorted(out)], dropped


def build_sample_table(cohort_specs: dict, bamfiles_root: str, check_files: bool = True):
    """{(cohort, sample_id): SampleRecord} plus one log line per cohort.

    ``check_files=False`` skips existence checks (used by tests and dry docs)."""
    table: dict[tuple[str, str], SampleRecord] = {}
    notes: list[str] = []
    for cohort, spec in cohort_specs.items():
        spec = spec or {}
        platform = spec.get("platform", "unknown")
        if spec.get("samples"):
            recs, dropped = samples_from_sheet(cohort, spec["samples"], spec.get("min_coverage"),
                                               platform, check_files)
            src = f"sheet {spec['samples']}"
            if dropped:
                src += f"; {dropped} below min_coverage={spec['min_coverage']} left out"
        else:
            if spec.get("min_coverage") is not None:
                raise SampleSheetError(f"{cohort}: min_coverage needs a sample sheet (samples:)")
            d = Path(bamfiles_root) / cohort
            recs = samples_from_dir(cohort, d, platform)
            src = f"{d}/"
        for r in recs:
            table[(cohort, r.sample_id)] = r
        n_cram = sum(r.is_cram for r in recs)
        notes.append(f"{cohort}: {len(recs)} samples ({n_cram} CRAM, {len(recs) - n_cram} BAM) from {src}")
    return table, notes
