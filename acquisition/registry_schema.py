"""Schema and helpers for the CxSV background-genome registry.

The registry is one TSV at ``<datasets>/registry/sample_registry.tsv`` holding one
row per *file* that belongs to the background set.  It is the single source of
truth for what has been acquired; every acquisition script appends to it and
every consumer (the Snakemake pipeline, QC summaries, reports) reads from it
rather than globbing the filesystem.

Design notes
------------
* One row per file, not per sample.  A sample may contribute a FASTQ row and
  later an alignment row; they are tied together by ``individual_id``.
* ``individual_id`` is the normalised Coriell/1KGP identifier (``NA#####`` /
  ``HG#####``).  The 1kGP-LRSC release uses ``GM#####`` for the same people, so
  it is normalised to ``NA#####`` on ingest - without this the same individual
  is silently counted twice in background frequency.
* ``coverage_x`` is always accompanied by ``coverage_source`` so a measured
  depth is never confused with an estimate derived from raw yield.
* ``status`` drives idempotent re-runs: only ``verified`` rows are considered
  complete.
"""

from __future__ import annotations

# Column order is load-bearing: readers index by name, but a stable order keeps
# the file diffable in git.
COLUMNS: list[str] = [
    "file_id",           # stable key: <cohort>/<basename>
    "individual_id",     # normalised NA#####/HG#####
    "source_sample_id",  # identifier as published (e.g. GM18501)
    "cohort",            # dataset subdirectory name
    "accession",         # ENA/SRA run or bucket key
    "platform",          # ONT | PACBIO
    "chemistry",         # R9/R10/Revio/SequelII, or "" if unknown
    "library",           # LSK110 / LSK114 / HiFi ...
    "basecaller",        # e.g. dorado-0.5.0-sup, guppy-sup, "" if unknown
    "data_type",         # fastq | unaligned_bam | aligned_bam | aligned_cram
    "reference",         # T2T-CHM13v2 | GRCh38 | none
    "coverage_x",        # numeric mean coverage
    "coverage_source",   # how coverage_x was obtained - never leave blank
    "rel_path",          # path relative to the datasets root
    "size_bytes",
    "md5",               # published or computed checksum, "" if none
    "status",            # pending | downloading | downloaded | verified | failed | excluded
    "date_added",        # ISO date
    "notes",
]

#: Values of ``status`` that mean the file is present and checked.
COMPLETE_STATUSES = frozenset({"verified"})

#: Values of ``coverage_source`` in use, with what each one actually measures.
COVERAGE_SOURCES = {
    "cramino": "mean coverage reported by cramino on the published alignment",
    "measured_chr20_depth": "samtools coverage meandepth over chr20:10-11 Mb",
    "mosdepth": "mosdepth mean over the whole genome",
    "ena_base_count": "ENA base_count / 3.1 Gb - raw yield, an upper bound on aligned depth",
    "hprc_index": "coverage column of the HPRC release-2 sequencing index",
}


def normalise_individual_id(raw: str) -> str:
    """Return the canonical 1KGP identifier for ``raw``.

    The 1kGP-LRSC release labels Coriell samples ``GM#####`` while ENA and the
    1KGP metadata use ``NA#####``.  Collapsing them is what stops one person
    being counted as two independent background observations.

    >>> normalise_individual_id("GM18501-ONT-chm13-R9.bam")
    'NA18501'
    >>> normalise_individual_id("HG00105_PMDV_FINAL.haplotagged.bam")
    'HG00105'
    """
    import re

    # NB: do not use \b here. Filenames separate fields with underscores
    # ("HG00105_PMDV_FINAL..."), and "_" is a word character, so \b never
    # matches between the digits and the underscore. Use explicit
    # non-alphanumeric lookarounds instead.
    m = re.search(r"(?<![A-Za-z0-9])(NA|GM|HG)(\d{5})(?!\d)", raw)
    if not m:
        raise ValueError(f"no 1KGP-style identifier found in {raw!r}")
    prefix, digits = m.group(1), m.group(2)
    return ("NA" if prefix == "GM" else prefix) + digits


def header_line() -> str:
    """The registry TSV header row, newline-terminated."""
    return "\t".join(COLUMNS) + "\n"
