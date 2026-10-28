#!/usr/bin/env python3
"""Build a download manifest for an ENA study of long-read human WGS.

Writes one row per file to download, with the publisher's md5 so the transfer
can be verified rather than assumed complete.

Usage:
    python3 build_ena_manifest.py PRJEB89727 --out manifest_PRJEB89727.tsv

Why the coverage column is an *estimate*: ENA reports ``base_count``, the number
of bases submitted, so ``base_count / 3.1e9`` is raw yield per haploid genome
size - an upper bound on the aligned mean depth, since some reads will not map.
It is recorded as ``ena_base_count`` in ``coverage_source`` so it is never
mistaken for a measured depth. Real depth is computed after alignment.
"""

from __future__ import annotations

import argparse
import csv
import sys
import urllib.parse
import urllib.request

ENA_PORTAL = "https://www.ebi.ac.uk/ena/portal/api/filereport"

#: Haploid human genome size used to turn yield into approximate fold-coverage.
GENOME_SIZE_BP = 3.1e9

FIELDS = [
    "run_accession",
    "sample_accession",
    "sample_title",
    "instrument_model",
    "library_strategy",
    "base_count",
    "read_count",
    "submitted_ftp",
    "submitted_md5",
    "submitted_bytes",
    "submitted_format",
    "fastq_ftp",
    "fastq_md5",
    "fastq_bytes",
]


def fetch_filereport(accession: str) -> list[dict]:
    """Return the ENA filereport for ``accession`` as a list of row dicts."""
    query = urllib.parse.urlencode(
        {
            "accession": accession,
            "result": "read_run",
            "fields": ",".join(FIELDS),
            "format": "tsv",
            "limit": 0,
        }
    )
    with urllib.request.urlopen(f"{ENA_PORTAL}?{query}", timeout=300) as fh:
        text = fh.read().decode()
    lines = text.rstrip("\n").split("\n")
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:]]


def normalise_individual_id(raw: str) -> str:
    """Canonical 1KGP identifier; GM##### and NA##### are the same person."""
    import re

    m = re.search(r"(?<![A-Za-z0-9])(NA|GM|HG)(\d{5})(?!\d)", raw)
    if not m:
        return ""
    prefix, digits = m.group(1), m.group(2)
    return ("NA" if prefix == "GM" else prefix) + digits


def build(accession: str, coverage_floor: float) -> list[dict]:
    """Expand the filereport into one output row per downloadable file."""
    rows = fetch_filereport(accession)
    if not rows:
        raise SystemExit(f"ENA returned no runs for {accession}")

    # Prefer the submitter's own files (FASTQ or unaligned BAM as deposited);
    # fall back to ENA's generated FASTQ only when no submitted file exists.
    out: list[dict] = []
    for r in rows:
        urls = [u for u in r["submitted_ftp"].split(";") if u]
        md5s = [m for m in r["submitted_md5"].split(";") if m]
        sizes = [s for s in r["submitted_bytes"].split(";") if s]
        kind = "submitted"
        if not urls:
            urls = [u for u in r["fastq_ftp"].split(";") if u]
            md5s = [m for m in r["fastq_md5"].split(";") if m]
            sizes = [s for s in r["fastq_bytes"].split(";") if s]
            kind = "ena_generated"
        if not urls:
            print(f"  WARNING: no files listed for {r['run_accession']}", file=sys.stderr)
            continue

        bases = int(r["base_count"]) if r["base_count"] else 0
        cov = bases / GENOME_SIZE_BP
        for url, md5, size in zip(urls, md5s, sizes):
            out.append(
                {
                    "run_accession": r["run_accession"],
                    "individual_id": normalise_individual_id(r["sample_title"]),
                    "sample_title": r["sample_title"],
                    "filename": url.rsplit("/", 1)[-1],
                    "url": "https://" + url if not url.startswith("http") else url,
                    "md5": md5,
                    "size_bytes": size,
                    "file_source": kind,
                    "format": r["submitted_format"] or "FASTQ",
                    "instrument_model": r["instrument_model"],
                    "base_count": bases,
                    "est_coverage_x": round(cov, 2),
                    "coverage_source": "ena_base_count",
                    "passes_coverage_floor": "yes" if cov >= coverage_floor else "no",
                }
            )
    return out


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("accession", help="ENA study accession, e.g. PRJEB89727")
    p.add_argument("--out", required=True, help="output manifest TSV")
    p.add_argument(
        "--coverage-floor",
        type=float,
        default=17.0,
        help="fold-coverage below which a sample is flagged (default 17, not dropped)",
    )
    a = p.parse_args()

    rows = build(a.accession, a.coverage_floor)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    total = sum(int(r["size_bytes"]) for r in rows if r["size_bytes"])
    passing = sum(1 for r in rows if r["passes_coverage_floor"] == "yes")
    indiv = {r["individual_id"] for r in rows if r["individual_id"]}
    missing_md5 = sum(1 for r in rows if not r["md5"])
    print(
        f"{a.accession}: {len(rows)} files, {len(indiv)} individuals, "
        f"{total / 1e12:.2f} TB, {passing} files >= {a.coverage_floor}x, "
        f"{missing_md5} without md5 -> {a.out}"
    )


if __name__ == "__main__":
    main()
