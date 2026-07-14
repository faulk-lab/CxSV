#!/usr/bin/env python3
"""
build_public_reference.py
===========================
Pool every "reference"-role cohort's population-level CxSV BED
(results/<cohort>/final/cxsv_population.bed) into one merged master list at
results/_public_reference/master_cxsv.bed — the benign/common-variant
reference that query (patient) cohorts get filtered against
(workflow/rules/private_filter.smk).

This is a standalone operational script, not a Snakemake rule: it reads
across however many independent cohort result trees are registered with
role: reference in config/config.yaml, which doesn't fit a single cohort's
DAG. Re-run it any time a reference cohort is added, or an existing one gets
more samples — cohorts already completed are untouched; only the merged
master list is rebuilt.

Usage
-----
    python scripts/build_public_reference.py
    python scripts/build_public_reference.py --config config/config.yaml
"""

import argparse
import subprocess
import sys
from pathlib import Path

import yaml


def load_reference_cohorts(config_path: str) -> tuple[list[str], str, str]:
    """Return (reference_cohort_names, results_root, master_bed_path) from config.yaml."""
    with open(config_path) as fh:
        cfg = yaml.safe_load(fh)

    cohorts = cfg.get("cohorts", {})
    reference_cohorts = [name for name, spec in cohorts.items()
                         if (spec or {}).get("role") == "reference"]

    results_root = cfg.get("results_root", "results")
    master_bed   = cfg.get("public_reference", {}).get(
        "master_bed", "results/_public_reference/master_cxsv.bed"
    )
    return reference_cohorts, results_root, master_bed


def main():
    parser = argparse.ArgumentParser(
        description="Pool all reference-cohort CxSV BEDs into one master list."
    )
    parser.add_argument("--config", default="config/config.yaml",
                        help="Path to config.yaml (default: config/config.yaml)")
    args = parser.parse_args()

    reference_cohorts, results_root, master_bed_path = load_reference_cohorts(args.config)

    if not reference_cohorts:
        sys.exit(
            "ERROR: No cohorts with role 'reference' found in "
            f"{args.config} -> cohorts:.\n"
            "Register at least one, e.g.:\n"
            "  cohorts:\n"
            '    "1000g": { role: reference }\n'
        )

    print(f"Reference cohorts registered: {', '.join(reference_cohorts)}")

    found_beds = []
    for cohort in reference_cohorts:
        bed = Path(results_root) / cohort / "final" / "cxsv_population.bed"
        if bed.is_file():
            n_loci = sum(1 for _ in bed.open())
            print(f"  [ok]      {cohort}: {bed}  ({n_loci} loci)")
            found_beds.append(bed)
        else:
            print(f"  [missing] {cohort}: {bed} not found — run it first "
                  f"(bash scripts/run_pipeline.sh --cohort {cohort}); skipping for now")

    if not found_beds:
        sys.exit(
            "\nERROR: None of the registered reference cohorts have been run yet.\n"
            "Run each with: bash scripts/run_pipeline.sh --cohort <name>"
        )

    master_bed = Path(master_bed_path)
    master_bed.parent.mkdir(parents=True, exist_ok=True)

    n_pooled = 0
    pooled_path = master_bed.parent / ".pooled.tmp.bed"
    with pooled_path.open("w") as pooled:
        for bed in found_beds:
            with bed.open() as fh:
                for line in fh:
                    if line.strip():
                        pooled.write(line)
                        n_pooled += 1

    sort_proc = subprocess.run(
        ["sort", "-k1,1", "-k2,2n", str(pooled_path)],
        capture_output=True, text=True, check=True
    )
    sorted_path = master_bed.parent / ".sorted.tmp.bed"
    sorted_path.write_text(sort_proc.stdout)

    with master_bed.open("w") as out:
        subprocess.run(
            ["bedtools", "merge", "-i", str(sorted_path), "-c", "4", "-o", "distinct"],
            stdout=out, check=True
        )

    n_merged = sum(1 for _ in master_bed.open())
    pooled_path.unlink()
    sorted_path.unlink()

    print(f"\nPooled {n_pooled} loci from {len(found_beds)}/{len(reference_cohorts)} "
          f"reference cohort(s) -> {n_merged} merged loci.")
    print(f"Master reference list: {master_bed}")


if __name__ == "__main__":
    main()
