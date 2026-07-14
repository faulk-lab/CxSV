# CxSV — Long-Read Complex Structural Variant Detection Pipeline

A Snakemake pipeline that detects complex structural variants (CxSVs) from
long-read sequencing data (PacBio / Oxford Nanopore) using an ensemble of
three callers: **Sniffles2**, **CuteSV**, and **SVIM**.

Clone this repo, drop BAMs into `bamfiles/<cohort>/`, and run one script.

The pipeline distinguishes two kinds of cohorts (see
[Public vs. Private CxSVs](#public-vs-private-cxsvs-cohorts) below):
**reference** cohorts (e.g. 1000 Genomes and other non-pathogenic samples)
build a pooled master list of common/benign CxSVs, and **query** cohorts
(e.g. patient samples) get filtered against that list, surfacing only CxSVs
*not* seen in the healthy-population reference — candidate pathogenic calls.

---

## Directory Structure

```
CxSV/
├── config/
│   ├── config.yaml          ← main pipeline config: cohorts, paths, params
│   ├── samples.txt           ← sample IDs (only needed for the optional 1000G download)
│   └── remote_1000g.yaml     ← optional: 1000G-ONT download settings
├── workflow/
│   ├── Snakefile
│   ├── rules/                ← rules split by pipeline phase
│   │   ├── indexing.smk
│   │   ├── calling.smk
│   │   ├── filtering_merging.smk
│   │   ├── cxsv_discovery.smk
│   │   ├── annotation.smk
│   │   ├── summary.smk
│   │   ├── plotting.smk
│   │   └── private_filter.smk    ← query cohorts only, see below
│   └── scripts/              ← scripts invoked by Snakemake rules
│       ├── vcf_to_bed.py, classify_cxsv.py, extract_cxsv_vcf.py,
│       │   annotate_cxsv.py, summarize_sv.py
│       ├── filter_vcf_by_bed.py, filter_table_by_bed.py  ← private-filter helpers
│       ├── plot_utils.R      ← shared VCF parsing / palette / theme helpers
│       └── plot_*.R          ← one script per figure
├── scripts/                  ← operational scripts (not called by Snakemake)
│   ├── run_pipeline.sh       ← runs the pipeline for one --cohort
│   ├── build_public_reference.py ← pools reference cohorts into the master CxSV list
│   ├── download_1000g_bams.sh← OPTIONAL: fetch 1000G-ONT BAMs
│   └── fetch_sample_list.py  ← OPTIONAL: discover 1000G-ONT sample IDs
├── bamfiles/
│   ├── 1000g/                 ← example reference cohort (gitignored)
│   └── pos_control/           ← example query cohort (gitignored)
├── resources/                 ← annotation BEDs (gitignored — see below)
├── hs1.fa, hs1.fa.fai          ← reference genome (gitignored — see below)
└── results/                   ← pipeline outputs (gitignored)
    ├── 1000g/, pos_control/, ...   ← per-cohort outputs (incl. logs/)
    └── _public_reference/          ← pooled master_cxsv.bed
```

---

## Quick Start

The pipeline assumes you run it **from the project root** — every path in
`config/config.yaml` is relative to it.

### 1. Register a cohort

Add it to `config/config.yaml → cohorts:`, tagged `role: reference` (healthy/
non-pathogenic samples, e.g. 1000 Genomes) or `role: query` (patient samples
to be filtered against the reference):

```yaml
cohorts:
  "1000g":        { role: reference }
  patient_batchA: { role: query }
```

### 2. Get BAMs into `bamfiles/<cohort>/`

Either:

- **Use your own long-read BAMs** — copy or symlink any `*.bam` file(s) into
  `bamfiles/<cohort>/` (coordinate-sorted). Any filename works — no naming
  scheme required. The sample ID used throughout the pipeline's output is
  just the filename with `.bam` stripped (e.g. `patient42.bam` → sample
  `patient42`), or
- **Download the 1000G-ONT test set** (optional, separate step — see
  [Optional: 1000G-ONT Download](#optional-1000g-ont-download) below).

### 3. Get the reference and annotation resources

- `hs1.fa` (CHM13v2.0 / T2T-CHM13v2) — download from the
  [T2T consortium](https://github.com/marbl/CHM13) and index with
  `samtools faidx hs1.fa`.
- `resources/*.bed` — see [Annotation Resources](#annotation-resources) below.

### 4. Run the pipeline, per cohort

```bash
conda activate snakemake_env
bash scripts/run_pipeline.sh --cohort 1000g
bash scripts/run_pipeline.sh --cohort patient_batchA
```

This dry-runs the DAG, then runs the full pipeline with `--keep-going`. To
resume an interrupted run:

```bash
bash scripts/run_pipeline.sh --cohort 1000g --resume
```

Or drive Snakemake directly (equivalent, since `workflow/Snakefile` is
auto-discovered from the project root):

```bash
snakemake --config active_cohort=1000g -n --cores 1     # dry run
snakemake --config active_cohort=1000g --cores 32       # full run
snakemake --config active_cohort=1000g --cores 1 \
  --cluster "sbatch -p short -c {threads} --mem=32G -t 4:00:00" \
  --jobs 20                                             # SLURM example
```

If a **query** cohort's results include private-CxSV outputs, you'll need
the reference master list built first — see the next section.

---

## Public vs. Private CxSVs (Cohorts)

Every cohort runs through the identical per-sample/per-cohort pipeline
(indexing → calling → filtering/merging → CxSV discovery → annotation →
summary → plots) into its own `results/<cohort>/`, completely independent of
every other cohort — adding a new patient batch never touches an earlier
one's results, and growing the reference set never forces a full patient
re-run.

**`role: reference`** cohorts (1000 Genomes, other non-pathogenic datasets)
are meant to be pooled into one master list of CxSVs expected in a healthy
population. Once you've run every reference cohort you want included:

```bash
python scripts/build_public_reference.py
```

This reads `config/config.yaml → cohorts:`, pools every `role: reference`
cohort's `results/<cohort>/final/cxsv_population.bed`, and merges them
(`sort` + `bedtools merge`) into `results/_public_reference/master_cxsv.bed`.
Re-run it any time you add a reference cohort or more samples to an existing
one — already-completed cohort runs are untouched.

**`role: query`** cohorts (patient samples) get three extra outputs once the
master list exists, from `workflow/rules/private_filter.smk`:

| File | Description |
|------|--------------|
| `final/private_cxsv_population.bed` | Cohort's CxSV loci **not** in the reference (candidate pathogenic) |
| `final/public_overlap_cxsv_population.bed` | Cohort's CxSV loci that **are** in the reference (likely benign; QC) |
| `final/private_cxsv_master_table.tsv` | Annotated master table, private subset only |
| `per_sample/{sample}_private_cxsv.vcf` | **Primary per-patient deliverable** — that sample's CxSV calls with anything matching the reference removed |

Matching uses a slop window, not exact overlap (`public_reference.match_window`
in `config.yaml`, default **1000 bp**) — CxSV loci are already clustered
regions (see `cxsv_params.cluster_distance`, 50 kb), so this only needs to
absorb boundary differences between independently-run cohorts, not merge
genuinely distinct loci. It sits between SURVIVOR's 300 bp per-call merge
tolerance (`survivor_params.max_dist`) and the much wider 50 kb CxSV
clustering window. Widen it if you're seeing likely-benign loci show up as
private just because of a few kb of breakpoint disagreement between cohorts;
narrow it if distinct nearby loci are getting incorrectly filtered out as
"seen in the reference."

```bash
# 1. Run every reference cohort
bash scripts/run_pipeline.sh --cohort 1000g

# 2. Build (or rebuild) the pooled master list
python scripts/build_public_reference.py

# 3. Run a query cohort — private-CxSV outputs are produced automatically
bash scripts/run_pipeline.sh --cohort pos_control
```

---

## Optional: 1000G-ONT Download

This is entirely separate from running the pipeline — it just populates
`bamfiles/<cohort>/` from the public 1000G-ONT S3 bucket. Skip it if you're
using your own BAMs.

```bash
python scripts/fetch_sample_list.py                 # writes config/samples.txt
bash scripts/download_1000g_bams.sh                 # downloads into bamfiles/1000g/
bash scripts/download_1000g_bams.sh --cohort other  # or into bamfiles/other/
bash scripts/run_pipeline.sh --cohort 1000g         # then run the pipeline as usual
```

Bucket/prefix/URL-template settings live in `config/remote_1000g.yaml` — edit
them there if the bucket layout changes. Estimated size: ~118 GB per BAM.

---

## Software Dependencies

| Tool      | Version tested | Installation               |
|-----------|-----------------|-----------------------------|
| Snakemake | ≥ 7.0           | `conda install snakemake`  |
| Sniffles2 | ≥ 2.2           | `conda install sniffles`   |
| CuteSV    | ≥ 2.1           | `conda install cutesv`     |
| SVIM      | ≥ 2.0           | `conda install svim`       |
| SURVIVOR  | ≥ 1.0.7         | Build from source (GitHub) |
| samtools  | ≥ 1.17          | `conda install samtools`   |
| bcftools  | ≥ 1.17          | `conda install bcftools`   |
| bedtools  | ≥ 2.31          | `conda install bedtools`   |
| Python    | ≥ 3.10          | via conda                  |
| R         | ≥ 4.3           | `conda install r-base`     |

### Required R packages

```r
install.packages(c(
  "ggplot2", "dplyr", "tidyr", "scales",
  "patchwork", "ggVennDiagram"
))
```

---

## Annotation Resources

`resources/` is gitignored — several of these files are hundreds of MB and
exceed GitHub's size limits, so they're fetched/regenerated locally instead
of committed. Paths are set in `config/config.yaml → annotation:`.

The pre-generated `resources/` files (already in the naming/format the
pipeline expects) are available for direct download here:
[arcticsynology.synology.me](https://arcticsynology.synology.me:5001/d/s/192nJPzqq5CpT5HWgZE0M7G3v7RR9FyA/VtNjg_kikWCKrrMjMe8bgTxjbr77P_uH-j7XAlKzMWA0).
Download and extract into `resources/` before running the pipeline.

Alternatively, regenerate them yourself from the original sources:

| Resource | Source | Format |
|----------|--------|--------|
| `genes.final.bed`   | UCSC Table Browser → knownGene → BED         | chr start end gene_name |
| `exons.filtered.bed`| UCSC Table Browser → knownGene exons          | chr start end gene_name |
| `repeats.final.bed` | UCSC rmsk table → BED                          | chr start end rep_class |
| `segdup.final.bed`  | UCSC genomicSuperDups → BED                    | chr start end |

> **Note:** All resources must use the same chromosome naming convention as
> the reference (e.g. `chr1`, not `1`, for hs1/CHM13).

---

## Output Files

Everything below lives under `results/<cohort>/` for the cohort you ran.

### Per-Sample (`per_sample/`)

| File                      | Description                       |
|---------------------------|------------------------------------|
| `{sample}_merged_sv.vcf`  | SURVIVOR merged long-read VCF     |
| `{sample}_cxsv_only.vcf`  | CxSV subset of merged VCF         |
| `{sample}_sv_summary.txt` | SV counts per caller + CxSV count |
| `{sample}_private_cxsv.vcf` | **query cohorts only** — `{sample}_cxsv_only.vcf` with anything matching the reference master list removed; see [Public vs. Private CxSVs](#public-vs-private-cxsvs-cohorts) |

### Final Call Sets (`final/`)

| File                     | Description                             |
|--------------------------|-------------------------------------------|
| `long_read_filtered.vcf` | Population-level merged long-read VCF   |
| `consensus.vcf`          | SVs supported by ≥2 callers             |
| `cxsv_only.vcf`          | CxSV-only VCF                           |
| `cxsv_master_table.tsv`  | 15-column annotated CxSV table          |
| `cxsv_population.bed`    | Population-level CxSV BED — for `reference` cohorts, this is what `build_public_reference.py` pools |
| `sv_counts.tsv`          | SV counts by type                       |
| `cxsv_count.txt`         | Total CxSV locus count                  |
| `private_cxsv_population.bed` | **query cohorts only** — cohort's CxSV loci not in the reference |
| `public_overlap_cxsv_population.bed` | **query cohorts only** — cohort's CxSV loci that are in the reference (QC) |
| `private_cxsv_master_table.tsv` | **query cohorts only** — annotated master table, private subset only |

### Pooled Reference (`results/_public_reference/`)

| File                | Description                                          |
|---------------------|-------------------------------------------------------|
| `master_cxsv.bed`   | Merged CxSV loci across every `role: reference` cohort, built by `scripts/build_public_reference.py` |

### Plots (`plots/`)

| File                           | Description                              |
|---------------------------------|-------------------------------------------|
| `sv_count_by_type.pdf`         | Bar plot: SV counts by type              |
| `sv_size_distribution.pdf`     | Violin + box: size on log scale          |
| `cxsv_venn.pdf`                | 3-way Venn: caller overlap               |
| `breakpoint_resolution.pdf`    | Density: CI width per caller             |
| `genomic_context.pdf`          | Bar: exonic/intronic/promoter/intergenic |
| `complexity_class_heatmap.pdf` | Bar + tile: complexity class + SV types  |
| `chromosome_density.pdf`       | Bar: SVs per Mb per chromosome           |

`workflow/scripts/plot_cxsv_population_sharing.R` is a standalone extra
script (not part of `rule all`) — run it manually once you have a per-locus
sample-sharing BED; see the comment at the top of that script.

---

## CxSV Definition (Computational)

A locus is classified as a CxSV if, within a **50 kb window**, it satisfies
any of five criteria (see `workflow/scripts/classify_cxsv.py`):

- **C1** — ≥ 3 unique breakpoints
- **C2** — Nested SV (one interval wholly inside another)
- **C3** — Overlapping SV types (physically share coordinates)
- **C4** — Copy-number shift (DEL/DUP) + orientation change (INV/BND/TRA)
- **C5** — Chromoanagenesis cluster: ≥ 3 SVs in a 10–50 kb window

CxSVs are never filtered by genomic context — intergenic, repeat, and segdup
loci are retained (most CxSVs are expected there).

## Complexity Class Labels

| Class Label | SV Types Present | Notes |
|-------------|-------------------|-------|
| `DEL+INV` | DEL + INV | |
| `DUP+INV` | DUP + INV | |
| `DEL+DUP+INV` | DEL + DUP + INV | |
| `INS+DEL` | INS + DEL | |
| `BND_cluster` | BND + ≥1 other | ≥3 breakpoints |
| `Multi_breakpoint_cluster` | Any | ≥5 breakpoints, <10 |
| `Chromothripsis_like` | Any | ≥10 breakpoints in window |
| `Complex_other` | ≥2 types, none of above | |

---

## Troubleshooting FAQ

### `active_cohort '<name>' is not registered in config['cohorts']`.

**Cause:** `--cohort <name>` (or `--config active_cohort=<name>`) doesn't
match anything under `config/config.yaml → cohorts:`. **Fix:** add it there
first, e.g. `patient_batchA: { role: query }`, or check for a typo — running
`bash scripts/run_pipeline.sh` with no `--cohort` prints the registered list.

### `glob_wildcards` finds no samples — pipeline exits immediately.

**Cause:** BAM files not in `bamfiles/<cohort>/` or wrong extension. **Fix:**

```bash
ls bamfiles/<cohort>/*.bam     # Should list your BAMs
```

### Sniffles2 fails with "ERROR: Index file not found".

**Cause:** The `.bam.bai` index was not created by `index_bam`, or the BAM is
not coordinate-sorted. **Fix:**

```bash
samtools sort -o bamfiles/<cohort>/SAMPLE_sorted.bam bamfiles/<cohort>/SAMPLE.bam
samtools index bamfiles/<cohort>/SAMPLE_sorted.bam
```

### CuteSV produces an empty VCF.

**Cause 1:** Working directory (`results/<cohort>/callers/SAMPLE_cutesv_work/`)
already exists with stale temp files from a previous failed run. **Fix:**

```bash
rm -rf results/<cohort>/callers/SAMPLE_cutesv_work/
bash scripts/run_pipeline.sh --cohort <cohort> --resume
```

**Cause 2:** Coverage too low — CuteSV requires ≥5–10× for reliable calls.
Check with `samtools coverage bamfiles/<cohort>/SAMPLE.bam`.

### SVIM outputs `variants.vcf` but needs a very high QUAL threshold.

**Cause:** SVIM scores are not probability-based; the default min QUAL is
appropriate for most datasets but may be too strict for low-coverage data.
**Fix:** Lower `filter_params.svim_min_qual` in `config/config.yaml` (e.g. to
3) for low-coverage samples (<20×).

### SURVIVOR merge fails with "Parsing error" or exits silently.

**Cause 1:** One of the input VCFs contains no records. SURVIVOR requires
≥1 record per input VCF. **Fix:**

```bash
bcftools view -H results/<cohort>/callers/SAMPLE_cutesv_filtered.vcf | wc -l
```

If 0, lower filter thresholds or check the caller logs.

**Cause 2:** VCF headers are malformed or INFO fields non-standard. **Fix:**

```bash
bcftools reheader --fai hs1.fa.fai -o fixed.vcf input.vcf
```

### `annotate_cxsv.py` raises "bedtools: command not found".

**Fix:**

```bash
conda activate cxsv_pipeline
which bedtools   # Should return a path
```

### Master table has "NA" for all gene overlap columns.

**Cause 1:** `genes_bed` / `exons_bed` paths in `config/config.yaml` point to
non-existent or empty files. **Fix:**

```bash
wc -l resources/genes.final.bed
head -1 resources/genes.final.bed   # Should show: chr1\tSTART\tEND\tGENE_NAME
```

**Cause 2:** Chromosome naming mismatch (`chr1` vs `1`). **Fix:**
`sed 's/^/chr/' resources/genes.bed > resources/genes_fixed.bed`, or strip
the prefix from the VCF depending on your reference.

### R scripts fail with "there is no package called 'ggVennDiagram'".

**Fix:**

```r
install.packages("ggVennDiagram")
```

`plot_cxsv_venn.R` falls back to the `VennDiagram` package automatically if
`ggVennDiagram` is unavailable.

### `patchwork` not found (complexity class plot fails).

**Fix:** `install.packages("patchwork")`

### Pipeline reruns rules that already completed successfully.

**Cause:** Snakemake timestamps changed (e.g. after copying files, NFS
mount). **Fix:**

```bash
snakemake --touch --config active_cohort=<cohort> --cores 1   # Mark up-to-date
snakemake --config active_cohort=<cohort> --cores 32          # Rerun only missing outputs
```

### CxSV count is 0 even though the pipeline completed.

**Cause:** Cluster parameters are too strict for your data. **Diagnosis:**
check `results/<cohort>/final/cxsv_summary.tsv` — are there rows with
`IsCxSV == NO`? **Fix:** relax `cxsv_params.min_breakpoints` in
`config/config.yaml`, then rerun `classify_cxsv` and downstream rules.

### A `query` cohort's private-CxSV rules fail: reference BED not found.

**Cause:** `results/_public_reference/master_cxsv.bed` doesn't exist yet.
**Fix:** run every `role: reference` cohort first, then build it:

```bash
bash scripts/run_pipeline.sh --cohort 1000g
python scripts/build_public_reference.py
```

### `build_public_reference.py` exits with "No cohorts with role 'reference'".

**Cause:** No cohort in `config/config.yaml → cohorts:` has `role: reference`
(or the ones that do haven't produced `final/cxsv_population.bed` yet — the
script lists which are `[ok]` vs `[missing]`). **Fix:** register at least one
reference cohort and run it to completion first.

---

## Citation

If you use this pipeline, please cite the underlying tools:

- **Sniffles2:** Smolka et al., *Nature Methods* 2024
- **CuteSV:** Jiang et al., *Genome Biology* 2020
- **SVIM:** Heller & Vingron, *Bioinformatics* 2019
- **SURVIVOR:** Jeffares et al., *Nature Communications* 2017
- **bedtools:** Quinlan & Hall, *Bioinformatics* 2010
