# CxSV — Long-Read Complex Structural Variant Detection Pipeline

A Snakemake pipeline that detects complex structural variants (CxSVs) from
long-read sequencing data (PacBio / Oxford Nanopore) using an ensemble of
three callers: **Sniffles2**, **CuteSV**, and **SVIM**.

Clone this repo, drop BAMs into `bamfiles/`, and run one script.

---

## Directory Structure

```
CxSV/
├── config/
│   ├── config.yaml          ← main pipeline config (relative paths, see below)
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
│   │   └── plotting.smk
│   └── scripts/              ← scripts invoked by Snakemake rules
│       ├── vcf_to_bed.py, classify_cxsv.py, extract_cxsv_vcf.py,
│       │   annotate_cxsv.py, summarize_sv.py
│       ├── plot_utils.R      ← shared VCF parsing / palette / theme helpers
│       └── plot_*.R          ← one script per figure
├── scripts/                  ← operational scripts (not called by Snakemake)
│   ├── run_pipeline.sh       ← runs the pipeline against bamfiles/
│   ├── download_1000g_bams.sh← OPTIONAL: fetch 1000G-ONT BAMs
│   └── fetch_sample_list.py  ← OPTIONAL: discover 1000G-ONT sample IDs
├── bamfiles/                  ← put your BAMs here (gitignored)
├── resources/                 ← annotation BEDs (gitignored — see below)
├── hs1.fa, hs1.fa.fai          ← reference genome (gitignored — see below)
├── results/                   ← pipeline outputs (gitignored)
└── logs/                      ← per-rule logs (gitignored)
```

---

## Quick Start

The pipeline assumes you run it **from the project root** — every path in
`config/config.yaml` is relative to it.

### 1. Get BAMs into `bamfiles/`

Either:

- **Use your own long-read BAMs** — copy or symlink any `*.bam` file(s) in
  (coordinate-sorted). Any filename works — no naming scheme required. The
  sample ID used throughout the pipeline's output is just the filename with
  `.bam` stripped (e.g. `patient42.bam` → sample `patient42`), or
- **Download the 1000G-ONT test set** (optional, separate step — see
  [Optional: 1000G-ONT Download](#optional-1000g-ont-download) below).

### 2. Get the reference and annotation resources

- `hs1.fa` (CHM13v2.0 / T2T-CHM13v2) — download from the
  [T2T consortium](https://github.com/marbl/CHM13) and index with
  `samtools faidx hs1.fa`.
- `resources/*.bed` — see [Annotation Resources](#annotation-resources) below.

### 3. Run the pipeline

```bash
conda activate snakemake_env
bash scripts/run_pipeline.sh
```

This dry-runs the DAG, then runs the full pipeline with `--keep-going`. To
resume an interrupted run:

```bash
bash scripts/run_pipeline.sh --resume
```

Or drive Snakemake directly (equivalent, since `workflow/Snakefile` is
auto-discovered from the project root):

```bash
snakemake -n --cores 1                    # dry run
snakemake --cores 32                      # full run
snakemake --cores 1 \
  --cluster "sbatch -p short -c {threads} --mem=32G -t 4:00:00" \
  --jobs 20                               # SLURM example
```

---

## Optional: 1000G-ONT Download

This is entirely separate from running the pipeline — it just populates
`bamfiles/` from the public 1000G-ONT S3 bucket. Skip it if you're using
your own BAMs.

```bash
python scripts/fetch_sample_list.py      # writes config/samples.txt
bash scripts/download_1000g_bams.sh      # downloads BAMs into bamfiles/
bash scripts/run_pipeline.sh             # then run the pipeline as usual
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

### Per-Sample (`results/per_sample/`)

| File                      | Description                       |
|---------------------------|------------------------------------|
| `{sample}_merged_sv.vcf`  | SURVIVOR merged long-read VCF     |
| `{sample}_cxsv_only.vcf`  | CxSV subset of merged VCF         |
| `{sample}_sv_summary.txt` | SV counts per caller + CxSV count |

### Final Call Sets (`results/final/`)

| File                     | Description                             |
|--------------------------|-------------------------------------------|
| `long_read_filtered.vcf` | Population-level merged long-read VCF   |
| `consensus.vcf`          | SVs supported by ≥2 callers             |
| `cxsv_only.vcf`          | CxSV-only VCF                           |
| `cxsv_master_table.tsv`  | 15-column annotated CxSV table          |
| `cxsv_population.bed`    | Population-level CxSV BED for filtering |
| `sv_counts.tsv`          | SV counts by type                       |
| `cxsv_count.txt`         | Total CxSV locus count                  |

### Plots (`results/plots/`)

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

### `glob_wildcards` finds no samples — pipeline exits immediately.

**Cause:** BAM files not in `bamfiles/` or wrong extension. **Fix:**

```bash
ls bamfiles/*.bam     # Should list your BAMs
# Edit config/config.yaml → data_dir if files are elsewhere
```

### Sniffles2 fails with "ERROR: Index file not found".

**Cause:** The `.bam.bai` index was not created by `index_bam`, or the BAM is
not coordinate-sorted. **Fix:**

```bash
samtools sort -o bamfiles/SAMPLE_sorted.bam bamfiles/SAMPLE.bam
samtools index bamfiles/SAMPLE_sorted.bam
```

### CuteSV produces an empty VCF.

**Cause 1:** Working directory (`results/callers/SAMPLE_cutesv_work/`)
already exists with stale temp files from a previous failed run. **Fix:**

```bash
rm -rf results/callers/SAMPLE_cutesv_work/
bash scripts/run_pipeline.sh --resume
```

**Cause 2:** Coverage too low — CuteSV requires ≥5–10× for reliable calls.
Check with `samtools coverage bamfiles/SAMPLE.bam`.

### SVIM outputs `variants.vcf` but needs a very high QUAL threshold.

**Cause:** SVIM scores are not probability-based; the default min QUAL is
appropriate for most datasets but may be too strict for low-coverage data.
**Fix:** Lower `filter_params.svim_min_qual` in `config/config.yaml` (e.g. to
3) for low-coverage samples (<20×).

### SURVIVOR merge fails with "Parsing error" or exits silently.

**Cause 1:** One of the input VCFs contains no records. SURVIVOR requires
≥1 record per input VCF. **Fix:**

```bash
bcftools view -H results/callers/SAMPLE_cutesv_filtered.vcf | wc -l
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
snakemake --touch --cores 1     # Mark all outputs up-to-date
snakemake --cores 32            # Rerun only missing outputs
```

### CxSV count is 0 even though the pipeline completed.

**Cause:** Cluster parameters are too strict for your data. **Diagnosis:**
check `results/final/cxsv_summary.tsv` — are there rows with
`IsCxSV == NO`? **Fix:** relax `cxsv_params.min_breakpoints` in
`config/config.yaml`, then rerun `classify_cxsv` and downstream rules.

---

## Citation

If you use this pipeline, please cite the underlying tools:

- **Sniffles2:** Smolka et al., *Nature Methods* 2024
- **CuteSV:** Jiang et al., *Genome Biology* 2020
- **SVIM:** Heller & Vingron, *Bioinformatics* 2019
- **SURVIVOR:** Jeffares et al., *Nature Communications* 2017
- **bedtools:** Quinlan & Hall, *Bioinformatics* 2010
