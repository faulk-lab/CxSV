# CxSV — Long-Read Complex Structural Variant Detection

Detects complex structural variants (CxSVs) in long-read genomes (Oxford
Nanopore, PacBio HiFi) with a three-caller ensemble (**Sniffles2**, **cuteSV**,
**SVIM**). It then separates **public** CxSVs, seen in a background panel of
presumed-healthy genomes, from **private** ones in a patient: the candidates
for pathogenic events.

The repository holds two Snakemake workflows that share one calling module:

| Workflow | Runs where | Input | Product |
|---|---|---|---|
| **Panel builder** `workflow/panel/Snakefile` | lab fleet | background cohorts (hundreds–thousands of genomes) | a **panel release**: CxSV loci with carrier counts, denominators, and the calling fingerprint |
| **Patient scorer** `workflow/patient/Snakefile` | anywhere, incl. a restricted environment | one patient run (or a family) + a panel release | every patient CxSV labelled public/private with its panel carrier count |

The patient scorer never needs the background reads, only the panel release,
so it can run where patient data must stay.

Planned work is in [`docs/ROADMAP.md`](docs/ROADMAP.md). The code that acquired
and aligned the lab's background panel is in [`acquisition/`](acquisition/README.md).

---

## Why the two sides must be called identically

A patient CxSV is private only relative to a panel built the same way.
Different caller versions, parameters, reference sequence or discovery code
would make some patient calls look private purely because of that difference.
Two mechanisms prevent this:

1. **One parameter file.** Everything that decides which SVs and CxSV loci a
   genome has lives in `config/calling.yaml` (sections marked
   `[FINGERPRINTED]`), read by both workflows.
2. **A calling fingerprint.** Each workflow writes `calling_fingerprint.json`:
   - the fingerprinted parameters;
   - every reference contig's name, length and MD5;
   - the tool versions;
   - a hash of the code that defines CxSV loci (comment and blank lines
     ignored).

   The panel release carries the panel's fingerprint, and the patient workflow
   stops at `check_panel_compatibility` if its own differs, listing every
   difference. (`allow_fingerprint_mismatch: true` overrides this for
   exploratory work only.)

CxSVs are **discovered per sample**: only one genome's own SVs are clustered
and classified, in both workflows. The panel then counts carriers per locus.
(Pooling SVs across many genomes before clustering, as the first version of
this pipeline did, chains SVs from different people into artificial loci; at
about 90 samples the pooled list already covered ~12% of the genome.)

---

## Layout

```
config/
  calling.yaml            shared + [FINGERPRINTED] calling/discovery parameters, reference, annotation
  panel.yaml              panel builder: cohorts, output root
  patient.yaml            patient scorer: runs, panel_release, panel_match
  sites/                  site overlays (paths, cohort sheets), e.g. faulk_lab_panel.yaml
envs/
  cxsv-pipeline.yaml      the one pinned environment for both workflows
  cxsv-pipeline.lock.yaml exact solved versions (conda env export)
workflow/
  panel/Snakefile         panel builder entry point
  patient/Snakefile       patient scorer entry point
  lib/                    plain-Python modules (unit-tested): sample discovery, fingerprint
  rules/
    common.smk            shared: inputs, wildcards, CRAM reference cache env
    indexing.smk          shared: reference cache, CRAM M5 check, BAM/CRAM indexes
    calling.smk           shared, fingerprinted: Sniffles2 / cuteSV / SVIM
    per_sample.smk        shared, fingerprinted: filter, merge, CxSV discovery, annotation
    sample_report.smk     shared: per-sample text report
    fingerprint.smk       shared: tool versions + calling fingerprint
    panel_cohort.smk      panel: per-cohort population call sets and QC plots
    panel_aggregate.smk   panel: carrier table and panel release
    patient_private.smk   patient: compatibility check, public/private labels, plots
  scripts/                scripts invoked by the rules (Python, R)
scripts/
  run_pipeline.sh         `run_pipeline.sh panel|patient [--extra-config FILE] [--resume] [--dry-run]`
acquisition/              background-panel acquisition and alignment (lab fleet)
tests/
  unit/                   pytest: sample discovery, fingerprint, panel/patient scripts, config guard
  synthetic/              seeded artificial genome, panel + patient end to end (runs in CI)
  positive_control/       generator for the lab positive control: design, haplotypes,
                          read simulation and splicing, scoring, structure figure
  lab/                    lab-only regression on the real positive control (never in CI)
docs/
  ROADMAP.md
  POSITIVE_CONTROL.md     the planted-CxSV catalogue, how it is built, latest regression result
```

---

## Quick start

Run everything **from the project root**; relative paths in the configs are
relative to it.

### 1. Environment

```bash
mamba env create -f envs/cxsv-pipeline.yaml
conda activate cxsv-pipeline
```

Caller versions are part of the panel's identity: use the same environment
(or the same container, roadmap item 2) for the panel and every patient.

### 2. Reference and annotation

- `reference`: T2T-CHM13v2 with UCSC `hs1` contig names (`chr1`..`chrM`); index
  it with `samtools faidx`. The RefSeq copy (`NC_*` names) is **not** compatible.
- Annotation BEDs: see [Annotation Resources](#annotation-resources).

Point `reference`, `ref_cache` and `annotation` at them in a site overlay, e.g.
`config/sites/<site>_panel.yaml`.

### 3. Build the panel

Declare background cohorts under `panel_cohorts:` in the site overlay. Each
cohort is either a directory (`bamfiles/<name>/*.bam|*.cram`) or a sample sheet:

```yaml
panel_cohorts:
  my_ont_cohort:               # directory mode: bamfiles/my_ont_cohort/
    platform: ONT
  my_sheet_cohort:
    samples: /data/registry/my_sheet_cohort.tsv
    min_coverage: 17           # keep rows with coverage_x >= 17
    platform: PACBIO           # default when the sheet has no platform column
```

A sample sheet is tab-separated with a header. Required columns are
`sample_id` and `alignment`; optional ones are `index`, `individual_id`,
`platform` and `coverage_x`. `individual_id` is what carriers are counted by,
so one person sequenced in two cohorts counts once. A sheet's indexes must
already exist, because the pipeline never writes into its directories.

```bash
bash scripts/run_pipeline.sh panel --extra-config config/sites/<site>_panel.yaml --dry-run
bash scripts/run_pipeline.sh panel --extra-config config/sites/<site>_panel.yaml
```

### 4. Score a patient

```yaml
# config/sites/<site>_patient.yaml
reference: /restricted/refs/hs1.fa
ref_cache: /restricted/refs/ref_cache
panel_release: /restricted/cxsv_panel/release
patient_runs:
  family_042:
    samples: /restricted/family_042/sheet.tsv
```

```bash
bash scripts/run_pipeline.sh patient --extra-config config/sites/<site>_patient.yaml
```

Do not override any `[FINGERPRINTED]` section in a patient overlay; the run
would stop at the compatibility check.

**CRAM inputs** must have been encoded against `reference`. The pipeline builds
an MD5-keyed htslib cache of it once (`ref_cache`, ~3.1 GB for hs1), points
every tool at it through `REF_PATH` (so no network access is needed), and
checks each CRAM's `@SQ` M5 tags against it before calling. Sniffles2 2.8.1,
cuteSV 2.1.4 and SVIM 2.0.0 give identical call sets from a matched BAM/CRAM
pair.

---

## Outputs

### Per sample — both workflows

`results/<panel|patient>/samples/<cohort>/<sample>/`

| File | Description |
|---|---|
| `callers/{sniffles,cutesv,svim/variants}.vcf` | raw calls; `callers/sniffles.snf` for joint genotyping |
| `merged_sv.vcf` | SURVIVOR merge of the three callers (SUPP_VEC bit order CuteSV, SVIM, Sniffles2) |
| `cxsv/cxsv_summary.tsv` | every SV cluster, `IsCxSV` and the criteria (C1–C5) that fired |
| `cxsv/cxsv_only.vcf` | merged calls inside this sample's CxSV loci |
| `cxsv/cxsv_master_table.tsv` | annotated CxSV loci (genes, exons, repeats, segdups, complexity class) |
| `cxsv/cxsv_loci.bed` | Chr, Start, End, Exon_Disrupt |
| `sv_summary.txt` | per-sample text report |

### Panel — `results/panel/`

| File | Description |
|---|---|
| `release/panel_cxsv_calls.tsv` | one row per (sample, CxSV locus): the table patient loci are matched against |
| `release/panel_cxsv_loci.tsv` | overlapping calls merged across samples, with carrier counts by platform (summary only) |
| `release/panel_summary.tsv` | samples / individuals / calls per cohort, per platform, total |
| `release/calling_fingerprint.json` | what a patient run must reproduce |
| `release/release_manifest.tsv` | sha256 of every release file |
| `panel/plots/cxsv_population_sharing.pdf` | carriers per locus and their genomic distribution |
| `cohorts/<cohort>/` | population VCF, consensus VCF, SV counts, stacked CxSV table, QC plots (SV types, sizes, chromosome density, breakpoint resolution, caller Venn, genomic context, complexity classes) |

### Patient — `results/patient/`

| File | Description |
|---|---|
| `panel_compatibility.txt` | fingerprint comparison with the panel release; written only when compatible (or overridden). The full report, including every difference on a mismatch, is `logs/panel_compatibility.log` |
| `samples/<run>/<sample>/private/cxsv_panel_status.tsv` | every CxSV locus + `Panel_Carriers`, `Panel_Carrier_Freq`, `Panel_Carriers_By_Platform`, `Panel_Status` |
| `samples/<run>/<sample>/private/private_cxsv_master_table.tsv` | private loci only |
| `samples/<run>/<sample>/private/private_cxsv.vcf` | the patient's SV records inside private loci |
| `samples/<run>/<sample>/plots/` | genomic context and complexity class plots |

A patient locus is **public** when at least `panel_match.min_carriers`
(default 1) panel individuals have a CxSV overlapping it after padding both
sides by `panel_match.window` bp (default 1000). The carrier count is kept on
every row, so a rare-but-seen locus can be told apart from one absent from the
panel. Breakpoint-level matching and callability weighting replace interval
overlap in roadmap item 3.

---

## Annotation Resources

`resources/` is gitignored — several of these files are hundreds of MB and
exceed GitHub's size limits, so they're fetched/regenerated locally instead
of committed. Paths are set in `config/calling.yaml → annotation:` (override per site).

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

## Tests

| Test | Data | Where it runs | Command |
|---|---|---|---|
| Unit | none | CI + anywhere | `pytest tests/unit` |
| Synthetic end to end | seeded artificial 1.2 Mb genome, 3 panel genomes + 1 patient | CI on every push | `bash tests/synthetic/run_synthetic.sh` |
| Positive control | sample 3559 with 12 planted heterozygous CxSVs (lab NAS only) | lab fleet | `bash tests/lab/pos_control_regression.sh [WORKDIR] [N_BACKGROUND]` |

The positive control is described in [docs/POSITIVE_CONTROL.md](docs/POSITIVE_CONTROL.md):
ten complex rearrangements at clinically relevant genes and two simple-SV
specificity controls, each with a breakend-level truth set. Its reads are
human-subject data and never leave the lab NAS; the repository holds only the
generator, the event catalogue and the truth coordinates.

![Planted positive-control events](docs/figures/pos_control_events.png)

---

## Troubleshooting

### `[panel] no cohorts/runs configured` / `[patient] no cohorts/runs configured`
Nothing under `panel_cohorts:` / `patient_runs:`. Both are empty in the
committed configs on purpose; declare them in a site overlay passed with
`--extra-config`.

### `sample sheet ... has N problem(s)`
A row points at a missing alignment or index, repeats a `sample_id`, or names a
file that is not `.bam`/`.cram`. The error lists the lines. Fix or regenerate
the sheet; index sheet files yourself (`samtools index FILE`).

### `check_panel_compatibility`: `this run is not comparable with the panel`
The patient run's calling fingerprint differs from the panel release's.
`results/patient/logs/panel_compatibility.log` lists each difference (parameter,
contig, tool version or code file). Use the panel's environment and do not
override `[FINGERPRINTED]` parameters; if the change is intended, build a new
panel release with it.

### `check_cram_reference` fails: `M5 ... is not in the pipeline reference`
The CRAM was encoded against a different reference (e.g. GRCh38, or CHM13 with
RefSeq `NC_*` names). Point `reference` at the FASTA it was made with, or
re-encode: `samtools view -C -T hs1.fa -o out.cram in.cram`.

### A caller fails on a CRAM: `Failed to populate reference`
The MD5 cache is missing or incomplete. Delete `<ref_cache>/md5_manifest.tsv`
and re-run; `build_ref_cache` rebuilds it.

### cuteSV produces an empty VCF
A stale `callers/cutesv_work/` from a failed run, or coverage too low (cuteSV
needs roughly 5–10x). Remove the work directory and resume with `--resume`.

### SURVIVOR merge fails or exits silently
SURVIVOR needs at least one record per input VCF. Check
`bcftools view -H <sample>/callers/<caller>_filtered.vcf | wc -l` and the caller log.

### Pipeline reruns rules that already completed
Timestamps changed (e.g. after copying files). `snakemake -s workflow/<panel|patient>/Snakefile --touch ...`
marks outputs up to date.

---

## Citation

If you use this pipeline, please cite the underlying tools:

- **Sniffles2:** Smolka et al., *Nature Methods* 2024
- **CuteSV:** Jiang et al., *Genome Biology* 2020
- **SVIM:** Heller & Vingron, *Bioinformatics* 2019
- **SURVIVOR:** Jeffares et al., *Nature Communications* 2017
- **bedtools:** Quinlan & Hall, *Bioinformatics* 2010
