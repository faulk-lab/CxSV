# CxSV roadmap

Goal (unchanged): separate **public** CxSVs, seen in a background panel of
presumed-healthy long-read genomes, from **private** CxSVs in a patient run.
Then rank the private ones with a score that flags likely pathogenic events.
Patient data must be processed inside a restricted environment, so the
pipeline and everything it needs travel there; patient reads never leave it.

**Decision (2026-10-28):** keep Snakemake as the workflow engine. The queued
items below restructure the pipeline around it; none of them changes the engine.

---

## Done in the current round (October 2026)

Background-panel data, on `/mnt/sargasso1/CxSV/datasets` and all aligned to
T2T-CHM13v2 (UCSC `hs1`):

| Cohort | Files | Platform | Status |
|---|---|---|---|
| `1kGP_ONT_1019_CHM13` (PRJEB89727) | 1,017 CRAMs (4.73 TB) | ONT | aligned in-house with minimap2 2.24-r1122; every CRAM's read count matches its source `.sra` exactly; 2 ENA-only runs excluded by PI decision |
| `1kGP_LRSC_ONT_CHM13` | 467 BAMs | ONT R9 + R10 | held 100 relocated, Tier A1 (247 R10) and A2 (120 R9) downloaded and verified |
| `1kGP_PacBio_CHM13` | 135 BAMs (in progress) | PacBio HiFi | Tier A3 download running; registry to be rebuilt when it finishes |

Code, pushed to feature branches and awaiting review:

- **`cram-input`**: the pipeline accepts BAM or CRAM, either dropped into
  `bamfiles/<cohort>/` or listed in a per-cohort sample sheet with an optional
  coverage floor (`min_coverage`). Adds an MD5-keyed htslib reference cache, so
  CRAMs decode on any node with no network access, and a per-CRAM `@SQ M5`
  check that stops a wrong-assembly CRAM before calling. Includes
  `config/faulk_lab_background.yaml` with the 17x floor. Tested: identical call
  sets from BAM and CRAM, a wrong-M5 CRAM rejected, and full-DAG dry runs of
  all three background sheets.
- **`acquisition-code`**: the `acquisition/` scripts that built the panel
  (manifests, resumable verified downloads, the claim-based multi-node
  alignment dispatcher, integrity records, verified `.sra` cleanup). Also
  `build_registry.py`, which rebuilds `registry/sample_registry.tsv`,
  `individuals.tsv` (one row per person, so an individual sequenced in two
  cohorts is counted once) and the per-cohort sample sheets from the files on
  disk.

Known gap carried forward: no single environment on the fleet has the whole
toolchain. Sniffles2, cuteSV/SVIM and Snakemake sit in three different conda
environments, and SURVIVOR, bcftools and the R plotting packages are missing.
Item 2 removes this.

---

## Queued: structural changes

**Start condition:** `cram-input` and `acquisition-code` are merged into
`main`, and the registry has been rebuilt after the PacBio download.

**Suggested order:** 1 → 6 (test scaffold) → 3 → 2 → 4 → 5. Item 6 grows
with each later item.

### 1. Two workflows, one shared calling module — **done (2026-10-28)**

> Implemented: `workflow/panel/Snakefile` and `workflow/patient/Snakefile`
> share `rules/{common,indexing,calling,per_sample,sample_report,fingerprint}.smk`
> and `config/calling.yaml`. The calling fingerprint (parameters, reference
> contig MD5s, tool versions, loci-defining code) is written by both; the
> patient workflow stops on a mismatch with the panel release. CxSVs are now
> discovered **per sample** in both workflows (decision 2026-10-28) and the
> panel counts distinct carrier individuals per locus. One pinned
> environment: `envs/cxsv-pipeline.yaml` (+ lock file). Verified on freya
> with chr20: BAM and CRAM give identical CxSV tables; the patient's table is
> byte-identical to the same sample's panel table; all its loci are public
> with 2 carriers; a changed `min_support` is rejected; the full 1,006-sample
> background panel resolves to a 13,603-job DAG.


- **Panel builder** (lab fleet): runs over the background cohorts and
  publishes a panel release (item 2).
- **Patient scorer** (restricted environment): one patient, or a trio, per
  run. Calls, genotypes against the panel, annotates and scores.
- Both import the **same** calling and genotyping rules, caller versions and
  parameters. Any difference between the two creates false "private" calls.
- **Done when:** both workflows run from one repo tag; a CI check fails if
  the calling module's parameters differ between them.

### 2. Panel release and an offline bundle for the restricted environment

- The panel builder's product is a small, versioned, checksummed **panel
  release** instead of genomes:
  - CxSV loci with their constituent breakpoints;
  - carrier counts and allele frequencies, split by platform (ONT R9, ONT
    R10, PacBio HiFi) and by superpopulation;
  - a **callability summary**: how many background samples had adequate
    coverage at each locus. Without it, "absent from the panel" cannot be
    told apart from "could not have been seen there".
- The **restricted-environment bundle** holds the container image(s), the
  `hs1` reference plus its MD5 cache, the annotation resources, the panel
  release, the workflow at a tagged commit and a checksum manifest. It needs
  no network access.
- Replace the three-environment setup with one pinned container per tool
  group, or one image, run under Apptainer/Singularity.
- **Done when:** the bundle, built from a tag, runs end to end on a clean
  machine with networking disabled, using a public sample as stand-in patient.
- **Open question for the restricted environment's IT:** which workflow
  engines and container runtime are allowed. Snakemake plus Apptainer is the
  assumption.

### 3. Joint genotyping and breakpoint-level matching

- Today `build_public_reference.py` pools each cohort's CxSV regions into one
  BED. `private_filter.smk` then calls a patient CxSV public if it lies
  within 1 kb (`match_window`). That has two problems:
  - no denominator: a CxSV seen once counts the same as one seen in 800
    people;
  - a different rearrangement at the same hotspot can mask a truly private
    event.
- Use Sniffles2 population mode: merge the per-sample `.snf` files the
  pipeline already writes into cohort genotypes. Match CxSVs on their
  constituent breakpoints and SV types, not on region overlap.
- In the patient, force-genotype the panel's known SVs (Sniffles2
  `--genotype-vcf`), so a common SV the callers missed is not mistaken for
  a private one.
- **Done when:** the panel release carries per-locus frequencies; the
  nested-inversion positive control stays private while a spiked-in common
  panel CxSV is classed public.

### 4. External SV catalogues as additional background

- Bundle published SV frequency resources the same way as the panel. They
  widen the ancestry and platform coverage beyond 1kGP.
- Candidates, each to be checked for licence, reference build (prefer CHM13
  natively) and version: CoLoRSdb, HGSVC3, HPRC release 2 (long-read), and
  gnomAD-SV (short-read; useful for large and balanced events).
- **Done when:** each catalogue is a versioned, checksummed bundle component;
  the scorer reports which catalogue, if any, contains each call.

### 5. Scoring as a separate, configurable module

- Start with an explicit points system modelled on the ACMG/ClinGen CNV
  scoring framework, not machine learning: labelled pathogenic CxSVs are few.
- Features:
  - privacy: frequency in the platform-matched panel plus external
    catalogues, weighted by callability;
  - genes and exons disrupted, weighted by ClinGen haploinsufficiency and
    triplosensitivity scores, OMIM, and HPO match to the patient's phenotype;
  - complexity class and size;
  - read support and agreement between callers;
  - repeat or segmental-duplication context, as an artefact penalty;
  - inheritance (de novo vs inherited) when parents are sequenced.
- Weights live in config, and each call's report lists every feature's
  contribution, so a clinician can see why a CxSV ranked where it did.
- **Done when:** known pathogenic CxSVs (positive control plus a curated
  literature set) rank in the top tier and the score breakdown is in the
  report.

### 6. Tests — **in progress (CI done 2026-10-28; lab positive control pending)**

> Done: `tests/unit` (pytest: sample discovery, fingerprint, panel/patient
> scripts, and a guard that no committed config overrides a fingerprinted
> section and both workflows include the same calling module);
> `tests/synthetic` (seeded, fully artificial 1.2 Mb genome with COMMON /
> RARE / PRIVATE CxSVs planted in 3 panel genomes + a patient: checks
> carriers 3/1/0, public/public/private, private VCF confined to the private
> event, and rejection of a changed fingerprinted parameter);
> `.github/workflows/ci.yml` runs both on every push. Passed on freya in 27 s.
> Pending: the real positive control is being rebuilt with a catalogue of
> planted CxSVs (docs/POSITIVE_CONTROL.md); its lab-only regression test
> follows it.


- CI on every push: a tiny synthetic dataset plus the nested-inversion
  positive control, run through both workflows.
- Before the first patient run, validate the restricted-environment bundle
  with a dry run, then a full run, on public data.
- **Done when:** CI is green on `main`, and the bundle validation is recorded
  with the release tag it used.

---

## On hold

### 7. A cluster scheduler (SLURM) for the lab fleet

Held for a separate discussion about deploying SLURM across the lab nodes.
The motivation: most of this round's operational failures came from
hand-built ssh dispatch (claim collisions, stale handles, crash recovery). A
scheduler would let Snakemake place jobs itself and retire
`acquisition/orchestrate_alignment.py`. No work starts on this until that
discussion has happened.
