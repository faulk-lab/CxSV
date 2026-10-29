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

**Suggested order:** 1 → 6 (test scaffold) → 3 → 2 → 4 → 5. Done: 1, 6, 3. Item 6 grows
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

### 3. Joint genotyping and breakpoint-level matching — **done (2026-10-28)**

> **Implemented** on branch `item3-breakpoint-matching`; details in
> [docs/JUNCTIONS.md](JUNCTIONS.md).
>
> - **Junction tables.** Every sample's merged SV calls become a junction
>   table (`rules/junctions.smk`, fingerprinted, with `junction_params`).
>   Duplicate and reciprocal breakend records are collapsed for matching only;
>   the C1-C5 criteria are unchanged.
> - **Panel junction sites.** The panel clusters all junctions across
>   individuals into sites in the shared junction-site schema. Clustering is
>   type-aware and leader-based, with 500 bp tolerance. Each site carries
>   distinct-individual carriers and allele frequencies by platform and
>   superpopulation (`panel_junction_sites.tsv`, `panel_junction_carriers.tsv`,
>   `panel_individuals.tsv`).
> - **Sniffles2 population mode.** The per-sample `.snf` files are jointly
>   genotyped (`sniffles_population_sites.tsv`) and compared with the clustered
>   sites.
> - **Patient matching.** Patient loci are matched junction by junction: a
>   panel individual carries a locus when it has >= 75% of its junctions.
>   Statuses are private / shared_junctions / public, and region overlap is
>   kept as a hotspot flag.
> - **Tests.** Unit tests (29) cover extraction, collapsing, clustering and the
>   matcher. The synthetic CI test gained a HOTSPOT DECOY (a different CxSV at
>   the PRIVATE locus in 2 panel genomes: PRIVATE stays private with 2
>   region-overlap carriers; COMMON stays public with 3) and a breakpoint check
>   of the SIMPLE deletion (3 carriers).
> - **E13.** The positive control gained E13, a CxSV spliced at regression
>   time into the control and 3 background genomes.
>
> **Result** (lab regression, 10 backgrounds, `RERUN=1`):
>
> - 40/41 junctions exact.
> - **E04 now PASSES as private**: 0 carriers, 4 region-overlap carriers
>   flagged as a hotspot.
> - **E13 is public with 3 carriers.**
> - E11/E12's junctions are matched on breakpoints with 0 carriers.
> - E01-E03 and E05-E11 still pass. E12 stays a known failure (classifier,
>   item 5).
> - The rerun gave byte-identical records.
>
> **Sniffles2 vs clustered frequencies** on that panel:
>
> - 954 of 964 Sniffles2 sites match a clustered site.
> - Summed per clustered site, carrier calls agree: DEL 1,241 vs 1,203; INS
>   1,473 vs 1,458.
> - Per site, Sniffles2 splits 129 clustered INS sites into several
>   sequence-distinct alleles, so per-Sniffles2-site INS carriers are lower.
> - Carrier counts are equal at 236/310 DEL sites and 209/257 one-to-one INS
>   sites.
> - Sniffles2 population mode reports no BND, so the 7 clustered BND sites
>   have no Sniffles2 counterpart.
>
> **Not implemented: forced genotyping of panel sites in the patient**
> (Sniffles2 `--genotype-vcf`). It cannot change a status, because it only
> adds patient junctions that are public by construction. The false-private
> error comes from panel callers missing the patient's junctions, which needs
> panel evidence rather than patient re-genotyping. It would also cost one more
> genome-wide Sniffles2 pass per patient. Reasoning in docs/JUNCTIONS.md.

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
- **Done when:** the panel release carries per-locus frequencies; every
  planted event in the positive control (docs/POSITIVE_CONTROL.md) stays
  private while a spiked-in common panel CxSV is classed public, and E11/E12
  (simple SVs) are matched on breakpoints rather than region overlap.

### 4. External SV catalogues as additional background — **acquisition done (2026-10-28); matching integration pending**

> Done: CoLoRSdb v1.2.0, HGSVC3 HGSVC2024v1.0, HPRC release 2 (MC v2.1) and
> gnomAD-SV v4.1 are downloaded with verified checksums to
> `/mnt/sargasso1/CxSV/datasets/external_catalogues/`, gnomAD-SV is lifted to hs1
> breakpoint by breakpoint (2.13% site loss), and each is converted to the shared
> junction-site schema with a manifest ([`acquisition/catalogues/`](../acquisition/catalogues/README.md);
> survey, licences and caveats in [`docs/EXTERNAL_CATALOGUES.md`](EXTERNAL_CATALOGUES.md)).
> None had to be excluded. Still to do: bundle them as panel-release components
> and make the scorer report the catalogue hits per call.

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
- Known issue to fix first: `resources/repeats.final.bed` holds UCSC
  RepeatMasker "joined" spans, which overstate repeat cover (e.g. *MECP2* is
  annotated LINE although only ~1.9% of it is soft-masked in hs1). Rebuild it
  from per-element RepeatMasker records (or the hs1 soft mask) before any
  repeat-based penalty is weighted.
- **Done when:** known pathogenic CxSVs (the planted gene events E01-E09 of the
  positive control, plus a curated literature set) rank in the top tier, the
  intergenic and simple-SV controls (E10-E12) rank below them, and the score
  breakdown is in the report.

### 6. Tests — **done (2026-10-28; bundle validation pending with item 2)**

> Done: `tests/unit` (pytest: sample discovery, fingerprint, panel/patient
> scripts, and a guard that no committed config overrides a fingerprinted
> section and both workflows include the same calling module);
> `tests/synthetic` (seeded, fully artificial 1.2 Mb genome with COMMON /
> RARE / PRIVATE CxSVs planted in 3 panel genomes + a patient: checks
> carriers 3/1/0, public/public/private, private VCF confined to the private
> event, and rejection of a changed fingerprinted parameter);
> `.github/workflows/ci.yml` runs both on every push. Passed on freya in 27 s.
> Lab positive control rebuilt (docs/POSITIVE_CONTROL.md): sample 3559 with
> 12 heterozygous planted events (10 CxSVs at MECP2, NF1, PMP22, SHANK3,
> NRXN1, AUTS2, SATB2/FOXP2, NSD2/EHMT1/ZIC2, CHD2, one intergenic; a simple
> DEL and INV as specificity controls) and a 37-junction truth set;
> `tests/lab/pos_control_regression.sh` scores it. First result: 36/37
> junctions exact, 10/12 events as expected; E04 (hotspot overlap -> item 3)
> and E12 (simple inversion called CxSV -> classifier fix) are listed known
> failures.


- CI on every push: the tiny synthetic dataset through both workflows. The
  real positive control is too large for CI and runs on the lab fleet
  (`tests/lab/pos_control_regression.sh`).
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
