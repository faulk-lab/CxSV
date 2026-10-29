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
Item 2 removes this. Since 2026-10-28 there is one env (`cxsv-pipeline`) and
one image built from its explicit package list (`container/Dockerfile`).

---

## Queued: structural changes

**Start condition:** `cram-input` and `acquisition-code` are merged into
`main`, and the registry has been rebuilt after the PacBio download.

**Suggested order:** 1 → 6 (test scaffold) → 3 → 2 → 4 → 5. Done: 1, 6, 3, 2, 4, 5 (item 2 and 5 branches
awaiting merge). Item 6 grows with each later item. Queued next: 8 (read-level mechanism axes).

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

### 2. Panel release and an offline bundle for the restricted environment — **done (2026-10-28; release tag pending the PI)**

> Done (branch `item2-panel-release`; docs/RELEASE_AND_BUNDLE.md):
>
> - **Versioned release.** Id `<name>-<YYYYMMDD>-<fingerprint12>`.
>   `release.json` carries the source commit, fingerprint, samples and
>   individuals by platform and superpopulation, tools and the container image
>   id / SIF sha256. `release_manifest.tsv` checksums every file.
> - **Callability.** Per-sample mosdepth CALLABLE masks ship in the release.
>   Each junction site gets the number of panel samples / individuals with
>   depth ≥ 10 (MAPQ ≥ 20) at both breakpoints. Every patient locus reports
>   `Panel_Callable` next to `Panel_Carriers`.
> - **Refusal.** The patient workflow refuses a release that does not verify
>   (checked at parse time, in a recorded rule, and inside the classifier).
> - **One pinned image.** `container/Dockerfile` installs
>   `envs/cxsv-pipeline.explicit.txt` (URL + MD5 per package). The lock had
>   drifted from the env and was regenerated. mosdepth 0.3.10 was added.
> - **Builder.** `scripts/make_bundle.sh image|bundle --commit` builds the
>   image, the SIF and the bundle from a `git archive` of the commit.
> - **Validated offline** (`tests/lab/validate_bundle.sh`):
>   - bundle `cxsv-bundle-20261028-26085cdd9977`: commit 26085cd, 7.5 GiB,
>     about 9 min to build;
>   - test release `TEST-posctrl10slice-20261028-ad454c604bdd`, built from the
>     10 regression slices in the bundle's SIF;
>   - the public stand-in HG00130 (not in the panel) ran end to end under
>     Singularity and Docker with networking verifiably off (loopback only,
>     TCP/DNS/HTTPS fail, positive control passes);
>   - records identical to a normal conda run;
>   - tampered releases were refused before any job ran.
> - **Positive-control regression** at 26085cd unchanged: 40/41 junctions
>   exact, E12 the only (known) failure.
>
> Not yet: a release **tag** (the PI's decision; the commit hash is
> recorded), a run on a separate clean machine (validated on freya with the
> container isolating the run), and a real (non-test) panel release.

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
  Met for a commit (26085cd) on freya with container isolation; the tag is
  the PI's to create.
- **Open question for the restricted environment's IT:** which workflow
  engines and container runtime are allowed. Snakemake plus Apptainer is the
  assumption. The bundle also ships a Docker archive, and `run_cxsv.sh` runs
  under Singularity, Apptainer or Docker.

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

### 4. External SV catalogues as additional background — **done (2026-10-29)**

> Acquisition (2026-10-28): CoLoRSdb v1.2.0, HGSVC3 HGSVC2024v1.0, HPRC
> release 2 (MC v2.1) and gnomAD-SV v4.1 were downloaded with verified
> checksums to `/mnt/sargasso1/CxSV/datasets/external_catalogues/`. gnomAD-SV
> was lifted to hs1 breakpoint by breakpoint (2.13% site loss). Each catalogue
> was converted to the shared junction-site schema with a manifest
> ([`acquisition/catalogues/`](../acquisition/catalogues/README.md); survey,
> licences and caveats in [`docs/EXTERNAL_CATALOGUES.md`](EXTERNAL_CATALOGUES.md)).
>
> **Integration (2026-10-29, branch `item5-scoring`):**
>
> - Rule `catalogue_hits` (`workflow/scripts/match_catalogues.py`) matches
>   every junction of every patient CxSV locus against each catalogue with
>   the item-3 matcher (`cxsv_matching.SiteIndex.query`, 500 bp, one
>   implementation; tolerance and fraction configurable in
>   `config/scoring.yaml → external_catalogues`).
> - Per locus and catalogue it reports junctions matched / total, max AF, max
>   group AF and whether the catalogue *contains* the locus (≥ 0.75 of its
>   junctions), plus `Catalogues_Containing` and `Catalogue_Max_AF`. A
>   per-junction side table gives the matched sites.
> - Each table is verified against its manifest's sha256 before use
>   (`catalogue_check.tsv`); a mismatch stops the run.
> - The bundle ships the four tables with their manifests (as since item 2),
>   and `bundle_site.yaml` points the scorer at them.
> - Catalogue frequency is a scoring feature (most severe catalogue; never
>   summed).
>
> **Result** on the positive control: no planted event is in any catalogue.
> Of the host's 16 other loci, 11 are contained in ≥ 1 catalogue (max AF
> 0.45–0.83).

- Bundle published SV frequency resources the same way as the panel. They
  widen the ancestry and platform coverage beyond 1kGP.
- Candidates, each to be checked for licence, reference build (prefer CHM13
  natively) and version: CoLoRSdb, HGSVC3, HPRC release 2 (long-read), and
  gnomAD-SV (short-read; useful for large and balanced events).
- **Done when:** each catalogue is a versioned, checksummed bundle component;
  the scorer reports which catalogue, if any, contains each call.

### 5. Scoring as a separate, configurable module — **done (2026-10-29)**

> Done (branch `item5-scoring`; [docs/SCORING.md](SCORING.md)):
>
> - **Scorer.** `workflow/lib/cxsv_scoring.py` + `workflow/scripts/score_cxsv.py`
>   + `rules/patient_scoring.smk`, weights and tiers in `config/scoring.yaml`.
>   It is an ACMG/ClinGen-style points system:
>   - privacy (panel carriers × callability; catalogues, max AF, never summed);
>   - genes broken / lost / gained, weighted by ClinGen HI/TS and gnomAD LOEUF;
>     exons; gene count; no gene;
>   - class and size; caller agreement and read support;
>   - a repeat / segdup artefact penalty;
>   - the PI's chromothripsis criteria 1 (interleaving), 2 (breakpoint tier),
>     5 (read-depth copy-number oscillation, sex-aware, PAR excluded),
>     7 (clustering vs a chromosome-intrinsic null, descriptive) and
>     8 (fold-back / BFB modifier);
>   - optional HPO match, inheritance (parents in the run) and OMIM (licence
>     key; off by default).
>
>   Each call gets a TSV row, long-format feature rows and a markdown report
>   with every feature's raw value, rule and points, the total and a tier.
>   Thresholds are operational choices, not validated cutoffs. Criteria 3, 4
>   and 6 are queued as item 8.
> - **Resources** (free only), versioned with URL, date and sha256 in
>   `resources/resources_manifest.tsv`
>   (`acquisition/resources/build_resources.py`):
>   - ClinGen gene dosage list (2026-10-28);
>   - gnomAD v4.1 constraint;
>   - HPO v2026-09-01;
>   - CHM13v2 PARs;
>   - CHM13v2 RepeatMasker 4.1.2p1.2022Apr14.
>
>   All are shipped in the bundle.
> - **Classifier fixes** found by the positive control (fingerprinted; panels
>   must be rebuilt):
>   - C1–C5 run on distinct junctions (duplicate / reciprocal records
>     collapsed, SVIM composites dropped), so a simple inversion is no longer
>     a CxSV (E12);
>   - `Chromothripsis_like` needs the extreme tier (> 5 breakends) plus
>     interleaving, and C5 no longer overrides type labels;
>   - new label `Complex_inversion`;
>   - the README class table lists every label.
> - **Repeats rebuilt.** `repeats.final.bed` now holds per-element
>   RepeatMasker records. Genome cover dropped from 2.82 Gb to 1.70 Gb. On
>   *MECP2*, repeat cover went from 83.8% to 50.6% (soft mask 50.4%), and exon
>   cover from 31.0% to 3.6%. The old file is kept as `.bak`.
> - **Result** (lab regression, fresh panel, `RERUN=1`):
>   - E01–E09 are all Tier 1 (1.37–2.74); E10 scores 0.33 (Tier 3) and the
>     public E13 −1.15 (Benign_likely);
>   - E12 is no longer a CxSV, and `known_failures.tsv` is empty;
>   - 40/41 junctions exact;
>   - the rerun was byte-identical.
>   - Class concordance went from 1/11 to 11/11 against the revised
>     expectations (4/11 against the original labels; 7 expectations revised
>     with reasons in `events.yaml`).
> - Not yet: the curated literature set of pathogenic CxSVs named in "done
>   when" (only the planted control so far); ClinGen's region-level dosage
>   list; calibration of the weights on real patients.

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

### 6. Tests — **done (2026-10-28; bundle validated 2026-10-28, release tag pending)**

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
  with a dry run, then a full run, on public data. Done for the test bundle
  `cxsv-bundle-20261028-26085cdd9977` (commit 26085cd; docs/RELEASE_AND_BUNDLE.md
  section 6); repeat with the first real panel release and its tag.
- **Done when:** CI is green on `main`, and the bundle validation is recorded
  with the release tag it used.

### 8. Read-level mechanism axes (PI criteria 3, 4, 6) — **queued (deferred by PI decision 2026-10-28)**

The PI's germline-chromothripsis criteria (germline-chromothripsis-pipeline-criteria)
have eight axes. Criteria 1, 2, 5, 7 and 8 are computed and scored since item 5
(docs/SCORING.md). The three axes below need read-level evidence (phased
SNVs, read sequence, reads spanning several breakpoints) that the per-sample
calling path does not keep. They are not implemented. As with every number
in this pipeline, the cutoffs are the doc's operational starting points, not
validated thresholds.

| Criterion | Doc's suggested cutoff | What it needs |
|---|---|---|
| **3. Single haplotype of origin** (interspersed heterozygosity) | all break-ends phase to one parental homolog; where CN loss exists, retained-heterozygosity segments interspersed with LOH; phrased so a balanced, copy-neutral event still passes (haplotype of the break-ends themselves) | read-backed SNV phasing of the patient (WhatsHap / margin / native ONT phasing; a small-variant caller such as Clair3 in the image), then assignment of each junction's supporting reads to a haplotype block; parents (trio) for parental origin |
| **4. Reconstructability** (derivative-chromosome walk) | graded score W = % of clustered breakpoints spanned by a consistent walk; ≥ 80% (tunable) via direct long-read span, local assembly contig or graph-inferred adjacency; gaps allowed only at annotated repeat / centromere / segdup; walk score gated by the BFB modifier (criterion 8) | reads or local-assembly contigs that span several breakpoints in cis (supporting-read names per junction from the callers, or a local assembly of each locus), a junction graph and a walk search; repeat/segdup tracks (present) for excused gaps |
| **6. Junction templating / microhomology** (primary discriminator) | per junction, from ±50–100 bp flanking sequence taken from the reads: short microhomology (~0–6 bp) with blunt joins leans NHEJ / chromothripsis; longer templated insertions homologous to a distant donor lean MMBIR / FoSTeS; leaning, not deterministic | the reads (or a consensus) across each junction, so the junction sequence and any inserted sequence are known; alignment of insertions to find donor loci; the reference for microhomology at both ends |

- **Done when:** each axis is a computed feature with its raw value in the
  score report, its weight in config/scoring.yaml (0 until calibrated), and a
  positive-control check (E06 templated insertions for criterion 6; E09 for 3
  and 4).

---

## On hold

### 7. A cluster scheduler (SLURM) for the lab fleet

Held for a separate discussion about deploying SLURM across the lab nodes.
The motivation: most of this round's operational failures came from
hand-built ssh dispatch (claim collisions, stale handles, crash recovery). A
scheduler would let Snakemake place jobs itself and retire
`acquisition/orchestrate_alignment.py`. No work starts on this until that
discussion has happened.
