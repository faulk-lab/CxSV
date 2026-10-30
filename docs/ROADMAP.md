# CxSV roadmap

*Regenerated 2026-10-29 after items 1–6 were merged into `main` (cae06a8).*

Goal (unchanged): separate **public** CxSVs, seen in a background panel of
presumed-healthy long-read genomes, from **private** CxSVs in a patient. Then
rank the private ones with a score that flags likely pathogenic events.
Patient data are processed inside a restricted environment, so the pipeline
and everything it needs travel there as an offline bundle; patient reads never
leave it. Workflow engine: Snakemake (decision 2026-10-28).

---

## Decisions waiting on the PI

Each item below has a recommendation; nothing in it starts until the PI says so.

1. **Start the full panel build** (item 12): **decided.** Scope (2026-10-29):
   all three cohorts, every verified file, no coverage floor. Dispatch
   (2026-10-30, after the measured estimate): all 5 usable nodes at once
   (freya, shadowfax, athena, betsy, keurig; princess excluded, since it cannot
   run the pinned SIF and is the sequencing host). Projected 37.5 h for the
   per-sample step (36–56 h; the NAS's ~630 MB/s read rate sets the floor),
   then the aggregate step on freya. Build preparation is on branch
   `panel-build-config`: [docs/PANEL_BUILD.md](PANEL_BUILD.md).
2. **Release tag for the first real panel release** (items 2, 6, 12). Tags are
   the PI's to create. *Recommendation:* a tag such as `panel-v1.0` on the
   commit the release was built from, then re-validate the offline bundle with
   that release.
3. **ClinGen region-level dosage** (item 9). The gene list gives TS 0 for
   *MECP2*, *PMP22* and *SHANK3* because their duplication syndromes are curated
   as regions (e.g. Xq28 *MECP2* duplication, 17p12 CMT1A), so duplications
   there are under-scored today. *Recommendation:* add the region list.
4. **Graded privacy penalty for partly matched loci** (item 10). A locus whose
   junctions are partly found in the panel or a catalogue (`shared_junctions`,
   below the 0.75 fraction) currently gets no privacy penalty; one real host
   locus near *EHMT1* (chr9:150,076,936–150,111,418) reaches Tier 1 this way.
   *Recommendation:* scale the privacy points by the matched fraction, tuned on
   the panel rather than on the positive control.
5. **Curated literature set of pathogenic CxSVs** (item 11, the open part of
   item 5's "done when"). To decide: which sources, and whether selected
   literature structures should also be re-planted into the positive control
   so detection, not only scoring, is tested.
6. **When to start item 8** (read-level mechanism axes, PI criteria 3, 4, 6), and
   whether phasing runs on patients only or on the whole panel.
7. **For the restricted environment's IT** (item 2): which workflow engine and
   container runtime are allowed. The bundle assumes Snakemake plus
   Apptainer/Singularity and also ships a Docker archive.
8. **SLURM for the lab fleet** (item 7): held for a separate discussion.
9. **Housekeeping on the NAS:**
   - the sargasso1 `#recycle` bin holds superseded test bundles (~24 GB) and
     earlier deletions, and frees space only when emptied;
   - `/mnt/sargasso2/CxSV-old` is scheduled for deletion (nothing to rescue,
     decision 2026-10-28).
10. **Commit identity.** Agent commits on freya use the global git identity
    (Chris Faulk). Keep it, or set a separate author name for agent commits.

---

## Status at a glance

| # | Item | Status | Where |
|---|---|---|---|
| 1 | Two workflows, one shared calling module | done 2026-10-28 | `main` |
| 2 | Versioned panel release + offline bundle | done 2026-10-28; tag and real release pending | `main`, docs/RELEASE_AND_BUNDLE.md |
| 3 | Breakpoint-level matching against the panel | done 2026-10-28 | `main`, docs/JUNCTIONS.md |
| 4 | External SV catalogues | done 2026-10-29 | `main`, docs/EXTERNAL_CATALOGUES.md |
| 5 | Configurable, explainable scoring + classifier fixes | done 2026-10-29; literature set open (item 11) | `main`, docs/SCORING.md |
| 6 | Tests (unit, synthetic CI, lab positive control, bundle validation) | done; repeat bundle validation with the first tagged release | `main`, docs/POSITIVE_CONTROL.md |
| 7 | SLURM for the lab fleet | **on hold** | – |
| 8 | Read-level mechanism axes (PI criteria 3, 4, 6) | queued | below |
| 9 | ClinGen region-level dosage | queued, awaiting decision 3 | below |
| 10 | Graded privacy penalty for partly matched loci | queued, awaiting decision 4 | below |
| 11 | Curated literature set of pathogenic CxSVs | queued, awaiting decision 5 | below |
| 12 | Full panel build and first real release | build prepared and estimated (branch `panel-build-config`); dispatch confirmed 2026-10-30; build not yet launched; tag awaits decision 2 | below, docs/PANEL_BUILD.md |

---

## Background panel and registry (as of 2026-10-29)

All data on `/mnt/sargasso1/CxSV/datasets`, aligned to T2T-CHM13v2 (UCSC `hs1`).
`acquisition/build_registry.py` rebuilds `registry/` from the files on disk.

| Cohort | Files | Platform | Notes |
|---|---|---|---|
| `1kGP_ONT_1019_CHM13` (PRJEB89727) | 1,017 CRAMs | ONT | aligned in-house with minimap2 2.24-r1122; read counts match the source `.sra`; 2 ENA-only runs excluded (PI decision) |
| `1kGP_LRSC_ONT_CHM13` | 467 BAMs | ONT R9 + R10 | held 100 plus Tier A1/A2 |
| `1kGP_PacBio_CHM13` | 135 BAMs | PacBio HiFi | Tier A3, all 270 files verified |

1,619 verified files from 1,582 individuals (80.8 TB); 1,065 individuals reach
17x in at least one file, but the panel has no coverage floor (PI, 2026-10-29). Superpopulation labels (from the 1kGP-LRSC metadata and IGSR):
AFR 434, SAS 314, EAS 299, EUR 277, AMR 245, unlabelled 13.

---

## Done

### 1. Two workflows, one shared calling module — done 2026-10-28

`workflow/panel/Snakefile` and `workflow/patient/Snakefile` share
`rules/{common,indexing,calling,per_sample,sample_report,fingerprint}.smk` and
`config/calling.yaml`. CxSVs are discovered **per sample** in both (decision
2026-10-28) and the panel counts distinct carrier individuals. A calling
fingerprint (parameters, reference contig MD5s, tool versions,
loci-defining code) is written by both; the patient workflow stops on a
mismatch. BAM and CRAM give identical results.

### 2. Versioned panel release and offline bundle — done 2026-10-28

- Release id `<name>-<YYYYMMDD>-<fingerprint12>`; `release.json` (commit,
  fingerprint, counts by platform and superpopulation, image id) and
  `release_manifest.tsv` (sha256 of every file). The patient workflow refuses a
  release that does not verify.
- **Callability:** per-sample mosdepth CALLABLE masks (depth ≥ 10, MAPQ ≥ 20).
  Every patient locus reports `Panel_Callable` next to `Panel_Carriers`, so
  "absent" can be told from "could not have been seen".
- One pinned image (`container/Dockerfile` from
  `envs/cxsv-pipeline.explicit.txt`), runnable under Docker, Singularity or
  Apptainer (`container/run_cxsv.sh`).
- `scripts/make_bundle.sh` builds image, test release and bundle from a commit;
  `tests/lab/validate_bundle.sh` validates it with networking verifiably off.
- Validated bundles: `cxsv-bundle-20261028-26085cdd9977` (item 2) and
  `cxsv-bundle-20261029-1dcabd9770f7` (item 5, with scoring resources), both in
  `/mnt/sargasso1/CxSV/bundles/`. A public stand-in patient run offline gave
  records identical to a normal run; tampered releases were refused.
- Not yet: a release tag and a real (non-test) panel release (item 12); a run
  on a separate clean machine (validated on freya with container isolation).

### 3. Breakpoint-level matching — done 2026-10-28

Every sample's merged SV calls become a junction table; the panel clusters
junctions across individuals into sites with carriers and frequencies by
platform and superpopulation. A patient locus is **public** only if a panel
individual carries ≥ 75% of its junctions (500 bp tolerance); statuses are
private / shared_junctions / public, and region overlap is kept as a hotspot
flag. Sniffles2 population mode runs alongside for comparison (it reports no
BND sites and splits insertions by sequence, so matching uses the clustered
sites). Forced genotyping of panel sites in the patient was not implemented:
it cannot change a status. Details: docs/JUNCTIONS.md.

### 4. External SV catalogues — done 2026-10-29

CoLoRSdb v1.2.0, HGSVC3 (HGSVC2024v1.0), HPRC release 2 (MC v2.1) and
gnomAD-SV v4.1 (lifted GRCh38→hs1 per breakpoint, 2.13% site loss), all openly
licensed, checksummed and converted to the shared junction-site schema under
`/mnt/sargasso1/CxSV/datasets/external_catalogues/`. Every patient locus reports
which catalogues contain it and their maximum AF. Catalogue frequencies are
never summed: the catalogues share individuals with each other and with the
panel. Details: docs/EXTERNAL_CATALOGUES.md.

### 5. Configurable, explainable scoring — done 2026-10-29

- ACMG/ClinGen-style points in `config/scoring.yaml`: privacy (panel carriers ×
  callability; catalogue max AF), ClinGen HI/TS and gnomAD LOEUF, exons, class,
  size, caller agreement and read support, a repeat/segdup penalty, the PI's
  chromothripsis criteria 1, 2, 5, 7 and 8, and optional HPO match,
  inheritance and OMIM (off by default). Each call's report lists every
  feature's raw value and points, the total and a tier. Thresholds are
  operational choices, not validated cutoffs.
- Classifier fixes from the positive control (fingerprinted; panels must be
  rebuilt): criteria run on distinct junctions (a simple inversion is no
  longer a CxSV); `Chromothripsis_like` needs > 5 breakends plus interleaving;
  new label `Complex_inversion`.
- `resources/repeats.final.bed` rebuilt from per-element RepeatMasker records
  (genome cover 2.82 → 1.70 Gb).
- Free resources only, versioned in `resources/resources_manifest.tsv`.
  Details: docs/SCORING.md.

### 6. Tests — done

- `tests/unit` (pytest) and `tests/synthetic` (artificial genome with COMMON /
  RARE / PRIVATE events and a hotspot decoy) run in CI on every push.
- Lab positive control (docs/POSITIVE_CONTROL.md): sample 3559 with 13 planted
  heterozygous events (10 CxSVs at *MECP2*, *NF1*, *PMP22*, *SHANK3*, *NRXN1*,
  *AUTS2*, *SATB2*/*FOXP2*, *NSD2*/*EHMT1*/*ZIC2*, *CHD2*, one intergenic; a
  simple DEL and INV; E13 planted in 3 background genomes as a public event).
- Latest regression (2026-10-29, 10 background slices): 13/13 events as
  expected, 40/41 junctions exact, `known_failures.tsv` empty; E01–E09 Tier 1,
  E10 Tier 3, E13 public and Benign_likely; rerun byte-identical.
- Bundle validation recorded for both test bundles; repeat it with the first
  tagged release.

---

## Queued

### 8. Read-level mechanism axes (PI criteria 3, 4, 6)

Deferred by PI decision 2026-10-28. They need read-level evidence that the
per-sample calling path does not keep. Cutoffs are the criteria document's
operational starting points.

| Criterion | Suggested cutoff | What it needs |
|---|---|---|
| **3. Single haplotype of origin** | all break-ends phase to one parental homolog; interspersed LOH where copy number is lost; a balanced event must still pass | read-backed SNV phasing (e.g. Clair3 + WhatsHap in the image), haplotype assignment of each junction's supporting reads; parents for parental origin |
| **4. Reconstructability** (derivative walk) | ≥ 80% of clustered breakpoints spanned by a consistent walk (tunable); gaps only at repeats/centromeres/segdups; gated by the BFB modifier | reads or local-assembly contigs spanning several breakpoints in cis, a junction graph and a walk search |
| **6. Junction microhomology / templated insertion** | from ±50–100 bp of read sequence: 0–6 bp microhomology leans NHEJ/chromothripsis; templated insertions from a distant donor lean MMBIR/FoSTeS | reads or consensus across each junction, alignment of inserted sequence to find donors |

**Done when:** each axis is a computed feature in the score report with a
weight in config (0 until calibrated) and a positive-control check (E06 for
criterion 6; E09 for 3 and 4).

### 9. ClinGen region-level dosage

Add `ClinGen_region_curation_list_GRCh38.tsv` (lifted to hs1) as a scoring
resource, so recurrent-region HI/TS (e.g. *MECP2* duplication, CMT1A, 22q13
deletion) is scored alongside the gene list. **Done when:** duplications of the
*MECP2* and *PMP22* regions receive triplosensitivity points, and the
resource is versioned in the bundle.

### 10. Graded privacy penalty for partly matched loci

Scale the privacy points by the fraction of a locus's junctions found in the
panel or a catalogue, instead of the all-or-nothing 0.75 cut. Tune on the
full panel (item 12), not on the positive control. **Done when:** the
*EHMT1*-region host locus and similar shared_junctions loci fall below Tier 1,
and no planted event changes tier.

### 11. Curated literature set of pathogenic CxSVs

No single public, breakpoint-resolved set of pathogenic germline CxSVs is
known to exist; one has to be assembled (sources and approach in the item's
design note, to be written once decision 5 is made). Candidate sources:
ClinVar via dbVar nstd102 (pathogenic translocations, inversions and complex
records; mostly imprecise coordinates), and published sequence-resolved case
series and reviews of germline complex rearrangements. **Done when:** a
versioned junction-level table of literature CxSVs in hs1 exists, each with
its source and phenotype, and the scorer ranks them in the top tier.

### 12. Full panel build and first real release

*Build preparation done (branch `panel-build-config`,
[docs/PANEL_BUILD.md](PANEL_BUILD.md)).*
- The coverage floor is removed.
- The full DAG resolves for all 1,619 files (25,347 jobs).
- Per-sample cost was measured on 5 full genomes in the pinned SIF: 0.7–1.9
  CPU-h and 20–47 min each.
- A fleet runner handles per-node shards, staging to local disk, ledgers and
  pause/resume.
- The estimate is 37.5 h on 5 nodes. The NAS sets the time: 80.8 TB at about
  630 MB/s is a 35.6 h floor.
- The PI confirmed "all 5 nodes now" on 2026-10-30.

Next: merge, `prepare_build.sh` on the merged commit, launch, then the aggregate
step.

Build the panel over the registry with fingerprint `cd7e1d194cd4`, publish the
first real release, tag it (PI), rebuild and re-validate the bundle with it,
and re-run the positive-control regression against the full panel.
**Done when:** a tagged, checksummed release with callability exists for every
included cohort, and the offline bundle built from that tag passes
validation.

---

## On hold

### 7. A cluster scheduler (SLURM) for the lab fleet

Held for a separate discussion about deploying SLURM across the lab nodes.
Most of this round's operational failures came from hand-built ssh dispatch
(claim collisions, stale handles, crash recovery); a scheduler would let
Snakemake place jobs itself and retire `acquisition/orchestrate_alignment.py`.
No work starts on this until that discussion has happened.
