# CxSV roadmap

*Regenerated 2026-10-29 after items 1–6 were merged into `main` (cae06a8). Updated 2026-11-03: panel-v1.0
released and tagged (item 12), items 9 and 10 done, offline bundle rebuilt and re-validated with panel-v1.0,
whole-genome positive control against the full panel (branch `validation-panel-v1.0`).*

Goal (unchanged): separate **public** CxSVs, seen in a background panel of
presumed-healthy long-read genomes, from **private** CxSVs in a patient. Then
rank the private ones with a score that flags likely pathogenic events.
Patient data are processed inside a restricted environment, so the pipeline
and everything it needs travel there as an offline bundle; patient reads never
leave it. Workflow engine: Snakemake (decision 2026-10-28).

---

## Decisions waiting on the PI

Each item below has a recommendation; nothing in it starts until the PI says so.

1. **Start the full panel build** (item 12): **decided; done 2026-11-03** (release
   `cxsv-panel-20261102-cd7e1d194cd4`, item 12). Scope (2026-10-29):
   all three cohorts, every verified file, no coverage floor. Dispatch
   (2026-10-30, after the measured estimate): all 5 usable nodes at once
   (freya, shadowfax, athena, betsy, keurig; princess excluded, since it cannot
   run the pinned SIF and is the sequencing host). Projected 37.5 h for the
   per-sample step (36–56 h; the NAS's ~630 MB/s read rate sets the floor),
   then the aggregate step on freya. Build preparation is on branch
   `panel-build-config`: [docs/PANEL_BUILD.md](PANEL_BUILD.md).
2. **Release tag for the first real panel release** (items 2, 6, 12). Tags are
   the PI's to create. *Decided 2026-10-29: `panel-v1.0`. Done 2026-11-03:*
   annotated tag `panel-v1.0` on da6cb90 (the release's source commit);
   `release.json` of `cxsv-panel-20261102-cd7e1d194cd4` records it
   (`source.git_tag`), data files and their manifest rows unchanged.
3. **ClinGen region-level dosage** (item 9). *Decided 2026-10-29: add the
   region list; done, see item 9.* The gene list gives TS 0 for
   *MECP2*, *PMP22* and *SHANK3* because their duplication syndromes are curated
   as regions (e.g. Xq28 *MECP2* duplication, 17p12 CMT1A), so duplications
   there are under-scored today. *Recommendation:* add the region list.
4. **Graded privacy penalty for partly matched loci** (item 10). A locus whose
   junctions are partly found in the panel or a catalogue (`shared_junctions`,
   below the 0.75 fraction) currently gets no privacy penalty; one real host
   locus near *EHMT1* (chr9:150,076,936–150,111,418) reaches Tier 1 this way.
   *Recommendation:* scale the privacy points by the matched fraction, tuned on
   the panel rather than on the positive control. *Decided 2026-10-29; done
   2026-11-03 and merged (PR #17), see item 10.*
5. **Curated literature set of pathogenic CxSVs** (item 11, the open part of
   item 5's "done when"). To decide: which sources, and whether selected
   literature structures should also be re-planted into the positive control
   so detection, not only scoring, is tested. *Decided 2026-11-04: approved,
   both. Done: resource `lit-v1.0`, scoring test (A) and four re-planted
   structures L01-L04 (B); docs/LITERATURE_SET.md, item 11 below.*
6. **ON HOLD.** **When to start item 8** (read-level mechanism axes, PI criteria 3, 4, 6), and
   whether phasing runs on patients only or on the whole panel.
7. **ON HOLD.** **For the restricted environment's IT** (item 2): which workflow engine and
   container runtime are allowed. The bundle assumes Snakemake plus
   Apptainer/Singularity and also ships a Docker archive.
8. **ABANDONED 2026-11-05 (PI).** **SLURM for the lab fleet** (item 7): not pursued;
   fleet dispatch stays ssh-based (`scripts/panel_build/`).
9. **Housekeeping on the NAS:**
   - the sargasso1 `#recycle` bin holds superseded test bundles (~24 GB) and
     earlier deletions, and frees space only when emptied;
   - `/mnt/sargasso1/CxSV/bundles/cxsv-bundle-20261103-dc7f582c01cc.SUPERSEDED`
     (10.5 GB; failed the full-genome positive control at `locus_depth`, see
     docs/RELEASE_AND_BUNDLE.md 6a): **deleted 2026-11-04** (PI-approved);
   - `cxsv-bundle-20261104-65ef970509da` (9.9 GB; bundle of the superseded
     64e10c1b1363 candidate): **deleted 2026-11-05** (PI-approved; its manifests
     are kept in `build/panel-v1.1/provenance/`);
   - the stand-in genome HG01868 (91 GB BAM, public, re-downloadable) sits in
     `datasets/standins/HG01868/` and in freya local scratch;
   - `/mnt/sargasso2/CxSV-old` is scheduled for deletion (nothing to rescue,
     decision 2026-10-28).
10. **Commit identity.** Agent commits on freya use the global git identity
    (Chris Faulk). Keep it, or set a separate author name for agent commits.
11. **Merge `validation-panel-v1.0`** (new, 2026-11-03): fix `33c075a` + docs,
    on top of dc7f582 (`item10-graded-privacy`, merged as PR #17, `b535758`).
    No conflicts: main differs from dc7f582 only by the merge commit. CI green.
    The validated bundle is built from `33c075a`.
12. **Oversized clusters (new, 2026-11-03; proposed item 13).** In the whole
    genome, clustering at 50 kb chains real events to multi-Mb artefact records,
    and clusters over 500 kb are dropped whole. The median panel genome loses
    334 Mb (about 11% of the genome) this way. The positive control's *NF1*
    event E02 is lost, and its span is inside such a cluster in 952/1,619 panel
    genomes. A fix changes the calling fingerprint and needs a new panel release
    (classification + aggregate only). *Decided 2026-11-04: approved. Done
    2026-11-04 (section 13, docs/OVERSIZED_CLUSTERS.md); the panel rebuild is
    decision 14.*
14. **Rebuild the panel with item 13 (new, 2026-11-04).** Re-run only the
    downstream per-sample rules and the aggregate from the panel-v1.0 caller
    outputs into a new build root, tag `panel-v1.1` (proposed). About 3-3.5 h on
    freya with `callable_mask` held (plan and cost: docs/OVERSIZED_CLUSTERS.md
    section 5). *Recommendation:* approve, with the name `panel-v1.1`.
    *Decided 2026-11-04: approved as panel-v1.1. That build
    (`cxsv-panel-20261104-64e10c1b1363`) was held (healthy Tier 1 +86%,
    docs/R10_EXCESS.md). Decided 2026-11-05: adopt `single_caller_bnd_rule:
    svim_lowq` as the default and re-release panel-v1.1 from the merge commit
    into `/mnt/sargasso1/CxSV/build/panel-v1.1-r2`; the 64e10c1b1363 build is
    superseded (build root kept as provenance, never tagged).*
    *Done 2026-11-05: `cxsv-panel-20261105-c31ecc9f8fa2` (docs/panel_v1.1/README.md).*
13. **Small fixes found in the 2026-11-03 validation.** *Approved 2026-11-04;
    done on branch `small-fixes-2026-11` (section "Small fixes 2026-11"
    below), merged (PR #21).* As proposed:
    - `plot_complexity_class.R`: give the heatmap explicit factor levels, so its
      row order no longer depends on the host locale (validation check 7).
    - `run_cxsv.sh`: point `TMPDIR` at the run directory; the container's
      `/tmp` is a 64 MB tmpfs.
    - Record the image id and SIF sha256 in `release.json`.
    - Make the patient workflow refuse or transcode CRAM 3.1 inputs, because
      mosdepth 0.3.10 writes empty masks for them.
    None changes the calling fingerprint. *Recommendation:* approve as one
    small branch, then re-run the bundle validation.
15. **Bundle carrying the panel-v1.1 tag (new, 2026-11-05).** The validated
    bundle `cxsv-bundle-20261105-8fc5bef725bd` embeds the release as it was
    before the tag was recorded (`source.git_tag` null). *Decided 2026-11-05:
    rebuild the bundle from tag `panel-v1.1` with the finalised release, then
    re-run the 8-check validation.*
16. **Alternatives to the long-span exclude rule (2026-11-05).** *Decided
    2026-11-05: the breakends rule and its scoring-side variant are
    abandoned* (24/63 literature events for +57% healthy Tier 1). Paired-breakend
    linkage was not adopted (docs on branch `item13-paired-linkage`, kept
    unmerged for the record). The sensitivity gap goes to item 15 instead.
17. **Library-chimerism QC (item 14, new 2026-11-05).** *Decided 2026-11-05:
    add it.* *Decided 2026-11-05 (after item 14 reported): flag in every patient
    report, conservative mode on by default, no panel exclusion.* Deployed
    2026-11-06 in bundle `cxsv-bundle-20261105-95126e14f036`
    ([library_qc/DEPLOYMENT.md](library_qc/DEPLOYMENT.md)).
18. **Large duplication–triplication events (item 15, new 2026-11-05).**
    Plan: [docs/DUP_TRP_PLAN.md](DUP_TRP_PLAN.md); the PI approved the section 6 defaults.
    Phase 1 report: [docs/CN_LINKAGE_PHASE1.md](CN_LINKAGE_PHASE1.md). Phase 2 waits on the PI.

---

## Status at a glance

| # | Item | Status | Where |
|---|---|---|---|
| 1 | Two workflows, one shared calling module | done 2026-10-28 | `main` |
| 2 | Versioned panel release + offline bundle | done 2026-10-28; real release panel-v1.0 tagged 2026-11-03; bundle re-validated with it 2026-11-03 | `main`, docs/RELEASE_AND_BUNDLE.md |
| 3 | Breakpoint-level matching against the panel | done 2026-10-28 | `main`, docs/JUNCTIONS.md |
| 4 | External SV catalogues | done 2026-10-29 | `main`, docs/EXTERNAL_CATALOGUES.md |
| 5 | Configurable, explainable scoring + classifier fixes | done 2026-10-29; literature set open (item 11) | `main`, docs/SCORING.md |
| 6 | Tests (unit, synthetic CI, lab positive control, bundle validation) | done; bundle validation repeated with panel-v1.0 and a whole-genome positive control added 2026-11-03 | `main`, `validation-panel-v1.0`, docs/POSITIVE_CONTROL.md |
| 7 | SLURM for the lab fleet | **abandoned** 2026-11-05 (PI) | – |
| 8 | Read-level mechanism axes (PI criteria 3, 4, 6) | **on hold** (decision 6) | below |
| 9 | ClinGen region-level dosage | done 2026-10-30; merged | `main`, docs/SCORING.md |
| 10 | Graded privacy penalty for partly matched loci | done 2026-11-03; merged (PR #17) | `main`, docs/SCORING.md section 6b |
| 11 | Curated literature set of pathogenic CxSVs | **done 2026-11-04**: `lit-v1.0` (154 events, 315 junctions on hs1), tests A and B; "top tier" not met on main (7/63 Tier 1), proposals for item 13 | PR #19, `item11-literature-set-2`, docs/LITERATURE_SET.md |
| 12 | Full panel build and first real release | **done 2026-11-03**: release `cxsv-panel-20261102-cd7e1d194cd4`, tag `panel-v1.0`; bundle validated offline (check 7: one locale-dependent plot, data identical); positive control run against it | `main` (tag), `validation-panel-v1.0`, below |
| 13 | Oversized-cluster blind spot | **done 2026-11-04** (fix + panel evaluation); first panel-v1.1 build held (healthy Tier 1 +86%) and superseded; fix `svim_lowq` adopted 2026-11-05; panel-v1.1 released 2026-11-05 as `cxsv-panel-20261105-c31ecc9f8fa2` (healthy Tier 1 0.251, bundle validated, tag pending; docs/panel_v1.1/README.md) | `main` (`de23a0e`), `item13-oversized-clusters-2`, docs/OVERSIZED_CLUSTERS.md |
| – | Small fixes 2026-11 (8 fixes from the panel-v1.0 validation) | done 2026-11-04; merged (PR #21) | `main`, below |
| 14 | Library-chimerism QC (panel and patient) | **done 2026-11-06**: 24/1,619 panel libraries flagged; conservative mode on by default; no panel exclusion; bundle `cxsv-bundle-20261105-95126e14f036` validated 8/8, positive control unchanged | `main` (PR #28, #29), docs/LIBRARY_QC.md, docs/library_qc/DEPLOYMENT.md |
| 15 | Copy-number-anchored detection of large DUP-TRP/INV-DUP and other CN-linked events | **Phase 1 done** (decision 18; Phase 2 awaits the PI) | docs/DUP_TRP_PLAN.md, docs/CN_LINKAGE_PHASE1.md |

---

## Background panel and registry (as of 2026-11-03)

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

**Released panel:** `cxsv-panel-20261102-cd7e1d194cd4` (tag `panel-v1.0`), built from all 1,619 files
(item 12), at `/mnt/sargasso1/CxSV/build/panel-v1.0/results/panel/release` and inside the
validated bundle `/mnt/sargasso1/CxSV/bundles/cxsv-bundle-20261103-33c075a22071`.

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

*2026-11-03:* the bundle validation was repeated with panel-v1.0 on a whole-genome public
stand-in (HG01868), and a whole-genome positive-control driver against a full release was
added (`tests/lab/pos_control_full_genome.sh`, docs/POSITIVE_CONTROL.md).

### 9. ClinGen region-level dosage — done 2026-10-30, merged

`ClinGen_region_curation_list_GRCh38.tsv` (header date 2026-04-28, CC0) was
lifted to hs1 with the catalogue chain `hg38ToHs1.over.chain.gz` (md5
checked). 452 of 513 regions are kept, including 68 of the 69 with HI or TS
1–3; the one lost has no coordinates in the source. The result,
`scoring/clingen_region_dosage.tsv`, is in `resources_manifest.tsv`, so the
bundle ships it and its site config points to it.

New feature `region_dosage` (ClinGen 2A/2B style):
- a lost run that contains an HI region, or a gained run that contains a TS
  region, earns the full points (3: 1.00, 2: 0.45, 1: 0.15);
- a partial overlap earns those points × 0.50 (loss) or × 0.25 (gain) × the
  fraction covered;
- combined with `gene_dosage` as a maximum, never summed.

Done-when check:
- The rule works: unit tests, the synthetic CI (PRIVATE's duplication contains
  a stand-in TS-3 region, +0.10), and a whole-region gain of the real *MECP2*
  or CMT1A region scores 1.00.
- **The positive control's planted E01 (*MECP2*) and E03 (*PMP22*)
  duplications do not meet the rule:** they gain 28.6 of 101.7 kb (28%) and
  9.7 kb of 1.32 Mb (0.7%) of their regions. Partial credit (0.07, 0.00) is
  below the 0.90 their broken HI-3 genes already earn, so they get no region
  points and no total changes. 13/13 events PASS, `known_failures` stays empty,
  and all 28 loci are unchanged (docs/SCORING.md section 6).
- Kleefstra (9q34.3) and Phelan-McDermid (22q13.3) are not curated as regions.

Open for the PI:
- ClinGen 2A's "causative gene fully contained" path (not parseable from the
  list);
- 2F (a loss inside a benign HI-40 region);
- whether to plant a whole-region CMT1A or *MECP2* duplication in the control
  (item 11: not planted; *MECP2* overlaps E01 and no sequence-resolved CMT1A
  complex duplication was obtained; L01-L04 were planted instead).

### 10. Graded privacy penalty for partly matched loci — done 2026-11-03, merged (PR #17)

*Done:* partly matched (`shared_junctions`) loci get the public penalty
scaled by a linear ramp of their matched junction fraction f (0 at f <= 0.55,
full at f >= 0.70; `weights.panel_privacy.graded`, on by default). Tuned by
leave-one-family-out over all 1,619 panel genomes
(`calibration/graded_privacy/`): Tier-1 calls per healthy genome 0.84 -> 0.47,
with at most +1 point chance-demotion risk for private Tier-1 events. On the
10-background regression the EHMT1 host locus drops from Tier 1 (1.00) to
Tier 3 (0.00); no planted event changes tier, score or privacy verdict, on the
regression or scored against the full panel-v1.0. Against the full panel, from the
offline bundle, the whole-genome positive control has 0 Tier-1 host calls
outside the planted events (healthy genomes 0.46 per genome), and the *EHMT1*
host locus is public (14 carriers, Tier 3). Partial catalogue matches
are still not penalised (they cannot be calibrated on the panel without
leakage). Details: docs/SCORING.md section 6b.

Original plan: scale the privacy points by the fraction of a locus's junctions found in the
panel or a catalogue, instead of the all-or-nothing 0.75 cut. Tune on the
full panel (item 12), not on the positive control. **Done when:** the
*EHMT1*-region host locus and similar shared_junctions loci fall below Tier 1,
and no planted event changes tier.

### 12. Full panel build and first real release — done 2026-11-03

*Done.*
- **Build:** launched 2026-10-29 22:47 CDT on freya, shadowfax, athena, betsy and
  keurig, in the pinned SIF (sha256 `9f4e049f…`); calling fingerprint
  `cd7e1d194cd4` unchanged throughout. betsy hung mid-build and was dropped; its
  156 finished samples were kept after a byte-identical 4-sample spot-check.
  The fleet runner fixes were outside the fingerprint. Details:
  [docs/PANEL_BUILD.md](PANEL_BUILD.md).
- **Release** `cxsv-panel-20261102-cd7e1d194cd4`, not a test release:
  1,619 samples, 1,582 individuals (ONT-1019 1,017, LRSC 467, PacBio 135),
  1,304,862 junction sites, 14,428 CxSV loci, a CALLABLE mask per sample, and
  1,634 files (2.43 GB). Annotated tag `panel-v1.0` on da6cb90, recorded in
  `release.json` (`source.git_tag`); manifest sha256 `928d721a…`.
- **Bundle:** `cxsv-bundle-20261103-33c075a22071`, built from `33c075a` with
  this release, validated offline on the public stand-in HG01868, a 1kGP-LRSC
  whole genome that is not in the panel (docs/RELEASE_AND_BUNDLE.md 6a). 7 of
  8 checks pass as scripted. Check 7 flags one plot whose row order follows the
  host locale; all data outputs are identical, and with that plot rendered
  under the image's locale there are 0 unexpected differences.
- **Positive control,** whole genome against this release from the bundle
  (docs/POSITIVE_CONTROL.md):
  - 40/41 junctions exact, 1 partial;
  - E01, E03-E09 private Tier 1; E10 Tier 3; E11/E12 not CxSV;
  - **E02 not called** (oversized cluster; called with item 13: private,
    Tier 1, 1.95 points, docs/OVERSIZED_CLUSTERS.md);
  - E13 private (expected: its carriers exist only in the regression's spliced
    backgrounds);
  - 0 Tier-1 host calls.
- **Found on the way:** the patient `locus_depth` rule failed on a caller's
  off-reference BND mate contig (`chrN`). Fixed in `33c075a` with PI approval;
  outside the fingerprint.

Open at the release, fixed since on `small-fixes-2026-11` (below):
- `release.json` lacks the image id and SIF sha256 (`runtime` is null;
  the bundle records them); panel-v1.0's `release.json` stays as released;
- mosdepth 0.3.10 cannot read CRAM 3.1 and writes an empty mask without
  failing;
- SVIM reads CRAMs without `--reference` (depends on the MD5 cache);
- `plot_cxsv_venn.R` lets sample names collide.

---

### Small fixes 2026-11 — done 2026-11-04, awaiting merge

PI-approved 2026-11-04 (decision 13). Branch `small-fixes-2026-11`, one commit
per fix, each with a test. None changes the calling fingerprint
(`cd7e1d194cd4`): no fingerprinted rule, script or parameter changed.

1. `plot_complexity_class.R`: heatmap rows/columns and panel-A ties in an
   explicit order, independent of the locale (bundle validation check 7).
   Test: `tests/synthetic/check_plot_scripts.sh` (C, POSIX, C.UTF-8,
   en_US.UTF-8 give one PDF).
2. `run_cxsv.sh`: `TMPDIR`/`TMP`/`TEMP` = `/work/.tmp`, also mounted as
   `/tmp`, under Singularity and Docker. Test: `tests/lab/check_container_env.sh`.
3. Release tooling: `release.json` records the SIF sha256 and image id (and
   their source) for future releases; a container build without them is
   refused. panel-v1.0's `release.json` is not rewritten. Tests: unit.
4. Callability vs CRAM 3.1: `callable_mask.sh` transcodes CRAM 3.1 to a
   temporary CRAM 3.0 (or fails, `panel_callability.cram31`), and fails on
   mosdepth decode errors or a 0-bp mask. Test:
   `tests/synthetic/check_cram31_callability.sh` with the fixture
   `tests/fixtures/cram31/`.
5. SVIM and the reference: `REF_PATH`/`REF_CACHE` set in the container
   environment (`run_cxsv.sh`, `lib.sh`); new `check_svim_output` fails a
   sample with SVIM decode errors or contigs cut off. Caller rules and outputs
   unchanged. Tests: unit, `check_container_env.sh`.
6. `plot_cxsv_venn.R`: sample name `<cohort>/<sample>` from the input path,
   so call ids no longer collide across samples. Test: `check_plot_scripts.sh`.
7. Panel-build runner: at start, samples the ledger leaves open get a `done`
   event (finished during a restart) or a `retry` event (re-queued); a live
   instance without a pid file is adopted, not started twice. Tests: unit.
8. Docs: this list; the dc7f582 bundle deletion; item 13 in progress.

## Queued

### 13. Oversized-cluster blind spot (done 2026-11-04; panel rebuild = decision 14)

`bedtools cluster -d 50000` chains records into clusters; a cluster over
`max_cluster_size` (500 kb) is dropped whole (`FILTERED_oversized_cluster`). In
whole genomes, multi-Mb artefact records (for example a 25.7 Mb Sniffles2 DEL
across the chr17 pericentromere) chain real events into such clusters:

- the median panel genome has 12 oversized clusters covering 334 Mb (5th-95th
  percentile 158-636 Mb; 1,616/1,619 genomes have at least one);
- the positive control has 539 Mb, and its *NF1* event E02 is lost there;
- E02's span is inside an oversized cluster in 952 panel genomes;
- the 10-slice regression has only 4.4 Mb, so it could not show this.

Privacy is unaffected (junction-level matching), but CxSV *detection* in these
spans is lost in patients and panel alike. **Proposal:** keep records longer
than `max_cluster_size` out of cluster linkage (or split clusters at them),
check against the positive control and the panel, and rebuild classification
and aggregate. The callers do not need to rerun, but the fingerprint changes,
so this means a new panel release. **Done when:** E02 is called in the whole
genome, the per-genome oversized span drops by an order of magnitude, and the
regression and full-panel positive control pass.

**Result (2026-11-04; docs/OVERSIZED_CLUSTERS.md).** 99.5% of the oversized
span comes from 48,773 records longer than 500 kb: 79% breakend pairs that
SURVIVOR writes as INV intervals, 89% single-caller, 62% recurrent in >= 100
panel genomes, and only about a quarter with a read-depth change. Such records
now take no part in CxSV discovery (`cxsv_params.span_link_max` 500 kb,
`long_span_rule: exclude`); they stay junctions and merged calls and are listed
per sample in `cxsv/long_span_records.tsv`. No locus that panel-v1.0 classified
changes. Blind span per panel genome: median 334 → 1.4 Mb (all four platforms;
what remains are dense SV runs in satellite and acrocentric sequence). Over 100
panel genomes, re-run downstream only (caller outputs byte-identical), healthy
LOFO Tier 1 0.196 → 0.237 per genome (97 genomes; three QC-outlier genomes
carry most of the increase), provisional until the rebuild. Positive control:
E02 private Tier 1; every other planted event and the host's 0 Tier-1 calls
unchanged. The first version merged in `de23a0e` clustered such records at
their breakends instead; that rule raised healthy Tier 1 to 0.371 and gave the
host 3 Tier-1 calls, and is replaced (kept as `long_span_rule: breakends`;
abandoned 2026-11-05 by the PI, decision 16).
Also fixed: `calling_fingerprint` was not re-made after a code change in an
existing build tree.

**panel-v1.1 candidate (2026-11-04; docs/PANEL_BUILD.md, docs/R10_EXCESS.md).**
Rebuilt downstream-only from `main` 65ef970 as `cxsv-panel-20261104-64e10c1b1363`:
VERIFIED, junction sites identical to v1.0, CxSV loci 14,428 -> 17,578, blind span
334 -> 1.40 Mb median per genome. Its bundle passed offline validation (8/8) and the
positive control (E02 Tier 1, host 0 Tier 1). Full-panel healthy LOFO Tier 1 rose
0.457 -> 0.851 per genome (stop rule +50%), so the release is **on hold, not tagged**.
Diagnosis: 88% of the added Tier-1 calls are in 11 genomes already noisy in v1.0,
seven LRSC R10 genomes carrying 3,753-14,258 SVIM-only translocation records at
SVIM's quality floor that cluster with real deletions. Candidate discovery rule
`single_caller_bnd_rule: svim_lowq` (branch `r10-excess-diagnosis`, not merged)
gives 0.251 per genome (R10 0.320), with every positive-control event, the host,
HG01868 and literature test A unchanged; awaiting the PI's decision.
*2026-11-05: adopted as the default (config/calling.yaml); the 64e10c1b1363 build
is superseded and panel-v1.1 is re-released from the merge commit into
`/mnt/sargasso1/CxSV/build/panel-v1.1-r2` (stop rule: healthy Tier 1 must not
exceed v1.0's 0.457 per genome).*

**panel-v1.1 released (2026-11-05; docs/panel_v1.1/README.md).** `cxsv-panel-20261105-c31ecc9f8fa2`,
built downstream-only from `main` 8fc5bef into `/mnt/sargasso1/CxSV/build/panel-v1.1-r2`:
VERIFIED, junction sites identical to v1.0, 13,488 CxSV loci, blind span 1.40 Mb
median per genome. Full-panel healthy LOFO Tier 1 is 0.251 per genome (v1.0 0.457;
R10 0.320). The bundle `cxsv-bundle-20261105-8fc5bef725bd` passed offline validation
(8/8; HG01868 0 Tier 1), and the positive control passed (E02 and E07 Tier 1, host
0 Tier 1). Tag `panel-v1.1` pending.

**Follow-up: residual healthy-genome noise (new, 2026-11-05).** Three genomes keep
>= 10 Tier-1 calls under `svim_lowq` (147 of the panel's 407; without them, 0.160 per genome):
- HG01879 (ONT-1019, 73; v1.0 52) and HG00268 (ONT-1019, 42; v1.0 37): noisy
  before item 13; SVIM-only DEL at QUAL < 10 and single-caller INV, not BND/TRA.
- HG01867 (LRSC R10, 32; v1.0 5), newly noisy: it carries 2,642 SVIM-only TRA records
  at QUAL >= 10, above the `svim_lowq` cut, which v1.0 hid inside 2.29 Gb of oversized clusters.
Candidates: library-chimerism QC (supplementary-alignment rate) or coverage-relative
support; neither is tested (docs/R10_EXCESS.md).


### 8. Read-level mechanism axes (PI criteria 3, 4, 6) — on hold (decision 6)

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

### 11. Curated literature set of pathogenic CxSVs — done 2026-11-04

Details are in docs/LITERATURE_SET.md.

- **Resource:** `lit-v1.0` (`/mnt/sargasso1/CxSV/resources/literature/lit-v1.0`,
  built by `acquisition/literature/`). It holds 154 published germline events:
  63 at junction level from Nazaryan-Petersen 2018, Grochowski 2024,
  Bahrambeigi 2019, Eisfeldt 2024 and Collins 2017, and 91 ClinVar (nstd102)
  region records. Each event carries its source, DOI/PMID, phenotype label and
  the authors' classification. All 315 junctions lifted to hs1. The resource
  is checksummed and shipped by `make_bundle.sh --literature`.
- **Annotation:** `literature_match` (event / junction / region / none), off by
  default; it adds no points.
- **Test A (scoring):**

  | Discovery rule | Events called (of 63) | Tier 1 |
  |---|---|---|
  | v1.0 | 8 | 6 |
  | breakends (de23a0e) | 24 | 18 |
  | main 65ef970 (exclude) | 8 | 7 |
  | interval rule, no size filter | 37 | 32 |

  39 of the 63 events contain junctions longer than 500 kb; exclusion loses
  28 events that another rule calls. 26 events (22 of them DUP-TRP/INV-DUP)
  are not called under any rule, because clustering at 50 kb splits them.
- **Test B (whole-genome positive control, add-ons L01-L04):**
  - 46/46 L-series junctions (48/48 with E02) are called exactly, by bundle
    33c075a and by main 65ef970.
  - L02 (1.04), L03 (1.04) and L04 (1.74) are Tier 1 under both codes.
  - L01 (*PLP1* DUP-TRP/INV-DUP) is not called by either code.
  - On main, L03 is split into 13 fragments because 20 of its 40 junctions
    are excluded.
  - E02 is rescued on main.
- **Done when — not fully met:** the scorer does not rank the literature events
  in the top tier, because of discovery, not weights. The proposals (paired-
  breakend linkage of long records, interleaved inversion pairs, locus
  footprint, biallelic loss) are in docs/LITERATURE_SET.md section 7; the
  first two belong to item 13.

---

### 14. Library-chimerism QC — done 2026-11-06 (decision 17)

Done: metric, per-platform calibration, conservative mode (default on) and the
exclusion study (docs/LIBRARY_QC.md); deployed in bundle
`cxsv-bundle-20261105-95126e14f036` (docs/library_qc/DEPLOYMENT.md, which also lists
three small follow-ups: patient `qc_group`, a report wording fix, release-id
collisions of test releases). The original plan:

In 11 panel genomes (7 LRSC R10, 1 R9, 3 ONT-1019), thousands of SVIM-only
low-QUAL TRA/BND records and an elevated supplementary-alignment fraction point
to chimeric libraries (docs/R10_EXCESS.md). `svim_lowq` removes most of their
Tier-1 noise. Three residual genomes (HG01879, HG00268, HG01867) still hold 36%
of healthy Tier 1. Plan:
1. Define a per-library chimerism metric from the alignments and the call set,
   independent of scoring.
2. Calibrate it per platform on the panel.
3. Report it in every patient run, with a conservative scoring mode for
   flagged libraries.
4. Measure what flagged genomes contribute to the panel before deciding on
   exclusion.

### 15. Copy-number-anchored detection of large DUP-TRP/INV-DUP — Phase 1 done, Phase 2 awaits the PI (decision 18)

None of the 26 published DUP-TRP/INV-DUP events (MECP2, PLP1; 0.42–7.4 Mb) is
called today: their junctions are > 500 kb or > 50 kb apart, and copy number is
not used in discovery. Plan: [docs/DUP_TRP_PLAN.md](DUP_TRP_PLAN.md).
Phase 1 (prototype outside the pipeline, `calibration/cn_linkage/`; 2026-11-06):
[docs/CN_LINKAGE_PHASE1.md](CN_LINKAGE_PHASE1.md). It passes tests A–E: DUP-TRP literature
26/26 called, grid 97.5% at ≥ 20x and ≥ 100 kb, PC L01 Tier 1, healthy Tier 1 +0.01 per
genome on the 100-genome subset, HG01868 0 Tier 1. Recommendation: go for Phase 2.

---

## Abandoned

### 7. A cluster scheduler (SLURM) for the lab fleet — abandoned 2026-11-05 (PI)

Not pursued; fleet dispatch stays ssh-based. Original note: held for a separate discussion about deploying SLURM across the lab nodes.
Most of this round's operational failures came from hand-built ssh dispatch
(claim collisions, stale handles, crash recovery); a scheduler would let
Snakemake place jobs itself and retire `acquisition/orchestrate_alignment.py`.
No work starts on this until that discussion has happened.
