# Copy-number track: copy-number-anchored CxSV loci (roadmap item 15)

The junction track (docs/JUNCTIONS.md, `classify_cxsv.py`) builds a CxSV locus from
SV junctions that lie close together. It misses large duplication–triplication
events (DUP-TRP/INV-DUP; Carvalho et al. 2011) whose junctions are hundreds of kb
apart or partly undetected: in the literature set none of 26 DUP-TRP events was
called on `main` before this track (docs/CN_LINKAGE_PHASE1.md, test A).

The copy-number track is a second discovery track in the shared calling module. It
**adds loci and removes none**: every junction-track locus is kept unchanged.

- Plan and decisions: [DUP_TRP_PLAN.md](DUP_TRP_PLAN.md) (section 6 defaults approved by the PI).
- Phase 1 evaluation (tests A–E, prototype in `calibration/cn_linkage/`):
  [CN_LINKAGE_PHASE1.md](CN_LINKAGE_PHASE1.md).
- This document: the pipeline implementation (Phase 2).

## 1. Steps (workflow/rules/cn.smk)

| Rule | Script | Output (per sample unless noted) |
|---|---|---|
| `cn_params_json` | – | `results/cn_params.json`: `cn_params` + the `cxsv_params` / `junction_params` the linkage re-uses |
| `cn_depth` | `scripts/cn_depth.sh` | `cn/depth.regions.bed.gz` (mosdepth, 1 kb bins, MAPQ ≥ 20, fast mode), `cn/depth.sha256` |
| `panel_cn_pon` (panel) | `scripts/cn_pon.py` | `panel/cn/pon/<stratum>.npz`, `panel/cn/pon/pon_members.tsv` |
| `cn_segment` | `scripts/cn_segment.py` | `cn/segments.tsv`, `cn/stats.json` |
| `panel_cn_segments` (panel) | – | `panel/cn/panel_cn_segments.tsv` (every panel library's segments) |
| `sample_cn_linkage` | `scripts/cn_linkage.py` | `cxsv/clusters.bed`, `cxsv/cxsv_summary.tsv`, `cxsv/cxsv_only.vcf` (junction-track + copy-number loci); `cn/cn_candidates.tsv`, `cn/cn_loci.tsv`, `cn/cn_only.tsv` |

The junction-track rules now write `cxsv/junction_clusters.bed`,
`cxsv/junction_summary.tsv` and `cxsv/junction_only.vcf`. `sample_cn_linkage`
copies them and appends the copy-number loci, so every rule downstream
(`cxsv_master_table`, junction tables, panel aggregation, patient scoring) reads
the same file names as before.

Library code: `workflow/lib/cxsv_cn.py` (binning, sex, ploidy, PoN, HMM
segmentation, nested scan, refinement) and `workflow/lib/cxsv_cnlink.py`
(candidates, record attachment, classes, public filter). Both are ports of the
Phase 1 prototype (`calibration/cn_linkage/cnlib.py`, `cnlink.py`); the parameter
values are the Phase 1 values (`tests/unit/test_cn_track.py` checks that
`config/calling.yaml` equals the library defaults).

### 1a. Segmentation

1. Depth in 1 kb bins, averaged to 10 kb; genome median over autosomes.
2. Sex from the chrX / autosome median ratio (`sex_ratio_male` 0.75). Ploidy track:
   autosomes 2; chrX 2 (female) or 1 (male, PAR 2); chrY 1 in males (PAR 0).
   A genome without a chrX track is treated as diploid.
3. Ratio to the stratum PoN median. Bins are masked where the PoN median is below
   0.25, its relative MAD above 0.30 or fewer than 4 PoN genomes have data.
4. Viterbi HMM over integer copy states 0–6; segments with at least 3 bins and a
   copy change ≥ 0.6; boundaries refined on 1 kb bins (± 15 kb).
5. Nested scan inside long gain segments (≥ 60 kb): a sub-segment ≥ 30 kb at
   z ≥ 6 and ≥ 0.4 copies above its host becomes a nested TRP (PC L01: 120 kb TRP
   inside a DUP).

### 1b. Linkage and classes

Adjacent abnormal segments (gap ≤ 500 kb) form a candidate. SV records of the
sample attach at segment boundaries (± 10 kb, widened over masked bins; an INV's
other end may lie up to 100 kb outside), with no distance limit between records.
Classes: `DUP-TRP/INV-DUP`, `CN_nested_gain`, `CN_linked_gain`, `CN_linked_loss`,
`CN_linked_gain_loss`. Decisions:

| Decision | Meaning | Result |
|---|---|---|
| `locus` | linked change with ≥ 1 attached junction | written as a CxSV locus (`CxSV_Criteria` starts with `CN_linked`) |
| `simple_cnv` | one segment explained by one DEL/DUP record | nothing (a simple CNV) |
| `public_cn` | every segment also seen (≥ 50 % reciprocal) in ≥ 1 other panel individual | nothing (polymorphic) |
| `cn_only` | linked change with no junction | listed in `cn/cn_only.tsv`; the patient report lists them untiered |

A copy-number locus whose id equals an existing junction-track locus is not
written again (`cn/cn_loci.tsv`, `written`).

## 2. Panel-of-normals and strata

One PoN per stratum; the stratum names are the library-QC groups of item 14:
`ONT_1kGP_1019`, `ONT_LRSC_R9`, `ONT_LRSC_R10`, `PacBio_HiFi`.

- **Members:** every panel library of the stratum except the library-QC flag list
  (`cn_params.pon_exclude` = `calibration/library_qc/flagged_v1.1.tsv`, 24
  libraries). `pon_members.tsv` lists every library with sex, chrX ratio, genome
  median and `in_pon`.
- **chrX** pools both sexes after ploidy normalisation; **chrY** uses males only.
- **Content:** per-bin median, relative MAD and genome count only (no per-genome
  data). Panel genomes are segmented against the full PoN of their stratum; the
  public filter compares only with other individuals.
- **Stratum of a sample** (`cn_stratum_of`): sample-sheet column `cn_stratum`,
  then `qc_group`, then the run's `cn_strata_table`, then
  `patients.<sample>.qc_group`, the cohort spec's `qc_group`,
  `cn_params.chemistry_stratum[<sheet chemistry>]` and finally
  `cn_params.default_stratum[<platform>]` (ONT → `ONT_1kGP_1019`, PACBIO →
  `PacBio_HiFi`).
- **Patient runs** resolve the stratum at run time (`cn_segment.py --stratum auto`)
  when the site sets none of the above: the library-QC group detected from the
  alignment (`qc/library_qc.tsv` `qc_group`; header, then read cues;
  lib/cxsv_libgroup.py, docs/LIBRARY_QC.md 4a), so the PoN and the library-QC
  calibration use the same group. If the release has no PoN for that group, the
  platform default is used. A stratum set by the site must have a PoN; otherwise
  the run stops with an error listing the release's strata. `cn/stats.json` and
  the report's summary line give the stratum and where it came from.
- **1kGP panel:** its sample sheets have no chemistry column, so the panel build
  config sets `cn_strata_table: calibration/cn_linkage/panel_strata_v1.tsv`
  (1,619 libraries; 1,017 / 220 / 247 / 135). `make_panel_strata.py` derives it
  from the lab registry's `chemistry` column and checks it against the item-14
  groups (`lofo_r2_per_genome.tsv`): all 1,619 agree.

## 3. Release and patient runs

The panel release gains three components (docs/RELEASE_AND_BUNDLE.md):
`cn/pon/<stratum>.npz`, `cn/pon/pon_members.tsv` and `cn/panel_cn_segments.tsv`.
A patient run segments against `<release>/cn/pon/<stratum>.npz` and filters
against `<release>/cn/panel_cn_segments.tsv`. A release made before the track
has no `cn/` and is refused (and its calling fingerprint differs anyway).

In the patient workflow:

- `locus_depth` (`locus_segments.py --cn-loci cn/cn_loci.tsv`) adds the
  copy-number segment boundaries of each written copy-number locus as breakends,
  so a DUP or TRP segment bounded by an undetected junction is still measured for
  criterion 5.
- `score_cxsv.py --cn-only cn/cn_only.tsv` adds the section *Copy-number-only
  candidates (not scored, not tiered)* to the report and their count to the
  summary lines. They get no points and no tier.

## 4. Fingerprint

`cn_params` is a fingerprinted section of `config/calling.yaml`; `rules/cn.smk`,
the four scripts and both libraries are in `CALLING_CODE`. A panel release and a
patient run therefore agree on the whole track. The PAR intervals are part of
`cn_params` (`par`), not a separate file. Not fingerprinted (run configuration):
`cn_depth_precomputed` and `cn_strata_table`; the PoN itself is a checksummed
release file.

## 5. Panel builds

- **Sharded build** (`workflow/panel_build`): shards stop at
  `cxsv/junction_summary.tsv` and `cn/depth.regions.bed.gz`, because everything
  after the linkage needs the aggregate PoN. The aggregate run builds the PoNs,
  segments every library and runs the linkage, the junction tables and the
  aggregation.
- **Downstream-only rebuild** (`scripts/panel_build/rebuild_downstream.sh`):
  `aggregate.yaml` must set `cn_depth_precomputed` (the full-panel depth pass,
  `<root>/<cohort>/<sample>/depth.regions.bed.gz` + `depth.ok` + `depth.sha256`)
  and `cn_strata_table`; the script stops otherwise. `cn_depth.sh` verifies the
  sha256 and bin count of each precomputed file and hard-links it; a missing file
  fails the job rather than reading the alignment on the aggregate host.

### Full-panel depth pass (cn_depth_v1)

`/mnt/sargasso1/CxSV/build/cn_depth_v1`: frozen code in `code/` (`depth_pass.sh`,
`depth_worker.sh`, `cn_depth_build.sh`, `SHA256SUMS`); `tasks.tsv` (1,619
libraries, 80.77 TB, largest first); one claim file per task; `ledger.tsv`
(`cn_depth_build.sh ledger`). 100 libraries were reused from Phase 1 by hard link
(path, size and bin count verified, sha256 recorded). The settings are the
pipeline's (`cn_depth.sh`: 1 kb, MAPQ 20, fast mode, CRAM 3.1 → 3.0 transcode with
a record-count check).

## 6. Tests

- `tests/unit/test_cn_track.py`: config = library defaults; sex and PAR ploidy;
  PoN pooling, save/load round trip and masking of absent contigs; DUP-TRP on a
  male chrX; an empty PoN gives no segments; `cn_pon.py` (flag exclusion,
  members table) and `cn_segment.py` end to end; the report's untiered section;
  copy-number boundaries in the locus-depth intervals.
- `tests/unit/test_cn_linkage.py`: the Phase 1 prototype (kept for the record).
- Synthetic CI: the chr20-only synthetic genomes have no chrX and a small panel,
  so the track runs end to end (depth, PoN, segmentation, linkage) and must leave
  the expected results unchanged.
