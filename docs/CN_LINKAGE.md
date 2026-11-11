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
2. Sex call from the chrX and chrY depth (section 2a): male, female or ambiguous.
   Ploidy track from the ploidy sex: autosomes 2; chrX 2 (female) or 1 (male,
   PAR 2); chrY 1 in males (PAR 0). An ambiguous call gets the legacy chrX rule
   (`sex_ratio_male` 0.75). A genome without a chrX track is treated as diploid.
3. Ratio to the stratum PoN median (chrX / chrY: the PoN of the genome's ploidy sex). Bins are masked where the PoN median is below
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
| `public_cn` | every segment also seen (≥ 50 % reciprocal, same direction) in ≥ `public_min_genomes` (1) other panel individuals; on chrX / chrY in ≥ `public_min_genomes_sexchrom` (2) individuals with the same sex call, sex-discordant libraries excluded (section 2b) | nothing (polymorphic) |
| `cn_only` | linked change with no junction | listed in `cn/cn_only.tsv`; the patient report lists them untiered |

**Locus coordinates (panel-v1.3, item g).** A `locus` candidate gives one CxSV locus
per contig that carries one of its segments. Its `Start`–`End` (and `Size`) is the union
of the candidate's segments on that contig and the attached breakends there
(`cxsv_cnlink.locus_span`); `cxsv_cnlink.check_locus` stops the step on a zero-length,
off-contig or out-of-segment locus. Breakends on a contig without a segment (a
translocation's mate) remain records of the locus (they are in `cxsv_only.vcf` and its
junctions) but make no locus of their own. panel-v1.2 took the span from the breakends
alone and wrote one locus per breakend contig: BOC_5883's 17p13.3 loss
(chr17:10000-2944000 and 3392000-5037000 at 1/2, joined by a translocation end at chr17:1)
was reported as `chr17:0-2`, and the mates of three candidates became 2 bp loci on chr3,
chr4, chr5, chr10, chr12 and Sniffles2's `chrN`. Re-running the linkage on the nine test
genomes (`calibration/cn_linkage/compare_cn_loci.py`, v1.2 public rule): BOC_5883 9 → 3
CN loci (8 → 2 written; chr17 now `chr17:0-5037000`), the other eight unchanged in number
and in which loci are written; every CN locus now spans its segments (for example BOC_7308
chr14:16046045-16046202 → chr14:16037000-16599000). The span feeds the scorer's `size`
feature (+0.05 at ≥ 100 kb, +0.15 at ≥ 1 Mb) and the subtelomeric context, so CN-locus
scores can rise by up to 0.15; the v1.3 design check re-scores.

A copy-number locus whose junction-derived span (`cn/cn_loci.tsv` `junction_span`, the
classifier's coordinates before the segments are added; the v1.2 id) or id equals an
existing junction-track locus is not written again (`cn/cn_loci.tsv`, `written`).

## 2. Panel-of-normals and strata

One PoN per stratum; the stratum names are the library-QC groups of item 14:
`ONT_1kGP_1019`, `ONT_LRSC_R9`, `ONT_LRSC_R10`, `PacBio_HiFi`.

- **Members:** every panel library of the stratum except the 24 library-QC flags
  (`cn_params.pon_exclude`, the list of `calibration/library_qc/flagged_v1.1.tsv`
  written into `config/calling.yaml`, so the PoN membership is fingerprinted). `pon_members.tsv` lists every library with its ploidy sex, sex
  call, chrX and chrY ratios, genome median, `in_pon` and `in_sex_reference`.
- **Autosomes** use every member. **chrX** has one reference per sex (`chrX|male`,
  `chrX|female` in the `.npz`); **chrY** uses males only (`chrY|male`). Members
  with an ambiguous sex call enter the autosomes only. (panel-v1.2 candidate
  `cxsv-panel-20261108-6001d99ea4d4` pooled both sexes on chrX; section 2a.)
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

### 2a. Sex-stratified chrX / chrY reference (panel-v1.2-r2)

**Why.** The panel-v1.2 candidate pooled both sexes in the chrX PoN. Where male
chrX depth per copy differs from female depth, males got spurious gains (state 2
at ploidy 1): ONT-1019 males 4.9 chrX gain segments per genome against 0.2 in
females. In the full-genome positive control from bundle
`cxsv-bundle-20261108-c40e26528d5f`, the PLP1 DUP-TRP L01 (a female host) was
then dropped as `public_cn`; 13 of the 14 panel genomes that matched it were
males with such gains.

**Sex call** (`cxsv_cn.call_sex`, thresholds in `cn_params.segmentation`):

| Ratio | Definition |
|---|---|
| chrX ratio | median non-PAR chrX 10 kb bin (depth > 0) / genome median |
| chrY ratio | mean non-PAR chrY 10 kb bin / genome median (most of chrY has no mappable depth, so the mean, not the median) |

- **male:** chrX ratio in `sex_x_male` [0.40, 0.60] and chrY ratio ≥ `sex_y_male_min` 0.10;
- **female:** chrX ratio in `sex_x_female` [0.80, 1.12] and chrY ratio ≤ `sex_y_female_max` 0.05;
- **ambiguous:** anything else, e.g. aneuploid-looking, mosaic or chrX-wide gains.
  Ploidy then follows the legacy chrX rule, the library is left out of the chrX/chrY
  reference, and its chrX/chrY segments count for no one in the public filter.

The thresholds come from the panel survey (`calibration/cn_linkage/sex_survey.py`
on the 1,619 cn_depth_v1 tracks). Panel medians:
- males: chrX 0.491 (robust sd 0.008), chrY 0.14–0.28;
- females: chrX 0.967 (robust sd 0.024), chrY ≤ 0.03.

**Concordance with the registry** ([sex/sex_concordance.tsv](cn_linkage/sex/sex_concordance.tsv)).
There are 0 discordant definite calls. 799 of 812 registry females are female, 790 of
794 registry males are male, and the 13 libraries without a registry sex are called
7 female, 5 male and 1 ambiguous.

The 18 ambiguous libraries are listed in
[sex/sex_ambiguous.tsv](cn_linkage/sex/sex_ambiguous.tsv):
- 13 registry females (11 with a chrX ratio of 0.50–0.76 and no chrY, i.e. possible chrX loss; 2 above the female range, including NA19238 at 1.17 with chrX-wide gains);
- 4 registry males: HG00246 with low chrY, NA19226 without chrY, HG01308 with chrX 0.61, and HG02372 with chrX 0.97 plus chrY (XXY-like);
- 1 library without a registry sex.

The legacy rule had called 13 of the low-chrX females "male".

**Design check** (`calibration/cn_linkage/sex_design_check.py`, every panel library and
the positive control re-segmented against the new PoNs built from cn_depth_v1;
[sex/design_summary_by_sex.tsv](cn_linkage/sex/design_summary_by_sex.tsv),
[sex/fig_sex_chrx_reference.png](cn_linkage/sex/fig_sex_chrx_reference.png)).

chrX gain segments per genome, pooled reference → sex-specific reference:

| | panel males | panel females |
|---|---|---|
| ONT-1019 | 4.89 → 0.01 | 0.22 → 0.56 |
| LRSC R10 | 3.08 → 1.78 | 0.31 → 0.96 |
| LRSC R9 | 3.52 → 1.39 | 1.15 → 2.31 |
| PacBio | 0.19 → 0.14 | 0.08 → 0.22 |
| all | 4.01 → 0.48 | 0.38 → 0.88 |

- **Male excess:** state 2 at ploidy 1 was 73 % of all panel chrX gains and is now 34 %.
- **Female gains rose:** about 2.3-fold. The female-only reference has a lower spread
  (relative MAD about 0.16 vs 0.20 pooled) and is centred on female depth (per copy
  about 3 % below males), so fewer bins are masked and small deviations reach the
  copy-change cut.
- **Autosomes:** identical in all 1,619 genomes.
- **chrY:** 0.11 → 0.14 segments per male genome.

L01 in the positive control: no panel genome now matches its segment under sex-matched
counting ([sex/l01_public_design.tsv](cn_linkage/sex/l01_public_design.tsv)), so it is no
longer public. Without sex matching NA19238 (ambiguous) would match.

The nested TRP level of L01 is not found against the female reference: its nested z is
5.73, below `nested_min_z` 6.0; it was 6.14 against the pooled reference. L01 is
therefore one 3-copy segment of 740 kb with its 3 INV records. `nested_min_z` is
unchanged.

### 2b. chrX / chrY public rule (panel-v1.3, decision D2)

**Why.** Under the v1.2-r2 rule one sex-matched panel library was enough to make a
chrX/chrY segment public. PC L01 (the PLP1 DUP-TRP in a female) stays private only while
the female chrX bound is 1.12: at 1.18 or 1.25 the LRSC library of NA19238 (chrX ratio
1.17, chrX-wide gains; its ONT-1019 library is 1.02) becomes female and makes L01 public
on its own.

**Rule** (`cxsv_cnlink.cn_public`, `cxsv_cnlink.sexchrom_discordant`; `cn_params.linkage`):

- a chrX / chrY segment is public only when ≥ `public_min_genomes_sexchrom` (2) other
  individuals of the sample's sex call carry it; autosomes keep `public_min_genomes` (1);
- counts are distinct individuals (two libraries of one individual count once);
- an individual with more than one panel library whose libraries disagree on sex
  (different `sex_call`, or chrX ratios differing by more than `sexchrom_discordance_x`
  0.10, or chrY ratios by more than `sexchrom_discordance_y` 0.10) counts for no one on
  chrX / chrY; its autosomal segments still count. All of its libraries are excluded,
  since depth alone cannot tell which one is wrong.

`cn_pon.py` records the decision per library in `pon/pon_members.tsv`, column
`sexchrom_public` (`yes`, or `no: discordant libraries of <ID>: <ratios>`). The PoN itself
is unchanged by this column. `cn_linkage.py --pon-members` reads that column; for a release
made before v1.3 (no column) it recomputes the same rule from `sex_call`, `x_ratio` and
`y_ratio`.

**Calibration on panel-v1.2-r2** (`pon_members.tsv`, 1,619 libraries): 35 individuals have
more than one library. The chrX ratio difference has median 0.018, p95 0.078, p99 0.135 and
maximum 0.146; the chrY ratio difference has maximum 0.046. At 0.10 two individuals are
excluded (4 libraries): NA19238 (ambiguous 1.171 vs female 1.024) and HG03499 (female
0.856 vs 0.970). At 0.05 four individuals would be excluded, at 0.15 only NA19238 (by its
sex calls). The effect on L01 and on public chrX/chrY segments per genome is measured in
the v1.3 design check (docs/PLAN_panel_v1.3.md, phase 2).

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

### After every patient run: check the sex call (interim, panel-v1.2)

Until roadmap item 16 (c) is done, the patient report prints only the ploidy sex and the
chrX ratio. An ambiguous sex call, and with it a possible sex-chromosome aneuploidy,
does not appear in the report. After every patient run, the analyst reads
`results/patient/samples/<cohort>/<sample>/cn/stats.json`:

- `sex_call` (`male`, `female` or `ambiguous`), `x_ratio`, `y_ratio`, and
  `pon_sex_reference` (the chrX/chrY reference used).
- **If `sex_call` is `ambiguous`:** flag the sample to the analyst as a possible
  sex-chromosome aneuploidy or mosaicism (for example 45,X, 47,XXY, 47,XXX or
  mosaic X loss), and treat every chrX and chrY copy-number call of that run as
  unreliable. In panel-v1.2 an ambiguous sample is segmented on chrX with the
  reference of its ploidy sex (`sex_ratio_male` 0.75) and is never counted for
  sex-matched panel frequency. A 45,X sample gets male ploidy and no chrX loss. An
  XXY sample gets female ploidy and its chrY is not segmented. Neither case is
  reported as a finding.
- `score_cxsv.py` infers sex separately from flank depth (`infer_sex`, cut-offs
  0.75 / 0.85; a sample-sheet sex overrides it). That call can disagree with
  `sex_call`. When the two disagree, `cn/stats.json` is the reference.

Panel evidence that this case is real: 18 of 1,619 panel libraries are ambiguous.
Nine of them are registry females with a male-range chrX ratio (0.50-0.66) and
no chrY, consistent with X loss in cell lines or a registry error. Three are
females whose chrX ratio exceeds the female range: NA18946 1.15; NA19238 LRSC
1.17, whose other library (ONT-1019) is normal at 1.02; and HG02805 1.26
(docs/cn_linkage/sex/sex_ambiguous.tsv).

## 4. Fingerprint

`cn_params` is a fingerprinted section of `config/calling.yaml`; `rules/cn.smk`,
the four scripts and both libraries are in `CALLING_CODE`. A panel release and a
patient run therefore agree on the whole track. The PAR intervals are part of
`cn_params` (`par`), not a separate file. Not fingerprinted (run configuration):
`cn_depth_precomputed` and `cn_strata_table`; the PoN itself is a checksummed
release file.

## 5. Panel builds

- **Sharded build** (`workflow/panel_build`): shards stop at
  `cxsv/junction_summary.tsv`, `cxsv/junction_only.vcf` and `cn/depth.regions.bed.gz` (`sv_summary.txt`
  moves to the aggregate run), because everything
  after the linkage needs the aggregate PoN. The aggregate run builds the PoNs,
  segments every library and runs the linkage, the junction tables and the
  aggregation.
- **Downstream-only rebuild** (`scripts/panel_build/rebuild_downstream.sh`):
  `aggregate.yaml` must set `cn_depth_precomputed` (the full-panel depth pass,
  `<root>/<cohort>/<sample>/depth.regions.bed.gz` + `depth.ok` + `depth.sha256`)
  and `cn_strata_table`; the script stops otherwise. `cn_depth.sh` verifies the
  sha256 and bin count of each precomputed file and copies it (`cp --reflink=auto`),
  then checks the copy's sha256 against the recorded one; a missing file fails the job
  rather than reading the alignment on the aggregate host. panel-v1.3 (item a): the
  source is never hard-linked, so a rebuild cannot change `cn_depth_v1` mtimes (in v1.2
  the shared inodes changed 1,619 of them and stopped the run once); the precomputed
  depth costs ~1,619 × the size of one 1 kb depth file of disk per build.

### Full-panel depth pass (cn_depth_v1)

`/mnt/sargasso1/CxSV/build/cn_depth_v1`: frozen code in `code/` (`depth_pass.sh`,
`depth_worker.sh`, `cn_depth_build.sh`, `SHA256SUMS`); `tasks.tsv` (1,619
libraries, 80.77 TB, largest first); one claim file per task; `ledger.tsv`
(`cn_depth_build.sh ledger`). 100 libraries were reused from Phase 1 by hard link
(path, size and bin count verified, sha256 recorded). The settings are the
pipeline's (`cn_depth.sh`: 1 kb, MAPQ 20, fast mode, CRAM 3.1 → 3.0 transcode with
a record-count check).

## 6. Tests

- `tests/unit/test_cn_sex.py`: sex calls (panel survey values, ambiguous cases),
  one chrX reference per sex with ambiguous members on autosomes only, a male-only
  depth artefact absorbed by the male reference, the legacy pooled PoN still read,
  sex-matched public counting, the panel table's sex column; panel-v1.3 D2: two
  sex-matched individuals on chrX/chrY (one on autosomes), distinct-individual counting,
  sex-discordant libraries (panel values of NA19238 and HG03499) excluded on chrX/chrY
  only, `pon_members.tsv` `sexchrom_public` written by `cn_pon.py` and read (or recomputed
  for an older table) by `cn_linkage.py`.
- `tests/unit/test_cn_locus_coords.py`: `check_locus` rejects zero-length, out-of-segment,
  mate-contig and past-the-end loci; the BOC_5883 structure end to end through
  `cn_linkage.py` gives one locus `chr17:0-5037000` and none on the mate contig; a CN locus
  with a junction-track locus's span is not written again.
- `tests/unit/test_cn_depth_precomputed.py`: precomputed depth is copied, the source's
  mtime, inode and link count are unchanged, a wrong sha256 or bin count fails.
- `tests/unit/test_cn_track.py`: config = library defaults; sex and PAR ploidy;
  chrX reference per sex, save/load round trip and masking of absent contigs; DUP-TRP on a
  male chrX; an empty PoN gives no segments; `cn_pon.py` (flag exclusion,
  members table) and `cn_segment.py` end to end; the report's untiered section;
  copy-number boundaries in the locus-depth intervals.
- `tests/unit/test_cn_linkage.py`: the Phase 1 prototype (kept for the record).
- Synthetic CI: the chr20-only synthetic genomes have no chrX and a small panel,
  so the track runs end to end (depth, PoN, segmentation, linkage) and must leave
  the expected results unchanged.
