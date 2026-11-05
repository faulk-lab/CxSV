# Oversized clusters: the long-record blind spot (roadmap item 13)

Status 2026-11-04: implemented and evaluated on panel-v1.0 data. *Update
2026-11-05:* the panel was rebuilt with this rule as panel-v1.1
(`…64e10c1b1363`); healthy Tier 1 rose to 0.851 per genome, mostly from SVIM-only
low-quality breakend records in 11 noisy genomes that this rule had exposed
(they were hidden inside oversized clusters before). That build is superseded;
panel-v1.1 is re-released with `single_caller_bnd_rule: svim_lowq` from a new
build root `panel-v1.1-r2` (docs/R10_EXCESS.md).
History: the first version (`de23a0e`, merged into `main` as PR #20 before the
evaluation finished) clustered long records at their breakends. The evaluation
below replaced that with leaving them out of discovery; the follow-up is branch
`item13-oversized-clusters-2` (on `main` 359f6d7). `de23a0e` alone is not the
recommended state: it changes the fingerprint and raises healthy Tier-1 rates
(section 4b). Every number labelled *provisional* below
rests on a 100-genome subset re-scored against the panel-v1.0 junction sites,
not on a rebuilt panel.

Code: `calibration/oversized_clusters/` (analysis), `workflow/scripts/vcf_to_bed.py`,
`workflow/scripts/extract_cxsv_vcf.py`, `workflow/lib/cxsv_junctions.py`,
`workflow/rules/per_sample.smk`, `workflow/rules/fingerprint.smk`,
`config/calling.yaml` (`cxsv_params.span_link_max`, `long_span_rule`). Tables:
`docs/item13/`. Figures: `docs/figures/item13/`.

## 1. The problem

Discovery writes every merged SV of a sample to `merged_sv.bed` and clusters it
with `bedtools cluster -d 50000` (`cluster_distance`). A DEL/DUP/INV is written
as one interval POS..END, so a record that spans 25 Mb joins every SV inside
those 25 Mb into one cluster. `classify_cxsv.py` drops any cluster longer than
`max_cluster_size` (500 kb) whole (`FILTERED_oversized_cluster`). In the
whole-genome positive control a 25.7 Mb Sniffles2 DEL across the chr17
pericentromere (16.7-42.4 Mb) put 224 records into one 25.8 Mb cluster and the
planted *NF1* event E02 was lost. Panel-wide, 1,616 of 1,619 genomes lost a
median 334 Mb this way.

## 2. Characterisation (all 1,619 panel-v1.0 genomes)

`characterise.py` re-implements the clustering in Python and reproduces the
released `clusters.bed` partition exactly in **1,619/1,619** genomes, then
describes every oversized cluster and every record longer than 500 kb
(`collect.py` → `docs/item13/*.tsv`).

**Mechanism** (`oversized_cluster_mechanisms.tsv`). 22,541 oversized clusters.
99.5% of their summed span is in the 18,165 clusters that contain at least one
record longer than 500 kb; dense runs of short SVs (2,889 clusters) and
clusters bridged only by 50-500 kb records (1,487) make up 0.3% and 0.2%.

**The long records** (`giant_record_composition.tsv`, `giant_records_by_platform.tsv`,
`giant_record_sites.tsv`, figure `fig_chaining_records.png`). 48,773 records
longer than 500 kb (median 8.5 Mb, up to 248 Mb; 18.6 per ONT 1kGP genome, 41-61
per LRSC / PacBio genome):

- 79% (38,653) are breakend pairs: the caller wrote a BND with a mate on the
  same chromosome and SURVIVOR rewrote it as `SVTYPE=INV` with `END` = mate
  position, so it is clustered as an interval although it joins two distant
  points and spans nothing;
- 89% have one supporting caller (Sniffles2 only 23,945, SVIM only 19,560);
  the 5,268 two-caller records are all SVIM + Sniffles2. cuteSV supports none
  (its default `--max_size` is 100 kb), so "all three callers" is impossible
  above 100 kb;
- they are recurrent: 62% lie at sites (ends within 100 kb) seen in >= 100 panel
  genomes, 22% at 10-99, 9% in one genome. The top site, chr1:16.4-145.2 Mb (a
  128.8 Mb breakend pair), is in 1,493 genomes; the E02 chainer
  chr17:16.7-42.4 Mb is a DEL in 645 genomes and a breakend pair in 709; the
  17q21.31 inversion polymorphism (chr17:46.45-47.15 Mb) is in 741;
- 56% have both breakends within 10 kb of a segmental duplication or satellite
  (25% one end);
- **read depth** (`depth_check.py`, `giant_record_depth_tests.tsv`): for 139
  DEL/DUP sampled across caller support, recurrence and cohort, read counts in
  20 kb windows inside vs flanking the record, each relative to a same-cohort
  genome without the call (110 informative). Only 23% of DELs (ratio <= 0.75)
  and 27% of DUPs (>= 1.25) change copy number. Caller count does not predict it
  (two-caller DEL 1 of 15, single-caller DEL 31%; DUP 28% vs 26%); rare records
  are supported more often than recurrent ones (DEL 44% vs 13%).

So the records that blind discovery are mostly recurrent reference/paralogue
artefacts, mostly without a copy-number change; a minority are real large
variants (e.g. 17q21.31), and caller support cannot tell them apart.

**What the oversized clusters hid**: 181,722 CxSV loci over the panel (112 per
genome) are found inside them once the long records stop chaining.

## 3. Design

Candidates were simulated on every panel genome with the unchanged classifier
(`characterise.py` VARIANTS; `rule_comparison_panelwide.tsv`):

| rule | blind span per genome, median (mean) | genomes with any | new CxSV loci / genome | of which hold a > 500 kb junction | loci lost |
|---|---|---|---|---|---|
| panel-v1.0 (interval) | 334 Mb (359) | 1,616 | – | – | – |
| long records as two breakends (`breakends`) | 1.64 Mb (2.80) | 1,284 | 122 | 20.8 | 0 |
| long records left out (`exclude`, **adopted**) | 1.40 Mb (2.70) | 1,245 | 112 | 0 | 0 |
| `breakends`, single-caller long records left out | 1.40 Mb (2.70) | 1,245 | 114 | 3.2 | 0 |
| also BND-notation pairs of any length as breakends | 1.55 Mb (2.67) | 1,270 | 125 | 20.9 | 1.3 |
| every record > 50 kb as breakends | 1.07 Mb (2.01) | 1,133 | 135 | 21.9 | 4.9 |

The threshold is `max_cluster_size` itself: a record longer than it would make
any cluster containing it oversized on its own, so under either rule **no
cluster that was classified before changes** (0 loci lost; also checked by
`test_clusters_that_are_not_oversized_are_unchanged`). Rules that touch
shorter records (last two rows) change classified loci and were not pursued.

Between the two minimal rules the choice was pre-declared as "the one with the
lower healthy-genome Tier-1 rate" (section 4b). `breakends` keeps the long
record as a junction of the loci at its ends; the scorer then counts genes and
read depth over the record's whole interval (up to 129 Mb), which inflates
`gene_dosage` / `gene_count` points and, in two noisy genomes, asked
`samtools bedcov` for 10-18 Gb of segments. `exclude` gave the lower rate and is
adopted:

> A DEL/DUP/INV longer than `cxsv_params.span_link_max` (= `max_cluster_size`,
> 500 kb) takes no part in CxSV discovery (`long_span_rule: exclude`). It stays a
> junction of the sample (`junctions/junctions.tsv`, so panel junction sites and
> privacy matching are unchanged) and a merged call (`merged_sv.vcf`), and is
> listed in `cxsv/long_span_records.tsv`. Before item 13 any cluster containing
> such a record was dropped, so no locus that panel-v1.0 classified loses
> anything. `extract_cxsv_vcf.py` no longer puts such records, nor
> inter-chromosomal records whose breakends lie outside every locus, into
> `cxsv_only.vcf` by interval overlap.

`long_span_rule: breakends` remains available (documented, unit-tested) but is
not used. What stays oversized are dense runs of short SVs in satellite,
acrocentric and subtelomeric sequence (median 1.4 Mb per genome; most common at
chr14:4 Mb, chr8:0, chr5:0, chr22:2, chr20:32, chr9:48 Mb); they are still not
classified.

Also fixed (needed for a correct rebuild): `calling_fingerprint` had no code or
config inputs, so in an existing work directory or build tree it was never
re-made after a code change. A release rebuilt in place, or a patient run
re-used with new code, would have carried the old fingerprint (the first
positive-control check here reported `COMPATIBLE` with the stale value). The
rule now takes the fingerprinted code files and `config/calling.yaml` as inputs
and the calling parameters as params.

## 4. Evaluation

**(a) Blind span, all 1,619 genomes** (re-clustered with the adopted rule;
`blind_span_per_genome.tsv`, `blind_span_by_stratum.tsv`, figure
`fig_blind_span.png`): median 334 → 1.40 Mb per genome. ONT 1kGP 280 → 0.57,
ONT LRSC R9 414 → 4.7, LRSC R10 509 → 7.1, PacBio HiFi 415 → 3.7 Mb; by coverage
<10x 206 → 0, 10-20x 272 → 0.53, 20-30x 323 → 1.18, >= 30x 440 → 5.2 Mb. The
simulation equals the real pipeline: on the 100-genome subset the re-run
workflow gave the same oversized span in every genome.

**(b) 100 panel genomes, provisional** (`eval_subset_100.tsv`: 40 ONT 1kGP, 20
LRSC R10, 20 LRSC R9, 20 PacBio, spread over each group's coverage range).
`downstream_rerun.sh` copied each genome's panel-v1.0 outputs and re-ran the
workflow with the branch code:

- Snakemake planned and ran only the 7 downstream rules per sample
  (`sample_vcf_to_bed`, `sample_cluster_sv`, `sample_classify_cxsv`,
  `sample_junctions`, `sample_extract_cxsv_vcf`, `sample_annotate_cxsv`,
  `sample_sv_summary`); no caller, filter, merge or callability job;
- caller VCFs and `merged_sv.vcf` byte-identical before/after (sha256), and the
  junction table identical except the `locus_id` column in 100/100 genomes;
- CxSV loci +126 per genome (median), oversized span median 381 → 3.2 Mb;
- cost: 7.3 s wall per genome (sum of rule benchmarks), at most 1.9 GB RSS;
  700 jobs in 97 s on 16 cores.

Healthy-genome Tier rates (`eval_scores.sh`: the item-10 LOFO approach, each
genome scored by the unchanged scorer against the released panel-v1.0 junction
sites with its family removed; read depth for non-public loci only;
`lofo_tier_rates_summary.tsv`, figure `fig_tier_rates.png`):

| genomes | rule | Tier 1 / genome | Tier 2 / genome | genomes with any Tier 1 |
|---|---|---|---|---|
| 97 (3 QC outliers removed) | panel-v1.0 | 0.196 | 2.84 | 17.5% |
| | breakends | 0.371 | 3.52 | 30.9% |
| | **exclude** | **0.237** | **3.14** | 20.6% |
| all 100 | panel-v1.0 | 1.97 | 12.8 | 20% |
| | **exclude** | 4.76 | 27.2 | 23% |

93 of 100 genomes keep exactly their Tier-1 count. Three genomes that were
already QC outliers under panel-v1.0 (HG04176 82 → 266 Tier-1 calls, HG03680
59 → 144, HG00268 37 → 43; the first two had 1.7 and 1.5 Gb in oversized clusters)
account for 275 of the 280 added Tier-1 calls. The breakends arm lacks HG04176
and HG03680 (their depth step was stopped after the segment set reached
10-18 Gb). These rates are provisional: the panel is not rebuilt (the junction
sites, and therefore privacy, do not change with clustering; region-overlap
carriers, which are not scored, still come from panel-v1.0 loci), and 100 genomes
is a subset.

**(c) Positive control, check only** (`pc_check.sh`: the bundle's whole-genome
run re-used with its Snakemake metadata, branch code bound over the bundle repo,
panel-v1.0 release; `pos_control_events_before_after.tsv`, figure
`fig_pos_control.png`). E02 is now called: private, `Complex_inversion`, Tier 1,
1.95 points (10-slice regression 1.97). E01 and E03-E12 keep class, privacy,
tier and score; E13 stays private (panel composition, as before). Host genome:
0 Tier-1 and 4 Tier-2 calls, as under panel-v1.0 (the breakends rule gave 3
Tier-1 host calls, each at an end of a > 500 kb record). Junctions 40/41 exact.
1,246 loci scored (1,065 before).

## 5. Rebuild plan (as proposed 2026-11-04; as done: docs/PANEL_BUILD.md, "Downstream-only rebuild")

The panel can be rebuilt from the panel-v1.0 caller outputs by re-running only
downstream per-sample rules and the aggregate.

1. Merge `item13-oversized-clusters-2` into `main` (main 359f6d7 already has
   `de23a0e`, item 11 and the small fixes); tag decision `panel-v1.1`.
2. New build root `/mnt/sargasso1/CxSV/build/panel-v1.1` (the release path is
   `<results>/release`, so a rebuild inside panel-v1.0 would overwrite its
   release). Per-sample upstream outputs (callers, filtered VCFs,
   `merged_sv.vcf`, callability, `.snf`) and the cohort / Sniffles2-population
   outputs are **symlinked** from panel-v1.0 (copying would move ~1 TB); only
   downstream outputs are written as new files. Update `aggregate.yaml`
   (`panel_results_root`, release description) and the fleet config's
   `expected_fingerprint`.
3. Dry run, then run, on freya:
   `snakemake -s <main>/workflow/panel/Snakefile --directory BUILD2/aggregate
   --configfile BUILD2/aggregate.yaml --forcerun sample_vcf_to_bed calling_fingerprint
   --allowed-rules <every rule except run_*, filter_*, merge_sv_callers,
   check_cram_reference, callable_mask, panel_sniffles_population>`.
   The work directory has no per-sample Snakemake metadata (as in panel-v1.0's
   `aggregate/`), so `--forcerun` is what makes the downstream rules run;
   `--allowed-rules` guarantees that no caller or mosdepth job can.
4. Positive-control check against the new release, bundle build and offline
   validation (as for panel-v1.0).

Planned jobs (dry runs of `main` 359f6d7 + this branch on the full panel-v1.0
tree, read-only; 0 files under `release/` or `panel/` changed), with
`--forcerun sample_vcf_to_bed calling_fingerprint`:

- **with `callable_mask` held** (`--allowed-rules` without the caller, filter,
  merge, `check_cram_reference`, `callable_mask`, `panel_sniffles_population`
  and `build_ref_cache` rules): **12,975 jobs** = 7 × 1,619 downstream
  per-sample jobs + `check_svim_output` × 1,619 (small fixes: new output) +
  `calling_fingerprint`, `panel_release_id` + 21 aggregate jobs
  (`cohort_cxsv_tables` ×3, `panel_sample_manifest`, `build_panel`,
  `build_panel_sites`, `plot_complexity_class_heatmap` ×3, `plot_cxsv_venn` ×3,
  `plot_genomic_context` ×3, `panel_site_callability`, `panel_sniffles_sites`,
  `plot_panel_sharing`, `panel_release_metadata`, `panel_release`, `all`). No
  caller, merge, mosdepth or Sniffles2 population job;
- **without holding it: 14,595 jobs**, adding `callable_mask` × 1,619 and
  `panel_callability_masks`. The small fixes made `callable_mask.sh` an input of
  `callable_mask`, newer than every panel-v1.0 mask ("Updated input files";
  `--rerun-triggers` without `mtime` does not stop it).

(Before the small fixes were merged, the same dry run planned 11,354 jobs.)

| step | measured basis | estimate |
|---|---|---|
| downstream per-sample rules | 7.3 s wall / genome; 100 genomes in 97 s on 16 cores | 1,619 genomes ≈ 30 min on 16 cores (NAS-bound, up to ~1 h) |
| `check_svim_output` (small fixes) | measured 0.9-3.9 s per genome (3 panel genomes, one per platform) | ≈ 1.1 CPU-h, ≈ 5-10 min on 16 cores |
| aggregate rules | panel-v1.0 aggregate logs (`rule_times.py`): manifest 0.41 h, cohort tables 0.37 h, build_panel_sites 0.85 h, site callability 0.06 h, Sniffles2 sites 0.54 h, plots ≤ 0.17 h, release 0.09 h | critical path ≈ 2.0 h |
| symlink farm + dry run | – | ≈ 15 min |
| **total, freya only** | | **≈ 3-3.5 h wall** |
| optional: `callable_mask` not held | panel-v1.0: 220 CPU-h, 76 job-wall-h, all ~80 TB of alignments re-read (NAS floor ~36 h) | avoid; before holding it, run small-fixes' `callable_mask.sh` on 6 genomes (2 ONT CRAM 3.0, 1 CRAM 3.1, 1 LRSC BAM, 1 PacBio) and require byte-identical masks (≈ 1-2 h) |
| positive-control check, bundle build + validation | panel-v1.0 validation: PC 25-35 min, bundle 90 s, validation 2.5 h | ≈ 3 h |

Artefacts that change: every genome's `cxsv/{merged_sv.bed, merged_sv_sorted.bed,
clusters.bed, cxsv_summary.tsv, cxsv_only.vcf, cxsv_master_table.tsv,
cxsv_loci.bed}`, the new `cxsv/long_span_records.tsv`, `junctions/junctions.tsv`
(`locus_id` column only), `sv_summary.txt`; cohort CxSV tables and plots;
`panel/panel_cxsv_calls.tsv`, `panel_loci.tsv`, `panel_cxsv_sharing.bed`,
`panel_sample_manifest.tsv`; the junction-site files are rewritten from the same
junctions (content expected identical apart from the version column);
`panel_site_callability`, Sniffles2 site comparison, `calling_fingerprint.json`
(new fingerprint, see below), `release.json`, the release
manifest and id. Unchanged: caller VCFs, `merged_sv.vcf`, callability masks
(if held), Sniffles2 population VCF. After the rebuild the item-10 graded
privacy calibration should be re-checked on the new panel (Tier rates in 4b
are provisional).

## 6. State of `main` and fingerprints

`de23a0e` (merged as PR #20) is correct as far as it goes (unit tests, the
LONGINV synthetic case and CI pass; no classified locus is lost) but it is
**not the recommended rule** and is incomplete:

- it links records longer than 500 kb at their breakends, which raises healthy
  Tier-1 rates (0.196 → 0.371 per genome on 97 panel genomes) and gives the
  positive-control host 3 Tier-1 calls (section 4);
- it does not list the long records anywhere (`long_span_records.tsv` is new);
- `calling_fingerprint` is still not re-made after a code change in an
  existing work directory or build tree, so a rebuild in place or a re-used
  patient run would carry a stale fingerprint.

`item13-oversized-clusters-2` (on `main` 359f6d7) fixes all three. Calling
fingerprints (computed with `cxsv_fingerprint.build` from each tree's
`workflow/` and `config/calling.yaml`, the hs1 MD5 manifest and the pinned
image's tool versions):

| code | fingerprint | code digest |
|---|---|---|
| panel-v1.0 (release) | `cd7e1d194cd4…` (reproduced) | `2320a190213d` |
| `de23a0e` = `main` 359f6d7 | `e6a91994aaa5…` | `f07d3d10b7da` |
| `item13-oversized-clusters-2` | `64e10c1b1363…` | `32b14645be1e` |

`main` as it stands is incompatible with panel-v1.0 (any patient run on it
needs a panel built with the same code), and no panel exists for `e6a91994aaa5`.
The next release should be built from `item13-oversized-clusters-2` merged.
