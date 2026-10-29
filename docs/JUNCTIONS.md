# Junctions, panel junction sites and breakpoint-level matching

Roadmap item 3. A patient CxSV locus is called **public** only when panel
individuals carry *the same rearrangement*, judged junction by junction.
Region overlap is kept only as an informational hotspot flag.

## 1. Junction model (per sample; fingerprinted)

A **junction** is one novel adjacency: two reference positions joined by a
structural variant. Every record of a sample's `merged_sv.vcf` (SURVIVOR merge
of cuteSV, SVIM and Sniffles2) becomes one junction
(`workflow/lib/cxsv_junctions.py`, `workflow/scripts/extract_junctions.py`,
rule `sample_junctions` in `workflow/rules/junctions.smk`):

| svtype | chrom2, pos2 | svlen |
|---|---|---|
| DEL, DUP, INV | chrom1, INFO/END | \|SVLEN\| (or END − POS) |
| INS | chrom1, POS | inserted length |
| BND (SURVIVOR `TRA` is written BND) | the mate, from the ALT brackets (or CHR2/END) | 0 |

Ends are ordered so that (chrom1, pos1) ≤ (chrom2, pos2) in the reference
`.fai` contig order, then by position. Positions are 1-based VCF positions.

**Collapsing (for matching only).** Some junctions appear twice in a merged
VCF. Sniffles2 writes each inversion breakpoint as two reciprocal breakend
records that SURVIVOR keeps (e.g. INV 51,213,665→51,223,040 and INV
51,223,041→51,213,667 at the *SHANK3* event E04), and a BND can appear with its
mate. Records of one sample are one junction when svtype and contigs agree and
both ends lie within `junction_params.collapse_distance` (300 bp, the same scale
as SURVIVOR's `max_dist`), and, for INS, when lengths are within
`ins_len_ratio` (0.5). Collapsing is leader-based: the best-supported record
(most callers) claims its neighbours, which are never compared with each
other, so it cannot chain. SUPP_VEC is OR-ed and genotype dosage is the maximum.
`classify_cxsv.py` and `annotate_cxsv.py` are unchanged: the C1–C5 criteria
still see every record. Fixing E12 there is roadmap item 5.

**Locus membership.** `cxsv/clusters.bed` says which cluster each record
belongs to, and `cxsv/cxsv_summary.tsv` says which clusters are CxSVs. A
junction's `locus_id` is `Chr:Start-End` exactly as in `cxsv_master_table.tsv`.
A translocation whose breakends sit in two CxSV loci belongs to both.

Output: `samples/<cohort>/<sample>/junctions/junctions.tsv`. It covers the whole
SV set, not only CxSV loci, because the panel must count carriers of every
constituent SV.

| column | meaning |
|---|---|
| chrom1 pos1 chrom2 pos2 svtype svlen | the junction |
| locus_id | CxSV locus/loci (comma-separated) or `.` |
| supp_vec | callers (cuteSV, SVIM, Sniffles2), OR-ed over collapsed records |
| dosage | 1 or 2, from the first called GT among Sniffles2, cuteSV, SVIM |
| n_records, record_ids | merged_sv.vcf records collapsed into it |

The junction table defines what is matched, so `rules/junctions.smk`, the
script, the library and `junction_params` are part of the calling fingerprint.
A patient run must build its junctions exactly as the panel did.

## 2. Shared junction-site schema

All site tables use this schema: the panel's clustered sites, the Sniffles2
population sites, and the external catalogues (roadmap item 4). Files are
tab-separated with a header row, 1-based positions, and UCSC hs1 contig names
chr1..chrM.

```
chrom1  pos1  chrom2  pos2  svtype  svlen  ac  an  af  source  version   [af_<group> ...]
```

- `svtype` ∈ {DEL, DUP, INV, INS, BND}; translocations are BND.
- DEL/DUP/INV: chrom2 = chrom1 and pos2 = END. INS: pos2 = pos1 and svlen is
  the inserted length. BND: (chrom2, pos2) is the mate.
- Each site is written once, with (chrom1, pos1) ≤ (chrom2, pos2) in hs1 .fai
  contig order, then position.
- `ac`, `an` and `af` may be `.` when a source does not provide them.
- Optional `af_<group>` columns hold per-group allele frequencies (platform,
  superpopulation).

Conventions agreed with the catalogue track (item 4,
`acquisition/catalogues/junction_sites.py` on branch
`item4-catalogue-acquisition`), where the schema itself is silent:

- `svlen` is unsigned.
- BND `svlen` is `.`. Only a catalogue that gives the outer span of an
  intra-chromosomal complex SV writes a number there; the panel's BNDs are
  breakend pairs.
- Orientation and row order follow hs1 `.fai` contig order, which is by
  length (… chr7, chrX, chr9, chr8 …), not numeric order.
- `af` is the source's AF when given, else ac/an. Sniffles2 population VCFs
  carry no INFO/AF.
- Sites of non-BND SVs shorter than 50 bp are dropped. 50 bp equals the
  pipeline's `sv_params.min_sv_len`, which the rules pass as `--min-svlen`.

The per-sample junction tables are internal and keep a numeric `svlen` (0 for
BND).

For `source=cxsv_panel`, `ac` is the sum of carrier dosages (het 1, hom 2,
capped by ploidy). `an` is the sum of ploidy over all panel individuals, with
chrX/chrY hemizygous in males; PAR is not modelled, and sex comes from the
registry `individuals.tsv`. `af` = ac/an. `version` is `panel_release_id` or
`fp-<first 12 hex of the calling fingerprint>`.

## 3. Panel junction sites (panel workflow)

`workflow/scripts/build_panel_sites.py` (rule `build_panel_sites`) pools every
panel sample's junction table and clusters the junctions across individuals
(`workflow/lib/cxsv_matching.py: cluster_sites`):

- **Type-aware:** svtype, chrom1 and chrom2 must agree.
- **Leader-based within `panel_sites.tolerance` (500 bp):** observations are
  visited in position order. Each joins the nearest open site whose first
  member has both ends within the tolerance; otherwise it starts a new site. A
  site never spans more than the tolerance at either end, and sites cannot
  chain. INS also need lengths within `ins_len_ratio`.
- **Carriers are distinct individuals:** one person sequenced in two cohorts
  counts once.

Release files:

| file | content |
|---|---|
| `panel_junction_sites.tsv` | shared schema, `source=cxsv_panel`, + `af_<platform>`, `af_<superpopulation>` |
| `panel_junction_carriers.tsv` | per site: envelope (min/max pos1, pos2, svlen), n_observations, n_carriers, carrier individual ids; the matcher reads this |
| `panel_individuals.tsv` | individual_id, sex, superpopulation, platforms, cohorts, samples (denominators and breakdowns) |
| `sniffles_population_sites.tsv` | Sniffles2 population-mode genotypes, same schema, `source=sniffles2_population` |
| `panel_frequency_comparison.tsv`, `..._summary.tsv` | per Sniffles2 site: matched clustered site and carriers on both sides |

**Sniffles2 population mode** (rule `panel_sniffles_population`) runs
`sniffles --input snf_inputs.tsv --vcf sniffles_population.vcf` on every panel
sample's `.snf` and re-genotypes each candidate in every sample.
`sniffles_population_sites.py` converts the result to the shared schema. Its
`an` counts only called (non-missing) genotypes. It also compares each site's
carriers with those of the clustered site it matches. The Sniffles2 sites are
kept alongside the clustered sites; matching uses the clustered sites, because
they come from the same three-caller merged calls as the patient's junctions.

### Sniffles2 population mode vs clustered sites (lab regression panel, 2026-10-28)

The panel is 10 ONT background genomes (34.4 Mb slices, including E13 in 3 of
them; workdir `cxsv_item3_regression_v2`). Sniffles2 2.8.1 population mode
genotyped 9,770 sample-sites: 9,328 called, 442 missing, 2,824 non-reference.

| | DEL | INS | INV | DUP | BND |
|---|---|---|---|---|---|
| Clustered sites | 402 | 477 | 13 | 3 | 7 |
| Sniffles2 sites (matched to a clustered site / total) | 310/314 | 634/640 | 8/8 | 2/2 | 0/0 |
| Clustered sites matched by >1 Sniffles2 site | 0 | 129 | 0 | 0 | - |
| Carrier calls, Sniffles2 vs clustered (summed per clustered site) | 1,241 vs 1,203 | 1,473 vs 1,458 | 36 vs 36 | 2 vs 2 | - |
| Sites with equal carrier counts (one-to-one matches) | 236/310 | 209/257 | 8/8 | 2/2 | - |

Readout:

- **Totals agree within 3%.** Summed carrier calls differ by 3% for DEL and
  1% for INS.
- **Per-site differences are small.** The mean absolute difference over all
  matched Sniffles2 sites is 1.65 carriers. Most of it comes from Sniffles2
  splitting one INS position into sequence-distinct alleles
  (`--combine-pctseq 0.7`), whereas the clustering uses position and length
  only.
- **No BND sites from Sniffles2.** It reported none, so translocation
  frequencies come from the clustered sites alone.
- **Unmatched clustered sites.** 191 of 902 clustered sites (keyed by
  coordinates: DEL 92, INS 80, BND 7, INV 5, DUP 1) have no Sniffles2
  counterpart. Most are cuteSV- or SVIM-only calls.

Matching therefore stays on the clustered sites, which are built from the same
three-caller calls as the patient's junctions. The Sniffles2 frequencies are
kept alongside as a second opinion.

## 4. Patient matching (patient workflow)

`workflow/scripts/classify_patient_cxsv.py` (rule `classify_patient_cxsv`;
settings in `config/patient.yaml: panel_match`, not fingerprinted):

1. A patient junction **matches** a panel site when svtype and contigs agree
   and each end lies within `tolerance` (500 bp) of the site's envelope. INS
   also need lengths within `ins_len_ratio`. The junction's carriers are the
   union of the carriers of the matched sites.
2. A panel individual **carries the locus** when it carries at least
   `min_junction_fraction` (0.75) of the locus's junctions.
3. Status:
   - `public`: at least `min_carriers` (1) individuals carry the rearrangement;
   - `shared_junctions`: some junctions occur in the panel, but no individual
     (or too few) carries the rearrangement. Examples are a different
     rearrangement, or a common simple SV, at the same place. **Not public**;
   - `private`: none of the locus's junctions occurs in the panel.
4. **Hotspot flag (information only):** `Panel_Region_Carriers` counts the panel
   individuals with any CxSV call within `hotspot_window` (1 kb) of the locus,
   which was the previous matching rule. `Panel_Hotspot = yes` when that count
   reaches `min_carriers` but the locus is not public, i.e. region overlap would
   have masked it.

New columns in `private/cxsv_panel_status.tsv`: Panel_Junctions,
Panel_Junctions_Seen, Panel_Carriers, Panel_Carrier_Freq,
Panel_Carriers_By_Platform, Panel_Carriers_By_Superpopulation,
Panel_Best_Junction_Fraction, Panel_Region_Carriers, Panel_Hotspot and
Panel_Status. `private/private_cxsv_master_table.tsv` and the private VCF hold
every locus that is not public. `private/junction_panel_status.tsv` lists every
patient junction, including simple SVs outside CxSV loci, with its panel sites,
carriers, carrier frequency and best site allele frequency.

### Choice of `min_junction_fraction` = 0.75

- **Small loci need every junction.** 0.75 requires all junctions for loci of
  ≤ 3 junctions, where one junction can be the whole difference between two
  rearrangements. From 4 junctions on, one unmatched junction is tolerated.
- **Caller artefacts make "all junctions" too strict.** In the synthetic test
  the patient's COMMON CxSV, which all 3 panel genomes carry identically,
  gained a cuteSV-only DEL inside its inversion (chr20:210,001–216,946,
  SUPP_VEC 100). No panel genome has that call, so the best fraction is 3/4 and
  the locus fails at 0.8.
- **Single-caller junctions cannot simply be ignored.** 13 of the 37 planted
  junctions of the positive control are recovered by Sniffles2 alone
  (SUPP_VEC 001; all 3 of E06's). Discounting single-caller junctions would
  drop true, distinguishing junctions and mask private events.
- **Leave-one-out.** On the 10 background genomes of the lab regression
  (pre-item-3 calls), each genome's CxSV loci (164 in total) were matched
  against sites built from the other 9. 105 had a full match in another
  genome, 135 had one at ≥ 0.75, 122 at ≥ 0.8, 154 at ≥ 0.67 and 161 at ≥ 0.5;
  80 of the 164 loci have ≤ 3 junctions. Most are clusters of common
  subtelomeric deletions.
- **Lowering further is the risky direction.** A private event sharing one
  50 kb cluster with common SVs would inherit their carriers. Those cases are
  `shared_junctions`, which stays reviewable.

### Forced genotyping (Sniffles2 `--genotype-vcf`): not implemented

The roadmap suggested force-genotyping panel sites near patient loci in the
patient. It is not implemented, for three reasons.

- **It cannot change a status.** Status depends on whether *panel individuals*
  carry the *patient's* junctions. A panel site that the patient's callers
  missed adds a patient junction that is public by construction. It can only
  turn private → shared_junctions, and cannot turn private → public or make a
  missed public event visible.
- **The false-private error mode comes from the panel side.** A patient
  junction becomes falsely private when panel genomes' callers miss it. The fix
  is to genotype the patient's junctions in the panel's evidence. That needs
  the panel reads, or a Sniffles2 mode that force-calls a VCF against `.snf`
  files, and Sniffles2 2.8.1's `--genotype-vcf` takes alignments only. The
  restricted environment has no panel reads.
- **Cost.** It needs one more genome-wide Sniffles2 pass per patient plus a
  BND/INV writer for Sniffles2's input format. That is about the runtime of
  the single-sample Sniffles2 call, for a column that is informational only.

What fills the gap: Sniffles2 population mode re-genotypes the panel from its
own `.snf` evidence, and `panel_frequency_comparison_summary.tsv` measures how
often per-sample calling misses a carrier.
