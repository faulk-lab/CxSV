# CxSV score report: HG01879

Explicit points system modelled on ACMG/ClinGen CNV scoring plus computed chromothripsis criteria 1, 2, 5, 7 and 8 (docs/SCORING.md). **Every threshold and weight is an operational, cohort-specific choice, not a validated cutoff.** This report is research output for review by qualified clinical scientists; it is not a diagnosis.

- loci scored: 3891; sex male (inferred, chrX/autosome depth 0.51)
- panel release: cxsv-panel-20261105-c31ecc9f8fa2; matching tolerance 500 bp, min junction fraction 0.75
- resources: repeats.final.bed (RepeatMasker_CHM13v2 4.1.2p1.2022Apr14, sha256 fb94ad9143cf); scoring/clingen_dosage.tsv (ClinGen_dosage 2026-10-28, sha256 c52a6c00e3b7); scoring/gnomad_constraint.tsv (gnomAD_constraint v4.1, sha256 0a2bcf947072); scoring/clingen_region_dosage.tsv (ClinGen_regions 2026-04-28, sha256 c57073a02580); scoring/clingen_region_dosage.liftover.tsv (ClinGen_regions 2026-04-28, sha256 9036d51d396e); scoring/hpo_genes_to_phenotype.txt (HPO_genes_to_phenotype v2026-09-01, sha256 507a17bff9c4); scoring/hp.obo (HPO_ontology v2026-09-01, sha256 93dace952fcb); scoring/par.bed …
- external catalogues: CoLoRSdb v1.2.0 (VERIFIED); HGSVC3 HGSVC2024v1.0 (VERIFIED); HPRC release2-mc-v2.1 (VERIFIED); gnomAD-SV v4.1 (VERIFIED)
- library QC (docs/LIBRARY_QC.md): **FLAG** (group ONT_1kGP_1019; svim_only_any_per_gb 112 (z 9.1 > 6)); reads with a supplementary alignment 0.430335, to another chromosome 0.134029; SVIM-only BND/TRA per Gb 1.31554 (z 1.68635), SVIM-only records per Gb 112.175 (z 9.10274)
  - **This library is an outlier for single-caller breakend calls; treat Tier-1 loci supported by one caller only with caution.**
  - **Conservative mode DEMOTED 46 Tier-1 locus/loci** (a Tier-1 locus needs more than 0.5 of its junctions called by >= 2 callers or >= 1 copy-number-changed segment(s); otherwise it is reported as Tier2 with its points unchanged); 27 Tier-1 locus/loci kept. The demoted loci are listed under "Library QC: loci demoted by conservative mode".
- tiers: Tier1_high_priority 27, Tier2 386, Tier3 681, Benign_likely 2797


## Library QC: loci demoted by conservative mode

FLAGged library, conservative mode on: a Tier-1 locus needs more than 0.5 of its junctions called by >= 2 callers or >= 1 copy-number-changed segment(s); otherwise it is reported as Tier2 with its points unchanged.

| Rank | Locus | Class | Total | Tier before | Tier reported | Reason |
|---:|---|---|---:|---|---|---|
| 2 | chrX:78429952-78464314 | Multi_breakpoint_cluster | 1.44 | Tier1_high_priority | Tier2 | 0/4 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 3 | chr4:116373352-116476374 | DEL+INV | 1.43 | Tier1_high_priority | Tier2 | 2/5 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 7 | chr6:77008299-77033578 | DEL+INV | 1.36 | Tier1_high_priority | Tier2 | 0/2 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 8 | chr5:113261516-113307405 | Complex_inversion | 1.36 | Tier1_high_priority | Tier2 | 0/2 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 9 | chr13:34388470-34433450 | DEL+INV | 1.35 | Tier1_high_priority | Tier2 | 0/3 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 10 | chr2:50210121-50254411 | Complex_inversion | 1.35 | Tier1_high_priority | Tier2 | 0/2 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 11 | chr12:23314692-23403229 | DEL+INV | 1.35 | Tier1_high_priority | Tier2 | 0/4 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 12 | chrX:32205210-32258441 | DEL+INV | 1.35 | Tier1_high_priority | Tier2 | 1/4 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 13 | chr18:33823336-33873532 | DEL+INV | 1.35 | Tier1_high_priority | Tier2 | 0/3 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 14 | chr5:87791020-87862639 | Complex_inversion | 1.34 | Tier1_high_priority | Tier2 | 0/5 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 15 | chr2:49889165-49985208 | DEL+INV | 1.34 | Tier1_high_priority | Tier2 | 0/4 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 18 | chr4:116543311-116551386 | DEL+INV | 1.31 | Tier1_high_priority | Tier2 | 0/2 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 19 | chr3:171896900-171910509 | DEL+INV | 1.31 | Tier1_high_priority | Tier2 | 0/2 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 20 | chr9:3270910-3290954 | DEL+INV | 1.31 | Tier1_high_priority | Tier2 | 0/2 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 21 | chr4:76514237-76530473 | DEL+INV | 1.31 | Tier1_high_priority | Tier2 | 0/2 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 22 | chr5:37222147-37274350 | Complex_other | 1.30 | Tier1_high_priority | Tier2 | 0/3 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 24 | chr6:109093420-109112839 | DEL+INV | 1.30 | Tier1_high_priority | Tier2 | 1/2 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 25 | chr9:3465448-3508526 | DEL+INV | 1.30 | Tier1_high_priority | Tier2 | 1/2 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 27 | chr13:34268235-34324370 | DEL+INV | 1.29 | Tier1_high_priority | Tier2 | 0/6 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 32 | chrX:146215302-146245938 | DEL+INV | 1.25 | Tier1_high_priority | Tier2 | 1/2 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 33 | chr12:23461519-23477336 | DEL+INV | 1.24 | Tier1_high_priority | Tier2 | 0/2 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 34 | chr6:80235223-80267486 | DEL+INV | 1.24 | Tier1_high_priority | Tier2 | 0/2 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 35 | chr5:146826634-146830352 | DEL+INV | 1.21 | Tier1_high_priority | Tier2 | 0/2 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 36 | chr15:54727633-54768927 | Complex_inversion | 1.21 | Tier1_high_priority | Tier2 | 0/2 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 37 | chr14:31256262-31306268 | DEL+INV | 1.20 | Tier1_high_priority | Tier2 | 0/2 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 39 | chr5:81259308-81288161 | CN_plus_orientation_change | 1.20 | Tier1_high_priority | Tier2 | 1/2 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 41 | chr17:62575269-62585774 | DEL+INV | 1.18 | Tier1_high_priority | Tier2 | 1/2 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 42 | chr4:142634576-142681672 | DEL+INV | 1.18 | Tier1_high_priority | Tier2 | 1/2 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 43 | chr8:117940171-117979330 | DEL+INV | 1.18 | Tier1_high_priority | Tier2 | 1/2 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 44 | chr2:166453980-166454425 | DEL+INV | 1.17 | Tier1_high_priority | Tier2 | 0/2 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 45 | chr2:165816303-165817600 | DEL+INV | 1.17 | Tier1_high_priority | Tier2 | 0/2 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 46 | chr12:23838652-23840216 | DEL+INV | 1.17 | Tier1_high_priority | Tier2 | 0/2 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 47 | chr2:183408157-183415373 | DEL+INV | 1.16 | Tier1_high_priority | Tier2 | 0/2 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 48 | chrX:73230825-73246228 | Complex_inversion | 1.15 | Tier1_high_priority | Tier2 | 0/2 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 49 | chrX:109202145-109246622 | Complex_inversion | 1.15 | Tier1_high_priority | Tier2 | 0/2 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 50 | chrX:21079397-21104410 | Complex_inversion | 1.15 | Tier1_high_priority | Tier2 | 0/2 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 51 | chr16:56552605-56582012 | DEL+INV | 1.15 | Tier1_high_priority | Tier2 | 1/3 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 52 | chrX:31905010-31953134 | Complex_inversion | 1.14 | Tier1_high_priority | Tier2 | 0/3 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 54 | chrX:106957279-106958560 | DEL+INV | 1.12 | Tier1_high_priority | Tier2 | 0/2 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 56 | chr2:50371552-50415721 | DEL+INV | 1.10 | Tier1_high_priority | Tier2 | 1/2 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 57 | chr9:14101089-14173973 | DEL+INV | 1.10 | Tier1_high_priority | Tier2 | 0/3 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 58 | chr9:86675087-86712048 | DEL+INV | 1.10 | Tier1_high_priority | Tier2 | 1/2 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 61 | chr12:65173837-65190630 | DEL+INV | 1.05 | Tier1_high_priority | Tier2 | 0/2 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 62 | chrX:118769156-118809659 | DEL+INV | 1.05 | Tier1_high_priority | Tier2 | 1/3 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 64 | chr12:62330429-62410448 | Multi_breakpoint_cluster | 1.04 | Tier1_high_priority | Tier2 | 0/6 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
| 69 | chr12:116065987-116105141 | DEL+INV | 1.02 | Tier1_high_priority | Tier2 | 1/2 junctions with >= 2 callers, 0 copy-number-changed segment(s) |
## Ranking

| Rank | Locus | Class | Panel | Genes | Total | Tier | Top contributions |
|---:|---|---|---|---|---:|---|---|
| 1 | chr2:183467164-183469636 | DEL+INV | private (0/1561) | NCKAP1 | 1.46 | Tier1_high_priority | gene_dosage +0.90; panel_privacy +0.26; exons +0.15 |
| 2 | chrX:78429952-78464314 | Multi_breakpoint_cluster | private (0/1561) | TBX22 | 1.44 | Tier2 | gene_dosage +0.90; ct_breakpoint_tier +0.25; panel_privacy +0.19 |
| 3 | chr4:116373352-116476374 | DEL+INV | shared_junctions (0/1561) | ANK2 | 1.43 | Tier2 | gene_dosage +0.90; ct_breakpoint_tier +0.25; exons +0.15 |
| 4 | chr15:55022044-55065025 | DEL+INV | private (0/1561) | TCF12 | 1.41 | Tier1_high_priority | gene_dosage +0.90; panel_privacy +0.26; exons +0.15 |
| 5 | chr12:45765561-45804889 | Complex_other | private (0/1561) | ARID2 | 1.41 | Tier1_high_priority | gene_dosage +0.90; panel_privacy +0.26; exons +0.15 |
| 6 | chr9:9670474-9861375 | DEL+INV | private (0/1561) | PTPRD | 1.36 | Tier1_high_priority | gene_dosage +0.30; ct_cn_oscillation +0.30; ct_breakpoint_tier +0.25 |
| 7 | chr6:77008299-77033578 | DEL+INV | private (0/1561) | MYO6 | 1.36 | Tier2 | gene_dosage +0.90; panel_privacy +0.26; exons +0.15 |
| 8 | chr5:113261516-113307405 | Complex_inversion | private (0/1561) | APC | 1.36 | Tier2 | gene_dosage +0.90; panel_privacy +0.26; exons +0.15 |
| 9 | chr13:34388470-34433450 | DEL+INV | private (0/1561) | MAB21L1,NBEA | 1.35 | Tier2 | gene_dosage +0.90; panel_privacy +0.25; exons +0.15 |
| 10 | chr2:50210121-50254411 | Complex_inversion | private (0/1561) | NRXN1 | 1.35 | Tier2 | gene_dosage +0.90; panel_privacy +0.25; exons +0.15 |

[excerpt: report header, library QC lines, demoted-loci section and the first ranking rows; full report on freya]
