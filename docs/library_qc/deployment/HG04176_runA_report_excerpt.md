# CxSV score report: HG04176

Explicit points system modelled on ACMG/ClinGen CNV scoring plus computed chromothripsis criteria 1, 2, 5, 7 and 8 (docs/SCORING.md). **Every threshold and weight is an operational, cohort-specific choice, not a validated cutoff.** This report is research output for review by qualified clinical scientists; it is not a diagnosis.

- loci scored: 1256; sex male (inferred, chrX/autosome depth 0.469)
- panel release: cxsv-panel-20261105-c31ecc9f8fa2; matching tolerance 500 bp, min junction fraction 0.75
- resources: repeats.final.bed (RepeatMasker_CHM13v2 4.1.2p1.2022Apr14, sha256 fb94ad9143cf); scoring/clingen_dosage.tsv (ClinGen_dosage 2026-10-28, sha256 c52a6c00e3b7); scoring/gnomad_constraint.tsv (gnomAD_constraint v4.1, sha256 0a2bcf947072); scoring/clingen_region_dosage.tsv (ClinGen_regions 2026-04-28, sha256 c57073a02580); scoring/clingen_region_dosage.liftover.tsv (ClinGen_regions 2026-04-28, sha256 9036d51d396e); scoring/hpo_genes_to_phenotype.txt (HPO_genes_to_phenotype v2026-09-01, sha256 507a17bff9c4); scoring/hp.obo (HPO_ontology v2026-09-01, sha256 93dace952fcb); scoring/par.bed …
- external catalogues: CoLoRSdb v1.2.0 (VERIFIED); HGSVC3 HGSVC2024v1.0 (VERIFIED); HPRC release2-mc-v2.1 (VERIFIED); gnomAD-SV v4.1 (VERIFIED)
- library QC (docs/LIBRARY_QC.md): **FLAG** (group ONT_1kGP_1019; svim_only_any_per_gb 122 (z 9.4 > 6); svim_only_bnd_per_gb 90.7 (z 6.5 > 6)); reads with a supplementary alignment 0.158034, to another chromosome 0.143308; SVIM-only BND/TRA per Gb 90.7272 (z 6.45431), SVIM-only records per Gb 122.007 (z 9.3759)
  - **This library is an outlier for single-caller breakend calls; treat Tier-1 loci supported by one caller only with caution.**
  - Conservative mode is on and demoted no locus: all 0 Tier-1 loci have multi-caller or copy-number support.
- tiers: Tier1_high_priority 0, Tier2 4, Tier3 24, Benign_likely 1228

## Ranking

| Rank | Locus | Class | Panel | Genes | Total | Tier | Top contributions |
|---:|---|---|---|---|---:|---|---|
| 1 | chr6:171735629-171881754 | Multi_breakpoint_cluster | public (12/1582) | DLL1,FAM120B | 0.57 | Tier2 | panel_privacy -1.00; gene_dosage +0.90; ct_breakpoint_tier +0.25 |
| 2 | chr9:146681626-146874916 | Multi_breakpoint_cluster | public (3/1582) | COL5A1 | 0.53 | Tier2 | panel_privacy -1.00; gene_dosage +0.90; ct_breakpoint_tier +0.25 |
| 3 | chr20:65251868-65304832 | Multi_breakpoint_cluster | public (198/1582) | EEF1A2,KCNQ2 | 0.53 | Tier2 | panel_privacy -1.00; gene_dosage +0.90; ct_breakpoint_tier +0.25 |
| 4 | chr11:1437321-1509740 | Multi_breakpoint_cluster | public (143/1582) | BRSK2 | 0.51 | Tier2 | panel_privacy -1.00; gene_dosage +0.90; ct_breakpoint_tier +0.25 |
| 5 | chr1:7034334-7079922 | Complex_other | public (2/1582) | CAMTA1 | 0.46 | Tier3 | panel_privacy -1.00; gene_dosage +0.90; ct_cn_oscillation +0.30 |
| 6 | chr2:1788018-1856196 | Complex_other | public (966/1582) | MYT1L | 0.35 | Tier3 | panel_privacy -1.00; gene_dosage +0.90; exons +0.15 |
| 7 | chr11:55535303-55628561 | Chromothripsis_like | public (476/1582) | OR4C11,OR4C6,OR4P4,OR4S2 | 0.30 | Tier3 | panel_privacy -1.00; complexity_class +0.30; ct_cn_oscillation +0.30 |
| 8 | chr19:11551563-11552946 | Complex_other | public (1241/1582) | ODAD3,PRKCSH | 0.20 | Tier3 | panel_privacy -1.00; gene_dosage +0.90; caller_agreement +0.10 |
| 9 | chr1:196039062-196230474 | Multi_breakpoint_cluster | public (365/1582) | CFH,CFHR2,CFHR4 | 0.19 | Tier3 | panel_privacy -1.00; ct_cn_oscillation +0.30; ct_breakpoint_tier +0.25 |
| 10 | chr2:242066471-242106830 | Multi_breakpoint_cluster | public (374/1582) | BOK,THAP4 | 0.15 | Tier3 | panel_privacy -1.00; ct_cn_oscillation +0.30; ct_breakpoint_tier +0.25 |

[excerpt: report header, library QC lines, demoted-loci section and the first ranking rows; full report on freya]
