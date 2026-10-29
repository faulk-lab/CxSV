# External SV catalogues as additional background (roadmap item 4)

Status (2026-10-28): **acquisition done.** All four candidate catalogues are
open, licence-compatible with clinical research, downloaded with checksums and
converted to the shared junction-site schema. **Matching integration (scorer
reports which catalogue contains each call) is not done yet.**

Code: [`acquisition/catalogues/`](../acquisition/catalogues/README.md), the
conversion rules are listed there. Data: `/mnt/sargasso1/CxSV/datasets/external_catalogues/`
(the README there gives the layout). Every number below was measured from the
downloaded files, except the entries marked *(source)*.

## Summary

| | CoLoRSdb | HGSVC3 | HPRC release 2 | gnomAD-SV |
|---|---|---|---|---|
| Version used | v1.2.0 | HGSVC2024v1.0 (Variant_Calls 1.0) | MC graph v2.1 (`hprc-v2.1-mc-chm13`) | v4.1 |
| Release date | 2025-02-05 (Zenodo) | files dated 2024-09-18 (FTP) | files dated 2026-07-29 (S3) | VCF dated 2024-04-17 (GCS); v4 SVs announced 2023-11 |
| Licence / terms | CC-BY-4.0 [1] | HGSVC3 data reuse statement: open, Fort Lauderdale principles [4] | Public domain; HPRC data-use best practices [7] | CC0 1.0 [10] |
| Builds offered | GRCh38, CHM13 [1] | GRCh38-NoALT, T2T-CHM13v2.0 [3] | GRCh38, CHM13v2.0, GRCh37 graphs [6] | GRCh38 only [9] |
| Build used | CHM13 (native) | T2T-CHM13v2.0 (native) | CHM13v2.0 (native) | GRCh38 → hs1 liftover |
| Individuals | 1,427 (max NS in the CHM13 SV VCF) | 65 (130 haplotypes) | 231 diploid (232 VCF samples − GRCh38 path) | 63,046 *(source)* [11] |
| Platform / method | PacBio HiFi reads; pbsv + Jasmine [1] | HiFi + ONT UL Verkko phased assemblies; PAV 2.4 + SV-Pop [3,5] | HiFi + ONT UL phased assemblies; Minigraph-Cactus graph, vcfwave-decomposed [6] | Illumina short-read WGS; GATK-SV [11,12] |
| Files used (bytes) | SV VCF 61,611,516 | insdel_sym 6,715,391; inv_sym 16,191; complex-SV BED 136,243 | wave VCF 2,480,200,481 | sites VCF 1,743,544,245; chain 2,678,042 |
| Published checksum | md5 (Zenodo) ✓ | md5 (MANIFEST) ✓ | none for the VCF (multipart ETag): sha256 computed; tbi md5 ✓ | md5 (downloads page = GCS hash) ✓; chain md5 (UCSC md5sum.txt) ✓ |
| AC / AN / AF | INFO AC, AN, AF, NS, AC_Hom/Het/Hemi, nhomalt | none in INFO; computed from phased GT | INFO AC/AN/AF (include GRCh38); recomputed without GRCh38 | INFO AC, AN, AF + genotype counts |
| Per-population AF | no | no (1kGP labels would allow it) | no | yes: afr ami amr asj eas fin mid nfe rmi sas (+ XX/XY) |
| SV types in source | DEL, INS, INV (≥20 bp) | INS, DEL, INV; 1,909 complex (CPX) in a BED | sequence alleles of all sizes (SNV → SV), INV flag | DEL, DUP, INS (incl. MEI), INV, CPX, CTX, BND, CNV |
| INS encoding | sequence-resolved REF/ALT | symbolic `<INS>` + SVLEN | sequence-resolved | symbolic `<INS>`/`<INS:ME:*>` + SVLEN |
| BND / translocation | none | none (complex SVs in BED, STRUCT_REF/QRY) | none | `<BND>`/`<CTX>` with CHR2, END2; CPX_INTERVALS |
| **Junction sites written** | **456,811** | **192,175** | **1,170,364** | **1,201,011** |
| DEL / DUP / INV / INS / BND | 266,191 / 0 / 854 / 189,766 / 0 | 83,016 / 203 / 375 / 106,673 / 1,908 | 386,875 / 0 / 106 / 783,383 / 0 | 619,051 / 261,392 / 12,899 / 294,521 / 13,148 |
| Liftover site loss | n/a | n/a | n/a | 2.13 % (26,101 of 1,227,167 sites; table below) |

Sizes on the NAS: CoLoRSdb 66 MB, HGSVC3 9.7 MB, HPRC 2.4 GB, gnomAD-SV 1.7 GB
(about 4.3 GB downloaded, far below the 50 GB confirmation threshold; 21 TB
free on `/mnt/sargasso1` before and after).

## Per catalogue

### CoLoRSdb v1.2.0
* Zenodo record 14814308 (concept DOI 10.5281/zenodo.11511512), versions v1.0.0
  (2024-06-10), v1.1.0 (2024-07-31), v1.2.0 (2025-02-05, latest) [1,2]. The README
  says the `pbsv.jasmine` VCF holds insertions, deletions and inversions, the
  inversions as symbolic `<INV>` [1]. v1.1.0 collapsed duplicate SVs that came from double Jasmine merging [1].
* Measured: 839,714 records, no FILTER values, max NS = 1,427 individuals.
  382,903 records are < 50 bp (pbsv reports from 20 bp) and are not converted.
* Aggregates HiFi genomes from several institutions (largest source: Children's Mercy
  Research Institute) [2], including HPRC genomes [15]. The v1.0.0 merge drew on
  1,381 samples for autosomes [13]. No per-population AFs are published.

### HGSVC3 (HGSVC2024v1.0)
* 65 individuals, mostly from 1kGP, 130 haplotype assemblies (Logsdon et al.
  2025 [5]). Variant calls by PAV 2.4.0.1 and merged with SV-Pop, on both references
  "using identical data, tools, and versions" [3].
* Data reuse statement: open, Fort Lauderdale principles [4]. The first global
  paper is published [5], so the producers' first-publication reservation is met.
* Measured: 188,500 SV VCF records (DEL 81,549, INS 106,675, INV 276) plus 1,909 complex SVs.
  1,478 complex SVs tile their span exactly with STRUCT_REF segments and are also
  decomposed into DEL/DUP/INV sites. The other 431 are written as the outer-span BND only.
* Of the 65 individuals, 19 are in our panel registry (8 in 1kGP-LRSC, 16 in ONT-1019, 1 in PacBio).

### HPRC release 2 (Minigraph-Cactus v2.1, CHM13)
* 232 assembled samples; HG00272 held out of the graph; GRCh38 and CHM13
  included as reference paths [6]. The resources README warns that the
  pangenomes "have not been fully QC'd, are not published, and may have known issues" [6].
  Preprint: HPRC2 [8].
* No SV-only call set is published; the SV content is the graph-decomposed VCF.
  We use the vcfwave-normalised VCF (`.wave.vcf.gz`, the most processed of the three) [6].
* Measured: GRCh38 removed before counting. 1,536,848 split alleles with ≥ 50 bp
  net change or the INV flag. 132,600 are carried only by GRCh38 (AC = 0 afterwards) and are dropped.
  232,570 same-record alleles are merged into a shared junction site (AC summed).
* HPRC data are public domain (not patented, not subject to copyright); the
  data-use page adds recommended best practices for users [7].
* 50 of the 231 individuals are in our panel registry (42 in ONT-1019); 4 are also in HGSVC3.

### gnomAD-SV v4.1
* 63,046 unrelated genomes; 1,199,117 high-quality SV sites of types DEL, DUP,
  INS, INV, mCNV, CTX, CPX [11]. GRCh38-native; no CHM13 release [9].
* Terms: free download, no publication restrictions; primary data CC0 [10].
* Measured: 2,154,486 records; 1,199,117 PASS (equal to the published count), and only
  PASS records are converted. 711 PASS multiallelic CNVs are not representable and are dropped.
  BND records are all non-PASS (UNRESOLVED), so the 13,148 BND sites are CPX outer spans and 99 CTX.
* **Liftover** with UCSC `hg38ToHs1.over.chain.gz` (md5 578ede17f56d43f595f23d174acfe592),
  each breakpoint lifted independently as a 1-bp interval:

| svtype | sites in | sites out | lost: unmapped end | lost: ends on two contigs | loss % | span changed > 50 % (kept) |
|---|---|---|---|---|---|---|
| DEL | 635,602 | 619,094 | 16,123 | 385 | 2.60 | 6,354 |
| DUP | 269,207 | 261,398 | 7,412 | 397 | 2.90 | 6,346 |
| INV | 12,966 | 12,900 | 63 | 3 | 0.51 | 33 |
| INS | 296,184 | 294,525 | 1,659 | 0 | 0.56 | 0 |
| BND | 13,208 | 13,149 | 59 | 0 | 0.45 | 0 |
| all | 1,227,167 | 1,201,066 | 25,316 | 785 | 2.13 | 12,733 |

  (Sites in = PASS sites after CPX decomposition. 55 exact duplicate rows are then collapsed,
  which leaves 1,201,011.)

## Findings that matter for integration

1. **No catalogue is excluded.** All four are open, need no access approval,
   and have licences compatible with clinical-research use (CC-BY-4.0
   requires attribution; HGSVC3 asks users to follow the Fort Lauderdale principles).
2. **The catalogues overlap the panel and each other.** HGSVC3 and HPRC are
   1kGP individuals (19 and 50 in our registry); CoLoRSdb contains HPRC samples.
   Frequencies must not be added across sources as if they were independent.
   Report per-catalogue hits instead.
3. **Only gnomAD-SV and the HGSVC3 complex-SV BED carry complex or
   translocation junctions** (BND rows). The long-read catalogues are almost
   entirely DEL/INS, so for CxSV junctions they mostly confirm simple components.
4. **Platform.** CoLoRSdb is HiFi-read calling; HGSVC3/HPRC are assembly-based.
   None is ONT read-based like the panel. gnomAD-SV is short-read: insensitive
   in repeats, and breakpoints are imprecise.
5. **HPRC v2.1 is not fully QC'd** by its producers' own statement [6]. Treat a
   hit there as supporting evidence and do not let it override other evidence.
6. **Liftover changed the span by more than 50 %** for 12,733 kept gnomAD DEL/DUP/INV
   sites (1.1 %). These were kept and flagged in the stats only; the matcher may want
   to treat them with a size-agnostic breakpoint window.

## References
1. CoLoRSdb v1.2.0, Zenodo record 14814308 and its README.md. https://zenodo.org/records/14814308 (doi:10.5281/zenodo.14814308)
2. CoLoRSdb workflows, https://github.com/juniper-lake/CoLoRSdb ; https://colorsdb.org/
3. HGSVC3 Variant_Calls 1.0 README and MANIFEST. https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/release/Variant_Calls/1.0/
4. HGSVC3 data reuse statement. https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC3/README_HGSVC3_datareuse_statement.md
5. Logsdon GA, Ebert P, Audano PA et al. Complex genetic variation in nearly complete human genomes. Nature 644, 430–441 (2025). doi:10.1038/s41586-025-09140-6
6. HPRC pangenome resources (release 2) README, https://github.com/human-pangenomics/hpp_pangenome_resources ; MC v2.1 README, https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/release2/minigraph-cactus/v2.1/README
7. HPRC Data Use: best practices (14 April 2025). https://humanpangenome.org/data-use/
8. HPRC2: A human pangenome reference with near-complete coverage of common genetic variation (preprint). doi:10.64898/2026.07.21.739710
9. gnomAD downloads (v4 structural variants). https://gnomad.broadinstitute.org/downloads
10. gnomAD policies / terms of use. https://gnomad.broadinstitute.org/policies
11. Structural variants in gnomAD v4 (gnomAD blog, Nov 2023). https://gnomad.broadinstitute.org/news/2023-11-v4-structural-variants
12. Collins RL et al. A structural variation reference for medical and population genetics. Nature 581, 444–451 (2020). doi:10.1038/s41586-020-2287-8
13. HiFi Solves CoLoRS collection overview. https://hifisolves.org/collections/consortium-of-long-read-sequencing-colors/overview
14. UCSC liftOver chain hg38ToHs1.over.chain.gz. https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/
15. Genome Research 35:559 (2025), "The additional diagnostic yield of long-read sequencing in undiagnosed rare diseases" (describes CoLoRSdb as >1,400 PacBio genomes including HPRC). https://genome.cshlp.org/content/35/4/559
