# Curated literature set of pathogenic germline complex SVs (roadmap item 11)

A breakpoint-level set of published germline complex rearrangements, lifted to
T2T-CHM13v2 (UCSC hs1 names), for testing the CxSV annotation, privacy and
scoring on known disease rearrangements, and for planting a few of them into
the positive control (detection). Full description, results and limitations:
[`docs/LITERATURE_SET.md`](../../docs/LITERATURE_SET.md).

| File | Role |
|---|---|
| `sources.tsv` | Every source: citation, DOI, PMID, file, URL, sha256, role (parsed / curated / context / not_obtained), how it was accessed |
| `curated/events.tsv` | One row per event: case, phenotype label, reported genes, inheritance, the authors' event type and classification |
| `curated/nazaryan2018_junction_types.tsv` | Junction type and sides for each of the 83 S2 Table junctions of Nazaryan-Petersen et al. 2018 |
| `curated/junctions_manual.tsv` | Junctions transcribed from Eisfeldt et al. 2024 (Table 2) and Bahrambeigi et al. 2019 (Table S2) |
| `curated/grochowski2024_table_s1.tsv` | DUP-TRP/INV-DUP segment coordinates of Grochowski et al. 2024 (Table S1) |
| `curated/cn_segments.tsv` | Copy-number changes reported for the curated events (used only by the scoring test) |
| `literature_set.py` | Junction model, stdlib .xlsx reader, derivative-chain → junctions, source parsers |
| `build_literature_set.py` | Parse + lift (item-4 `liftover.lift_points`) + write the versioned resource |
| `download_sources.py` | Re-fetch the source files and verify their sha256 |

## Sources

| source_id | Citation | Used |
|---|---|---|
| nazaryan2018 | Nazaryan-Petersen L, et al. Replicative and non-replicative mechanisms in the formation of clustered CNVs are indicated by whole genome characterization. *PLoS Genet* 2018;14:e1007780. doi:10.1371/journal.pgen.1007780 | S2 Table (83 junctions, 21 carriers), Table 1 (CMA, HGVS, phenotype, inheritance) |
| collins2017 | Collins RL, et al. Defining the diverse spectrum of inversions, complex structural variation, and chromothripsis in the morbid human genome. *Genome Biol* 2017;18:36. doi:10.1186/s13059-017-1158-6 | Additional file 3: derivative chromosomes of TL010, UTR22 (balanced germline chromothripsis) and TL009 (mosaic chromoanasynthesis) |
| eisfeldt2024 | Eisfeldt J, et al. A national long-read sequencing study on chromosomal rearrangements uncovers hidden complexities. *Genome Res* 2024;34:1774-1784. doi:10.1101/gr.279510.124 | Tables 1-2 (GRCh38 junctions; the Swedish GMS long-read pilot) |
| bahrambeigi2019 | Bahrambeigi V, et al. Distinct patterns of complex rearrangements and a mutational signature of microhomeology are frequently observed in PLP1 copy number gain structural variants. *Genome Med* 2019;11:80. doi:10.1186/s13073-019-0676-0 | Table S2: PLP1 complex gains (Pelizaeus-Merzbacher disease) |
| grochowski2024 | Grochowski CM, et al. Inverted triplications formed by iterative template switches generate structural variant diversity at genomic disorder loci. *Cell Genom* 2024;4:100590. doi:10.1016/j.xgen.2024.100590 | Table S1: 24 DUP-TRP/INV-DUP (23 MECP2, 1 PLP1) |
| clinvar_nstd102 | ClinVar structural variants distributed by dbVar as study nstd102 (GRCh38 variant calls, snapshot 2026-11-04) | Pathogenic / likely pathogenic complex substitutions and inversions (region level) |
| schuy2022 | Schuy J, Grochowski CM, Carvalho CMB, Lindstrand A. Complex genomic rearrangements: an underestimated cause of rare diseases. *Trends Genet* 2022;38:1134-1146 | context only (review; no breakpoint tables) |

Not obtained (named in the plan; see `sources.tsv` for why): Kloosterman et
al. 2011 *Hum Mol Genet* (germline chromothripsis), de Pagter et al. 2015
*Am J Hum Genet*, Kloosterman et al. 2012 *Cell Rep* Table S3, Gu et al. 2015
*Hum Mol Genet* (17p13.3), Carvalho et al. 2011 *Nat Genet* supplementary
tables. DECIPHER is not used (data-access agreement).

## Curation rules

1. **Facts only.** Coordinates, junction types, gene symbols, short phenotype
   labels and the authors' classification are recorded; no text is copied.
   Every hand-transcribed row names its table and row (`evidence`).
2. **Germline, published, ≥ 2 junctions.** An event is one rearrangement of one
   carrier with at least two novel junctions (or, for ClinVar, a complex
   substitution or inversion record). Simple translocations, single-junction
   events, rings and isodicentrics without a junction table are not included
   (Eisfeldt P1, P2, P7, P9, P10, P11). Ambiguous structures the authors
   describe as alternatives (Bahrambeigi BAB8924, BAB8940, BAB8955, BAB8960)
   and tandem duplications with only micro-insertions are not included.
3. **Junction = one novel adjacency with a side at each end** (`+` the
   sequence ending at the position is kept, `-` the sequence starting there):
   different chromosomes → BND, (+,-) → DEL, (-,+) → DUP, (+,+)/(-,-) → INV
   (`literature_set.svtype_of`), which is how Sniffles2/cuteSV/SVIM report it.
4. **Precision.** `bp`: sequence-resolved junction. `approximate`: array or
   mate-pair breakpoint, or a junction embedded in a long inverted repeat
   (DUP-TRP/INV-DUP junction 1). `region`: ClinVar interval ends.
5. **Nazaryan-Petersen sides.** S2 Table signs tell a direct join (equal signs)
   from an inverting one (opposite); deletion- vs duplication-like is taken from
   the case's Table 1 derivative or CMA segments, and head-to-head vs
   tail-to-tail from the derivative where it is given, otherwise by alternation
   (`basis` column). Three one-sided junctions of P1426_301 are not representable.
6. **Collins derivatives** are converted fragment by fragment
   (`derivative_junctions`); in UTR22 two telomeric fragments listed under the
   next derivative are re-attached (`fix_derivative_boundaries`, flagged).
7. **DUP-TRP/INV-DUP model** (Grochowski; Bahrambeigi BAB8938/BAB8964): the
   derivative is DUP1-TRP-DUP2 | TRP inverted | DUP1-TRP-DUP2, i.e. junction 1
   DUP2-end → TRP-end (head-to-head) and junction 2 TRP-start → DUP1-start
   (tail-to-tail), as for positive-control E01.
8. **Liftover** per breakpoint with the item-4 code: GRCh37 hg19 → hg38 → hs1
   (UCSC hg19ToHg38 + the item-4 hg38ToHs1 chain), GRCh38 → hs1. A junction is
   lost when an end does not lift or a DEL/DUP/INV lands on two contigs; a
   reversed chain block flips the side.
9. **Carrier sex / zygosity** are recorded only to synthesise read depth in the
   scoring test; when the paper does not give them, female (autosomal and
   unknown) or male (X-linked disease probands) is assumed and flagged.
