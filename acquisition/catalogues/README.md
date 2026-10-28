# External SV catalogues → junction-site tables

Roadmap item 4 (acquisition half). These scripts download published SV
frequency catalogues with checksums and convert each to one shared
**junction-site table**, so the matcher (a later step) can ask "is this
patient junction known in catalogue X, and how common is it?". Which
catalogues, why, and their licences: [`docs/EXTERNAL_CATALOGUES.md`](../../docs/EXTERNAL_CATALOGUES.md).

| Catalogue | Version | Native build | What is converted |
|---|---|---|---|
| CoLoRSdb | v1.2.0 | CHM13 | `CoLoRSdb.CHM13.v1.2.0.pbsv.jasmine.vcf.gz` |
| HGSVC3 | HGSVC2024v1.0 | T2T-CHM13v2.0 | `variants_T2T-CHM13_sv_insdel_sym`, `..._sv_inv_sym` VCFs and the PAV complex-SV BED |
| HPRC | release2-mc-v2.1 | CHM13v2.0 | `hprc-v2.1-mc-chm13.wave.vcf.gz` (Minigraph-Cactus graph VCF, vcfwave-normalised) |
| gnomAD-SV | v4.1 | GRCh38 → hs1 | `gnomad.v4.1.sv.sites.vcf.gz`, lifted with UCSC `hg38ToHs1.over.chain.gz` |

## Run it

```bash
conda env create -f envs/cxsv-catalogues.yaml      # once: bcftools/htslib 1.21, UCSC liftOver
conda activate cxsv-catalogues
bash acquisition/catalogues/run_all.sh \
    /mnt/sargasso1/CxSV/datasets/external_catalogues \
    /mnt/sargasso1/CxSV/refs/hs1.fa.fai  6
```

About 4.3 GB is downloaded. Converting all four takes under 10 minutes on
freya. `run_all.sh` prints the sizes and free space first. **Lab rule:** if a
future source pushes the total above about 50 GB, get the PI's explicit
go-ahead with those numbers before downloading.

| Script | Role |
|---|---|
| `sources.tsv` | Every file to fetch: URL, published checksum, and where that checksum is published |
| `download_catalogues.py` | Verified download to `<root>/<catalogue>/<version>/` and a `SOURCES.tsv` there (URL, UTC date, size, published checksum, computed md5 + sha256, `verified`). Exits non-zero on any mismatch |
| `junction_sites.py` | The schema: columns, contig order, orientation, validation, writer |
| `adapters.py` | Per-source record → sites logic (pure functions, unit-tested) |
| `liftover.py` | Breakpoint-wise liftover for GRCh38 sources and its loss bookkeeping |
| `convert_catalogue.py` | CLI: one catalogue → `<name>.<version>.junction_sites.tsv.gz` (+ `.tbi`), `.manifest.tsv`, and `.liftover_stats.tsv` if lifted |
| `make_summary.py` | All manifests → `<root>/external_catalogues_summary.tsv` |
| `run_all.sh` | All of the above, end to end |

Tests: `tests/unit/test_catalogue_junction_sites.py` (runs in CI; no data needed).

## The junction-site schema

Tab-separated, header row, 1-based, UCSC hs1 contig names:

```
chrom1  pos1  chrom2  pos2  svtype  svlen  ac  an  af  source  version  [af_<group> ...]
```

* `svtype` ∈ {DEL, DUP, INV, INS, BND}; translocations are BND.
* DEL/DUP/INV: `chrom2 = chrom1`, `pos2 = END` (VCF END, the last affected
  base); `pos1` is the VCF POS (the padding base, i.e. the last unaffected
  base on the left). INS: `pos2 = pos1` (the anchor base). BND: `(chrom2, pos2)` is the mate.
* Each site is written once, with `(chrom1,pos1) <= (chrom2,pos2)` in
  **hs1 .fai contig order** (by length: … chr7, chrX, chr9, chr8 …), then position.
  Rows are sorted the same way and the file is bgzipped and tabix-indexed
  (`tabix -s1 -b2 -e2 -S1`), so `tabix file.tsv.gz chr17:1-2000000` works.
* `ac`/`an`/`af` are `.` when the source does not give them.

Choices made here where the shared schema is silent (the matcher must know them):

1. **`svlen` is unsigned** (DEL svlen > 0). For BND it is the source's positive
   intra-chromosomal span when there is one (the outer span of a complex SV),
   otherwise `.`.
2. **`af`** is the source AF when given, else `ac/an`, else `.`.
3. **Minimum size 50 bp** (net length change) for every source, matching the
   pipeline. CoLoRSdb (pbsv reports from 20 bp) loses the most here.
4. **AC = 0 sites are dropped** (no carrier left; for HPRC these are alleles
   seen only in the GRCh38 haplotype).
5. **Complex SVs** (HGSVC3 PAV `CPX`, gnomAD `CPX`) become one **BND for the
   outer span** (chrom:start ↔ chrom:end, `svlen` = SVLEN) plus one typed site
   per component: gnomAD `CPX_INTERVALS` (DEL/DUP/INV; dDUP → DUP), HGSVC3
   `STRUCT_REF` segments (DEL; DUP/TRP/QUAD/MIX* → DUP; INV; INVDUP/INVTRP/INVQUAD
   → INV + DUP; NML skipped), the latter **only when the segment lengths tile
   the span exactly**. Component sites inherit the CPX's AC/AN/AF. Counts
   decomposed vs outer-span-only are in the manifest.
6. **Duplicates.** Alleles split from one multi-allelic record that give the
   same junction site are merged and their AC summed (HPRC graph VCF:
   different repeat sequences of equal length at one anchor). Across
   independent records only exact duplicate rows are collapsed.

Per-source specifics:

* **CoLoRSdb**: `bcftools norm -m -any`. pbsv sequence-resolved DEL/INS are
  classified from allele lengths; INV is symbolic with END. AC/AN/AF from INFO.
  No per-population AFs are published.
* **HGSVC3**: AC/AN computed with `bcftools +fill-tags` from the 65 phased
  genotypes (a `.` haplotype is not callable and is excluded from AN). Symbolic
  `<DEL>` has no END; END = POS + |SVLEN|.
* **HPRC**: the GRCh38 haplotype is a reference path in the graph, not an
  individual, so it is removed (`bcftools view -s ^GRCh38`, which recomputes
  AC/AN) before splitting alleles. Kept: alleles with |len(ALT) − len(REF)| ≥ 50,
  or vcfwave's INV flag (→ INV over the REF span when the net change is < 50 bp).
  Balanced complex alleles with < 50 bp net change and no INV flag are not
  represented. INFO/AF is not used (stale after removing GRCh38).
* **gnomAD-SV**: FILTER = PASS only (the 1,199,117 sites gnomAD calls high quality).
  CNV (multiallelic copy-number sites) cannot be represented and are dropped. CTX
  and BND → BND to (CHR2, END2). Per-ancestry AFs are kept as
  `af_AFR af_AMI af_AMR af_ASJ af_EAS af_FIN af_MID af_NFE af_RMI af_SAS`.
  **Liftover**: each distinct breakpoint is lifted alone as a 1-bp interval
  (UCSC liftOver, default `-minMatch=0.95`); INS lifts its single point. A site
  is lost if either end does not map, or if a DEL/DUP/INV's ends land on
  different contigs (BND ends may). Ends are re-oriented afterwards; `svlen`
  keeps the GRCh38 value. Loss by type: `gnomAD-SV.v4.1.liftover_stats.tsv`.
