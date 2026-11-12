# Patient intake (roadmap item 16k)

A patient CxSV is private only relative to a panel called the same way, and the
calling starts from the alignments. Every panel genome and the positive-control
host were aligned to T2T-CHM13v2 (UCSC **hs1** names) with one recipe; the
intake brings any patient input to an hs1 CRAM made with that recipe, so the
patient scorer (`workflow/patient/Snakefile`) compares like with like. Before
the intake, the nine monkey_meth test genomes (2026-11-10) were realigned with
ad-hoc scripts, and the reused BOC_3559 host CRAM carried the read-group sample
name `pos_control_3559` instead of the sheet's `BOC_3559`.

The intake is **not** part of the calling fingerprint (alignment is upstream of
calling). The recipe is fixed in code, not in config, so it cannot drift from the
panel by a site overlay.

## Recipe

```
samtools view -u -F 0x900 <input>            # BAM/CRAM: primary records only, each read once
  | samtools fastq -T MM,ML,MN,qs -          # FASTQ: intake_tool.py fastq-tags instead
  | minimap2 2.24-r1122 -y -ax map-ont --eqx -k 17 -K <batch> -R '@RG\tID:<sample>\tSM:<sample>\tPL:ONT\tDS:...' hs1.map-ont.k17.mmi -
  | samtools sort -O cram --reference hs1.fa  ->  <sample>.hs1.cram + .crai
```

| Element | Value | Why |
|---|---|---|
| aligner | minimap2 **2.24-r1122** | the panel's version; the script refuses any other (`exit 4`). It is in the container image (`envs/cxsv-pipeline.explicit.txt`, `minimap2-2.24-h7132678_1`); the intake needed no image change |
| preset, k | `-x map-ont`, `-k 17` (w = 10) | the panel index parameters; a prebuilt index is checked from its header (`MMI\2`, w, k, sequence count = `.fai` lines) |
| `--eqx`, `-a` | `=`/`X` CIGAR, SAM output | as the panel |
| `-y` + `-T MM,ML,MN,qs` | modified-base tags and dorado's qs copied into the alignments | only these four: `-T '*'` would also copy dorado's per-read `RG:Z`, which clashes with the new read group |
| `-K` | 4g (default) or 2g | memory and speed only, never the alignments |
| read group | `ID` = `SM` = sheet `sample_id`, `PL:ONT`, `DS` = sheet description + recipe | caller VCF sample columns then match the patient sheet |
| reference | the run's `reference` (bundle: `/bundle/reference/hs1.fa`) | CRAM encoded against hs1, whose MD5 cache ships in the bundle |

Inputs accepted, per sample (several inputs of one sample, e.g. flow cells, are
concatenated into one alignment):

| Input | Handling |
|---|---|
| unaligned BAM / dorado modBAM (`.bam`) | `samtools fastq -T MM,ML,MN,qs` |
| FASTQ (`.fastq`, `.fq`, optionally `.gz`) | comments reduced to valid `MM`/`ML`/`MN`/`qs` SAM tags (`intake_tool.py fastq-tags`): `minimap2 -y` copies the comment verbatim, and space-separated or non-SAM tokens (`runid=...`, `ch=...`) would make invalid SAM. `samtools import -T` cannot be used: it only splits tab-separated comments |
| BAM/CRAM from another aligner or reference naming | primary records only (`-F 0x900`), reads restored to their original orientation by `samtools fastq`, so supplementary/secondary records never duplicate a read; a CRAM made against another reference needs its FASTA (`input_reference`), otherwise the preflight stops with `cannot decode ... needs input_reference` |

Platform: ONT only (`map-ont`). The panel's PacBio genomes arrived aligned; a
PacBio intake would need the panel's PacBio recipe and is refused by the sheet
check.

## Checks (every one stops the sample; the CRAM is removed)

1. minimap2 is 2.24-r1122; the index is k = 17, w = 10, not homopolymer-compressed,
   with one sequence per `.fai` line.
2. Preflight: the first 2,000 primary records of every BAM/CRAM input decode.
3. **Read count:** input primary reads (sum over inputs, from `samtools fastq` or
   the FASTQ normaliser) == output primary mapped + unmapped (`samtools flagstat`).
4. Header: exactly one `@RG`, with `SM` = sample_id; exactly one minimap2 `@PG`,
   version 2.24-r1122, with `-x map-ont`, `-a`, `-y`, `--eqx`, `-k 17`.
5. Modified-base tags: if the first 2,000 input records carry `MM:Z`, the first
   2,000 output records must too.

The md5 of every input is computed from the same read pass and recorded in
`provenance.tsv` with its read count, the CRAM md5, minimap2's peak RSS and the
alignment wall time. `intake.json` holds the verdict, counts and recipe.

## Running

Intake sheet (tab-separated, header; `#` lines ignored):

| column | required | meaning |
|---|---|---|
| `sample_id` | yes | letters, digits, `.`, `_`, `-`; becomes `SM` and the patient sample id |
| `input` | yes | one path, or several separated by `,` |
| `individual_id` | no | default `sample_id` (carried into the patient sheet) |
| `platform` | no | `ONT` (default; the only accepted value) |
| `input_reference` | no | FASTA to decode CRAM inputs made against another reference |
| `description` | no | text for the read-group `DS` (basecaller, model, kit) |

Config (`config/intake.yaml` has the defaults):

```yaml
reference: /path/hs1.fa                      # the bundle sets /bundle/reference/hs1.fa
intake:
  sheet: /data/family_042/intake_sheet.tsv
  minimap2_index: null                       # or a prebuilt k17 map-ont hs1 index; null = build (9 GB)
  minibatch: "4g"                            # 2g on 60 GB nodes
  align: {threads: 26, sort_threads: 8, sort_mem: "2G", mem_mb: 72000}
```

Offline bundle (network off):

```bash
./run_cxsv.sh intake  --workdir /scratch/family_042 --config intake.yaml  --bind /data/family_042 --cores 32
# then, same work directory:
#   patient_runs: {family_042: {samples: /work/results/intake/patient_sheet.tsv}}
./run_cxsv.sh patient --workdir /scratch/family_042 --config patient.yaml --cores 32
```

Lab fleet: `snakemake -s workflow/intake/Snakefile --cores 32 --configfile <site>_intake.yaml`
with `intake.minimap2_index: /mnt/sargasso1/CxSV/refs/hs1.map-ont.k17.mmi`.
One sample can also be aligned with the script alone:
`REF_FA=hs1.fa REF_MMI=hs1.map-ont.k17.mmi bash workflow/scripts/intake_align.sh <sample> <out_dir> <input>...`.

Outputs, under `results/intake/`: `<sample>/<sample>.hs1.cram(.crai)`,
`intake.json`, `provenance.tsv`, `<sample>.flagstat.txt`, `<sample>.idxstats.tsv`,
`input_header.txt` (the inputs' headers without `@SQ`), `logs/`;
`patient_sheet.tsv` (sample_id, alignment, index, individual_id, platform; the
patient workflow's sheet format) and `intake_summary.tsv`.

## Memory and time

Measured on the nine monkey_meth genomes (ONT R10.4.1, 7.8-13.9 M reads, 26 or 20
minimap2 threads, 2026-11-10): minimap2 peak RSS 53-56 GB at `-K 4g`, 37-39 GB at
`-K 2g`; alignment 48-98 min. `samtools sort` adds `sort_threads x sort_mem`
(16 GB at the defaults), hence `mem_mb` 72000 (4g) or 54000 (2g). Building the
index takes a few minutes and about 20 GB.

## Tests

- `tests/unit/test_intake.py`: sheet parsing (all problems listed at once),
  FASTQ comment normalisation, index header check, flagstat parsing, read-count
  and header verdicts (wrong SM as in BOC_3559, minimap2 version, missing `-y`/
  `--eqx`/`-k 17`, a second minimap2 `@PG`), patient-sheet round trip, the recipe
  identical in the library and the script, minimap2 2.24 pinned in the image env.
- `tests/synthetic/check_intake.sh` (CI): four synthetic samples, a modBAM with
  per-read dorado RG and unmappable reads, a gzipped FASTQ with mixed comments,
  two BAMs from another reference naming with supplementary records, and a CRAM
  against another reference with `input_reference`. Every read must come out
  exactly once, with MM/ML identical per read, SM = sample_id, no dorado RG,
  hs1-style contigs; negatives: a k = 15 index, a CRAM without its reference, a
  PacBio row and a duplicate sample are refused.
