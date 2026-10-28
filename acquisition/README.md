# Background-panel acquisition

Code that built the CxSV **background panel**: public, presumed-healthy long-read
human genomes, aligned to T2T-CHM13v2 (`hs1`), whose CxSV calls form the
reference that patient CxSVs are filtered against (`role: reference` cohorts in
`config/config.yaml`; see `config/faulk_lab_background.yaml`).

These scripts are operational code for the Faulk lab fleet. Paths default to the
lab NAS (`/mnt/sargasso1/CxSV`), and node names (freya, shadowfax, athena, betsy,
princess, keurig) appear in `orchestrate_alignment.py`. The design notes in each
file's docstring say why each safeguard exists; most were added after a
specific failure.

## Data layout (lab NAS)

```
/mnt/sargasso1/CxSV/
├── code/                      ← deployed copies of the scripts below
├── refs/                      ← hs1.fa, .fai, hs1.map-ont.k17.mmi, ref_cache/ (pipeline)
├── resources/                 ← annotation BEDs used by the pipeline
└── datasets/
    ├── 1kGP_LRSC_ONT_CHM13/   ← 467 ONT BAMs: held 100 + R10 PMDV (A1) + R9 MINIMAP_2026 (A2)
    ├── 1kGP_ONT_1019_CHM13/   ← 1,017 CRAMs aligned in-house from PRJEB89727
    ├── 1kGP_PacBio_CHM13/     ← 135 PacBio HiFi BAMs (A3)
    ├── manifests/             ← what to fetch, with expected sizes / md5s
    ├── logs/                  ← download ledgers, per-sample integrity records, claims
    ├── quarantine/            ← outputs from faulty hardware, kept for comparison only
    └── registry/              ← sample_registry.tsv, individuals.tsv, sample_sheets/
```

## Workflow and scripts

| Step | Script | What it does |
|---|---|---|
| 1. Manifest | `build_ena_manifest.py` | ENA study → one row per file, with publisher md5 and a raw-yield coverage estimate (`base_count / 3.1 Gb`, recorded as `ena_base_count`) |
| | `probe_odp.py`, `make_odp_manifest.py` | One-off: probe the AWS SRA Open Data mirror for each run and build `manifest_PRJEB89727_odp.tsv` (the mirror served 52 MB/s per stream vs 0.65 MB/s from ENA) |
| 2. Download | `download_manifest_v2.py` | Resumable, parallel, verified download of any manifest (size, md5 where published, `samtools quickcheck` for BAMs); refuses to start below a free-space floor (`--min-free-tb`). Replaces the original `download_manifest.py` |
| | `relocate_held_bams.sh` | Two-pass copy-then-delete move of the 100 held 1kGP-LRSC BAMs from sargasso2; the delete pass acts only on ledger-verified copies |
| 3. Align | `align_sra_to_chm13.sh` | One `.sra` → sorted CRAM: `fasterq-dump \| minimap2 2.24-r1122 -x map-ont -a --eqx -K 10g \| samtools sort`, matching the recipe of the held BAMs; refuses any other minimap2 version |
| | `orchestrate_alignment.py` | Fleet dispatcher: atomic per-sample claim files, per-node memory/thread settings, busy checks, stale-claim reaper, drain files (`logs/drain.<node>`), node-error handling, per-sample integrity records (read count of `.sra` == primary alignments in CRAM) |
| | `rebuild_integrity.py` | Rebuilds `logs/align/integrity/<ACC>.tsv` from ground truth after the shared ledger was corrupted by concurrent SMB appends |
| 4. Clean up | `delete_verified_sra.py` | Deletes a `.sra` only when its CRAM has a matching integrity record, passes quickcheck and still has the verified size |
| 5. Register | `registry_schema.py` | Registry columns, status/coverage-source vocabularies, GM→NA identifier normalisation |
| | `build_registry.py` | Rebuilds `registry/sample_registry.tsv`, `individuals.tsv` and the per-cohort pipeline `sample_sheets/` from manifests, ledgers, integrity records and the files on disk. Idempotent; re-run after any acquisition step |
| Ops | `ops/betsy_stress_test.sh` | Node acceptance test: production minimap2 loop with output hashing against a reference node, per-core CPU and DRAM checks |

## Decisions recorded in the code

- **Reference:** T2T-CHM13v2 as UCSC `hs1` (chr1..chrM). The RefSeq copy
  (`NC_*` names) is not compatible with the held BAMs.
- **Coverage floor:** 17x. For the in-house CRAMs it applies to raw yield
  (`ena_base_count`), whose measured aligned:raw ratio is 0.960, so 17x raw ≈ 16.3x
  aligned. For the LRSC and PacBio BAMs it applies to measured depth (`cramino`,
  `measured_chr20_depth`). Files below the floor are kept and registered, and
  excluded from analysis with the pipeline's `min_coverage`.
- **One release per individual in PacBio:** the W release is preferred where both
  W and C exist.
- **Skipped:** the two PRJEB89727 runs available only as ENA FASTQ (HG02013,
  HG00377); recorded in the registry as `excluded`.
- **Storage:** alignments are kept as CRAM; `.sra` files are deleted once their
  CRAM is verified.
- **Never append to one shared file from several hosts on SMB.** O_APPEND is not
  atomic across clients, so every per-sample record is written as its own file
  via temp + rename.

## Re-running

```bash
python3 acquisition/build_registry.py --dry-run     # counts only
python3 acquisition/build_registry.py               # rewrite registry + sample sheets
bash scripts/run_pipeline.sh --cohort 1kGP_ONT_1019_CHM13 \
     --extra-config config/faulk_lab_background.yaml
```
