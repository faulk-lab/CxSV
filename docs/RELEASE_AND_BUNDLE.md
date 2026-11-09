# Panel releases and the offline bundle (roadmap item 2)

The patient scorer must run in a restricted environment, on a machine with no
network, from material whose origin can be checked. This document defines:

1. the **panel release**: what the panel builder hands over;
2. **callability**: how a release says where the panel could have seen a junction;
3. the **container image**: one pinned environment for both workflows;
4. the **offline bundle**: everything needed to score a patient with no network;
5. how to **build** and **run** each of them;
6. the **validation records** of the bundles. The current bundle is
   `cxsv-bundle-20261106-169503679a60` (commit 1695036, panel-v1.1 + library QC with the
   auto-detected group, conservative mode on; record in
   [library_qc/DEPLOYMENT.md](library_qc/DEPLOYMENT.md), which also covers the previous
   bundle `cxsv-bundle-20261105-95126e14f036`). The panel-v1.1 bundle
   `cxsv-bundle-20261105-8fc5bef725bd` (deleted 2026-11-06, PI-approved) is recorded in
   [panel_v1.1/README.md](panel_v1.1/README.md). Bundles from later commits need a release
   with their calling fingerprint (item 15 changed it).
   6a–6c below are the earlier bundles.

## 1. Panel release

The full lab panel (all registry files, no coverage floor) is built across the
fleet as described in [PANEL_BUILD.md](PANEL_BUILD.md). Its aggregate step is
the ordinary `workflow/panel` run, so the release below is the same.

A release is one directory. The patient workflow reads that directory and
nothing else (no reads, no per-sample VCFs).

| File | Content |
|---|---|
| `release.json` | identity and provenance (below) |
| `release_manifest.tsv` | `file  bytes  sha256` for **every** other file (relative paths, sorted) |
| `calling_fingerprint.json` | calling parameters, reference contigs (name, length, MD5), tool versions, calling-code digest: what a patient run must reproduce |
| `panel_junction_sites.tsv` | junction sites in the shared junction-site schema (`chrom1 pos1 chrom2 pos2 svtype svlen ac an af source version` + `af_<platform>` + `af_<superpopulation>`; docs/JUNCTIONS.md); `version` = the release id |
| `panel_junction_carriers.tsv` | the same sites with coordinate envelopes and the carrier individuals of each |
| `panel_junction_callability.tsv` | the same sites: `n_callable_samples`, `n_callable` (individuals), `n_individuals`, `callable_fraction`, `callable_<platform>`, next to `ac` |
| `panel_callability_summary.tsv` | sites callable in all / at least half / none of the individuals; callable bp per sample |
| `callability/` | `samples.tsv` + one CALLABLE mask per panel sample (`<cohort>__<sample>.callable.bed.gz`) |
| `panel_individuals.tsv` | sex, superpopulation, platforms, cohorts, samples per individual (the denominators) |
| `panel_cxsv_calls.tsv`, `panel_cxsv_loci.tsv`, `panel_summary.tsv` | CxSV loci per sample, merged loci with carriers by platform, counts per cohort/platform |
| `sniffles_population_sites.tsv`, `panel_frequency_comparison*.tsv` | Sniffles2 population-mode sites and their comparison with the clustered sites |
| `tool_versions.tsv` | versions of the fingerprinted tools |
| `cn/pon/<stratum>.npz`, `cn/pon/pon_members.tsv` | copy-number panel-of-normals per stratum (per-bin median, relative MAD, genome count of normalised 10 kb depth) and its member list; from panel-v1.2 (docs/CN_LINKAGE.md) |
| `cn/panel_cn_segments.tsv` | copy-number segments of every panel library (public filter of the copy-number track) |

**Release id:** `<release.name>-<YYYYMMDD>-<first 12 hex digits of the
calling fingerprint>`, e.g. `TEST-posctrl10slice-20261028-ad454c604bdd`. The
fingerprint part says which calling a patient must match. The name and date
tell apart releases that share a calling but differ in samples. A **test release**
(`release.test_release: true`) gets the suffix `-test`
(`cxsv-panel-20261106-c31ecc9f8fa2-test`), so it can never carry the id of a real
release built the same day from the same calling. Real releases keep the plain
form, so panel-v1.1 stays `cxsv-panel-20261105-c31ecc9f8fa2`. The id is fixed the
first time it is made (`results/panel/panel/release_id.txt`). Every score report
prints the id, `test release: yes/no` and the release's `release_manifest.tsv`
sha256.
`config/panel.yaml`:

```yaml
release:
  name: "cxsv-panel"       # letters, digits, . _ -
  date: null               # YYYYMMDD; null = the day the id is first made
  test_release: false      # true for panels that must never label a real patient
  description: ""
```

**`release.json`** (schema `cxsv-panel-release/1`) records:

- `release_id`, `name`, `date`, `test_release`, `description`, `created_utc`;
- `source`: the **git commit** the panel workflow ran from, with a dirty flag.
  It is read with git in a checkout, or from `workflow/COMMIT` in a
  `git archive` export (`.gitattributes` export-subst), which is how the code
  sits inside the bundle and the container;
- `calling_fingerprint`: sha256, code digest and parameters;
- `reference`: contig count, total length and MD5 per contig;
- `panel`: samples and individuals by cohort, **platform**, **superpopulation**
  and sex;
- `contents`: number of junction sites and CxSV loci;
- `junction_sites`: the clustering tolerance and which tables hold the
  carriers and frequencies;
- `callability`: the definition, its parameters and the site and per-sample summary;
- `tools`: fingerprinted tools plus mosdepth, snakemake and python;
- `runtime`: the container image tag, **image id**, **SIF sha256**, where
  those two came from (`ids_source`) and the environment-lock sha256 when
  built in the container; otherwise `CONDA_PREFIX`. The launchers pass the
  ids: the bundle's `run_cxsv.sh` from its `bundle.env`, the panel build's
  `scripts/panel_build/lib.sh` (`cxsv_container`) from the `bundle.env` next
  to the pinned SIF, or the SIF's own sha256 when it is not in a bundle (no
  image id then). A non-test release built in the container without either id
  is refused (`CXSV_ALLOW_MISSING_RUNTIME_IDS=1` overrides). panel-v1.0 was
  built before this (small fixes 2026-11): its `runtime.image_id` and
  `sif_sha256` are null, and its `release.json` stays as released. Its image
  is the SIF of `cxsv-bundle-20261029-1dcabd9770f7` (build `plan.json`; SIF
  sha256 `9f4e049fe5c5…`, image id `sha256:b6c1bded8b3a…` in that bundle's
  `bundle.env`, which is what `cxsv_container` now records for it);
- `files`: the file list.

Release tags are **not** created by the workflow. The commit hash is recorded
instead; tagging a release is the PI's decision. The first real release,
`cxsv-panel-20261102-cd7e1d194cd4`, carries the annotated tag `panel-v1.0`
(on its source commit da6cb90), written into `release.json` by
`write_release_metadata.py --finalise` (2026-11-03).

### The patient workflow refuses a release that does not verify

`workflow/lib/cxsv_release.py: verify` checks that every listed file exists
with the listed size and sha256, that no file is present that the manifest
does not list, and that `release.json` is listed. The check runs three times:

1. When the patient workflow is **parsed**: the run stops before any job
   starts, dry runs included.
2. As rule `verify_panel_release`, which writes the record
   `results/patient/panel_release_check.txt`: release id, test flag, source
   commit, fingerprint, manifest sha256, number of files verified, and the
   runtime image/SIF/bundle it ran in.
3. Inside `classify_patient_cxsv.py`, right before it reads the release. This
   catches a file changed during a run, even with its mtime kept.

The fingerprint check (`check_panel_compatibility`) follows as before.

## 2. Callability

Without callability, "absent from the panel" cannot be told apart from "the
panel could not have seen it there": no reads, or only multi-mapping reads.
Definitions are in `workflow/lib/cxsv_callability.py`; the rules are in
`workflow/rules/panel_callability.smk`.

- **Per panel sample:**
  `mosdepth -n -x --mapq <sv_params.min_mapq> --quantize 0:<min_depth>: -f hs1.fa`.
  This counts reads with MAPQ ≥ 20, the reads the SV callers use, and keeps
  the CALLABLE intervals (depth ≥ `min_depth`) as a gzipped BED. `gzip -n`
  makes the mask byte-reproducible.
- **CRAM 3.1** (`workflow/scripts/callable_mask.sh`): the pinned mosdepth
  0.3.10 cannot decode CRAM 3.1. It logs "Slice decode failure", exits 0 and
  leaves an empty mask. `panel_callability.cram31` chooses what happens:
  `transcode` (default) writes a temporary, record-count-checked CRAM 3.0 copy
  in the sample's `callability/` directory, uses it and deletes it (the input
  is never touched; this needs disk for one copy of the CRAM); `fail` stops
  with an error. In either case, decode errors in mosdepth's output or a mask
  with 0 callable bp fail the step. The panel-v1.0 build transcoded CRAM 3.1 at
  staging instead (`lib.sh cxsv_cram30`; 208 of 1,017 1kGP ONT CRAMs). The
  patient workflow has no mosdepth step: its CRAM readers (samtools 1.21,
  pysam 0.22.1 in the callers) decode CRAM 3.1. The synthetic patient run
  from a CRAM 3.1 copy gives the same CxSV loci, junctions, panel status and
  scores as from the BAM, and caller VCFs identical to a CRAM 3.0 copy (small
  fixes 2026-11). Between BAM and either CRAM, cuteSV's read support differs
  at 3 of 22 calls (e.g. RE 29 vs 28), reproducibly; the call set is the same.
- **Breakpoint callable:** a window of `flank` bp on at least **one** side of
  the breakpoint, `[pos-flank, pos)` or `[pos, pos+flank)`, lies inside one
  CALLABLE interval. One side is enough because the reads anchoring a junction
  come from either flank, and a homozygous deletion has no reads inside the
  deleted segment.
- **Junction callable in a sample:** both breakpoints callable.
- **Individual callable:** any of its samples is callable.
- **Locus callable:** callable at every junction of the locus.
- **Parameters** (`config/panel.yaml`, `panel_callability`): `min_depth: 10`
  (a heterozygous junction then has about 5 supporting reads in expectation,
  above `min_support` 3) and `flank: 100`. Not fingerprinted: they describe
  coverage, they do not decide which loci exist. They are recorded in
  `release.json`, and the patient side takes `flank` from there.
- **Panel side:** `panel_junction_callability.tsv` gives the number of panel
  samples and individuals callable at both breakpoints of every site.
- **Patient side:** the masks ship in the release, so callability can be asked
  at the patient's own junctions. This matters for private ones, which have no
  panel site by definition.
  - `cxsv_panel_status.tsv` gains `Panel_Callable` (individuals callable at
    every junction of the locus) and `Panel_Individuals`, next to
    `Panel_Carriers`.
  - `junction_panel_status.tsv` gains `panel_callable`.
  - A private locus with `Panel_Callable` close to `Panel_Individuals` is truly
    absent from the panel. A private locus with a low `Panel_Callable` was
    largely unobservable.
- A carrier can be below `min_depth` at a breakpoint. Carrier and callable
  counts are therefore reported side by side, never one capped by the other.
  In the test release, 15 of 902 sites are callable in no individual although
  they have carriers: low-MAPQ or low-depth calls.

Masks for a whole-genome panel are small, because the number of intervals is
set by the coverage dips, not by genome length. In the test release a 34 Mb
slice gives 0.2 to 5 kB per sample.

## 3. Container image

`container/Dockerfile`: ONE image for both workflows, for Docker, Singularity
CE and Apptainer.

- **Base:** `mambaorg/micromamba:2.10.0-debian12-slim`, pinned by digest.
- **Environment:** installed with
  `micromamba create -p /opt/cxsv --file envs/cxsv-pipeline.explicit.txt`.
  That is the explicit package list (URL + MD5 per package) of the freya env
  that built the panels. No solver runs, so the image holds exactly those
  package files.
- **Runtime settings:** no entrypoint, tools on `PATH`, `LC_ALL=C.UTF-8`, no
  bytecode writes. `CXSV_IMAGE` and `CXSV_ENV_LOCK_SHA256` are recorded by the
  workflows.
- **What it holds:** the environment only. The workflow code comes from the
  repository at a fixed commit, so the image is keyed by the lock
  (`cxsv-pipeline:<first 12 hex of sha256(Dockerfile, lock)>`), not by the commit.
- **Smoke test at build time:** every fingerprinted tool, snakemake, mosdepth
  and R/ggplot2.

The three env files must agree, and `tests/unit/test_env_lock.py` checks it:

- `envs/cxsv-pipeline.yaml`: the pins;
- `envs/cxsv-pipeline.lock.yaml`: the solved versions;
- `envs/cxsv-pipeline.explicit.txt`: the package files.

They had drifted: the committed lock predated a 2026-10-28 11:15
`mamba install minimap2=2.24 numpy` that re-solved the env (Python
3.11.18 → 3.11.9, R 4.4.3 → 4.4.0, about 60 libraries). It is now regenerated
from the env. The fingerprinted tools were the same in both.

mosdepth was added to the env at **0.3.10**: the newest build against htslib
1.21. 0.3.11+ need htslib ≥ 1.22, and htslib is part of the calling identity.

**Build:** `scripts/make_bundle.sh image --commit <rev>` does the following:

- builds from a `git archive` of the commit;
- `docker save` → `singularity build … docker-archive://…`;
- writes `image.env` with the image id, SIF sha256, Docker-archive sha256 and Dockerfile/lock sha256.

The output goes to `/mnt/sargasso1/CxSV/bundles/_images/img-<key>/` and is
reused while the Dockerfile and lock are unchanged.

| | |
|---|---|
| image | `cxsv-pipeline:ba820f7b2bd3` (key = sha256 of Dockerfile + lock), id `sha256:b6c1bded8b3a9357c63bbe9f44c0e0e60039a9d5fa92731094d93ebbdb2e15d0`, built from commit `144a4ff` (Dockerfile and lock unchanged since) |
| SIF | `cxsv-pipeline.sif`, 723,271,680 bytes, sha256 `9f4e049fe5c5eb89717e1f531d284836efc7c306e0103374872c8165c01ff967` |
| Docker archive | `cxsv-pipeline.docker.tar.gz`, 787,982,269 bytes, sha256 `e460187c70a6338bfb3b0362fade9c6f2dcb193f205c98be7c279191eeee8c8e` |
| env lock | `envs/cxsv-pipeline.explicit.txt` sha256 `f527840ef9206bc330005ca920bbf1e5e37d269b1a5d3bcfc7cfec5c579043ac` (291 packages) |
| base | `mambaorg/micromamba:2.10.0-debian12-slim@sha256:bd7f1077166c1c71e9375a0e54835f151d5f216b96301b75f888c53c46d8e78b` |
| build time | 89 s on freya (docker build 49 s with a warm package cache; docker save + SIF 40 s) |

The image reproduces the conda env's calling exactly: the test panel built in
the image (10 slices) has the same calling fingerprint (`ad454c604bdd`) as the
item-3 regression panel built in the conda env from the same slices. It also
has byte-identical `panel_junction_carriers.tsv`, `panel_cxsv_loci.tsv` and
`panel_cxsv_calls.tsv`, the same junction sites apart from the `version`
column, and identical merged-SV records (checked for HG00096 and HG00129).

### Running the workflows in the container

Both Snakefiles find their config files relative to themselves, and
`scripts_dir` resolves against the repository. The workflows therefore run
from a read-only copy of the code, with a separate writable working directory
(`--directory`). The exact command used to build the test release (panel
builder, Singularity, network off):

```bash
singularity exec --containall --net --network none --pwd /work \
    --bind <repo at commit>:/repo:ro --bind <panel workdir>:/work \
    --bind <crams>:<crams>:ro --bind /mnt/sargasso1/CxSV/refs:/mnt/sargasso1/CxSV/refs:ro \
    --bind /mnt/sargasso1/CxSV/resources:/mnt/sargasso1/CxSV/resources:ro \
    --bind /mnt/sargasso1/CxSV/datasets/registry:/mnt/sargasso1/CxSV/datasets/registry:ro \
    --env CXSV_IMAGE_ID=<image id> --env CXSV_SIF_SHA256=<sif sha256> \
    cxsv-pipeline.sif \
    snakemake -s /repo/workflow/panel/Snakefile --directory /work --cores 24 \
        --keep-going --printshellcmds --configfile /work/config/panel_site.yaml
```

Patient runs go through the bundle's `run_cxsv.sh` (section 5), which issues
the same `singularity exec` / `apptainer exec`, or
`docker run --rm --network none --user $(id -u):$(id -g) -v <bundle>:/bundle:ro -v <workdir>:/work …`.

## 4. Offline bundle

`/mnt/sargasso1/CxSV/bundles/<bundle id>/`, with bundle id
`cxsv-bundle-<YYYYMMDD>-<commit12>`:

| Path | Content |
|---|---|
| `run_cxsv.sh` | entry point (section 5) |
| `README.md` | run instructions, ids and checksums of this bundle |
| `bundle.json`, `bundle.env` | bundle id, commit, release id, image id, SIF sha256, Docker-archive sha256, base image, sizes by component, build times |
| `MANIFEST.sha256` | sha256 of every file (`sha256sum -c` format; covers `bundle.json`) |
| `bundle_manifest.tsv` | the same with sizes and components |
| `repo/` | the repository at the commit (`git archive`; `workflow/COMMIT` holds the hash) |
| `source/cxsv-<commit12>.tar.gz`, `source/cxsv-<commit12>.bundle` | the same commit as an archive and as a git bundle (`git clone cxsv-<commit12>.bundle`) |
| `container/` | `cxsv-pipeline.sif`, `cxsv-pipeline.docker.tar.gz`, `image.env`, `Dockerfile`, `cxsv-pipeline.explicit.txt` |
| `reference/` | `hs1.fa`, `hs1.fa.fai`, `ref_cache/` (htslib MD5 cache; mtimes kept, so Snakemake sees the cache as up to date) |
| `resources/` | `genes.final.bed`, `exons.filtered.bed`, `repeats.final.bed`, `segdup.final.bed`, `scoring/`; `literature/<version>/` = the curated literature set (item 11, annotation only), verified against its `MANIFEST.tsv` |
| `panel_release/<release id>/` | the panel release (verified before copying) |
| `catalogues/<name>/<version>/` | external catalogue junction-site tables (`.tsv.gz` + `.tbi`), each with its manifest, `SOURCES.tsv` and README; every table verified against its manifest's `output_sha256` before copying; `external_catalogues_summary.tsv` |
| `config/bundle_site.yaml` | site config with in-container paths (`/bundle/...`, `/work/...`) |

Nothing in the bundle refers to a path outside it. Inside the container the
bundle is `/bundle` (read-only), the run directory is `/work`, and the input
directories are mounted read-only at their own paths. The minimap2 index
(`hs1.map-ont.k17.mmi`, 9 GB) is not included: patients arrive aligned.

### Build

```bash
# 1. image (reused while container/Dockerfile and envs/cxsv-pipeline.explicit.txt are unchanged)
scripts/make_bundle.sh image  --commit <rev>
# 2. panel release, built in that image from the same commit (section 3 command)
# 3. bundle
scripts/make_bundle.sh bundle --commit <rev> --release <panel workdir>/results/panel/release
```

The builder takes code, Dockerfile and lock from `git archive <commit>`,
never from the worktree. For `git bundle` it creates a temporary local ref
and deletes it again; it pushes nothing and makes no tags.

**Literature set.** The bundle ships the newest *directory*
`<resources>/literature/lit-v*/` that holds a `MANIFEST.tsv`, or the set given
with `--literature DIR`; it is verified against its manifest and
`bundle_site.yaml` gets `literature_set.root`. If no set is found the build
stops before any copying, unless `--no-literature` is given (such a bundle
fails validation check 9). Selection: `scripts/lib/literature_select.sh`, tests
in `tests/unit/test_bundle_literature.py`.

*Bug fixed 2026-11-08:* the builder used `ls -d <resources>/literature/lit-v* |
sort -V | tail -1`, which on the lab NAS returned the file `lit-v1.0.build.log`
instead of the directory `lit-v1.0`, so the set was skipped without a message.
Bundles `cxsv-bundle-20261103-33c075a22071`, `cxsv-bundle-20261106-169503679a60`
(deployed) and `cxsv-bundle-20261108-c40e26528d5f` carry no literature set: their
patient runs have no `Literature_Match` / `Literature_Events` columns and no
`private/cxsv_literature_hits.tsv`. Scores and tiers are unaffected (the
literature match is annotation only).

**Validation** (`tests/lab/validate_bundle.sh`, 9 checks): integrity, release,
network off, offline Singularity run, offline Docker run, run outside the bundle
(with the lab copy of the same literature set), records identical, tampered
release refused, literature set shipped and verified (check 9, added with the fix
above; earlier records list 8 checks).

## 5. Running offline

```bash
B=/path/to/cxsv-bundle-...
$B/run_cxsv.sh verify-bundle            # every file against MANIFEST.sha256
$B/run_cxsv.sh verify-release           # the panel release against its manifest
$B/run_cxsv.sh netcheck                 # loopback only inside the container
$B/run_cxsv.sh patient --workdir /scratch/run1 --config patient.yaml \
    --bind /data/patient1 --cores 16 [--runtime singularity|apptainer|docker]
```

`patient.yaml` holds `patient_runs` (sample sheets with absolute paths under
the `--bind` directories) and `threads`. `run_cxsv.sh` mounts the bundle
read-only and isolates the container:

- **Singularity/Apptainer:** `--containall` (no `$HOME`, no host environment)
  and `--net --network none`.
- **Docker:** `--network none`, current uid. It loads the image from
  `container/cxsv-pipeline.docker.tar.gz` if needed and checks the image id
  against `bundle.env`.

`--no-netns` exists only for sites whose runtime forbids unprivileged network
namespaces. The bundle needs no network either way.

Inside the container, for every process (small fixes 2026-11):

- **Temporary files:** `TMPDIR`, `TMP` and `TEMP` are `/work/.tmp` (the run
  directory's `.tmp/`), and the same directory is mounted as `/tmp`. Python's
  `tempfile` and R's `tempdir()` follow `TMPDIR`. Without this, Singularity's
  `--containall` `/tmp` is a 64 MB tmpfs.
- **CRAM references:** `REF_PATH` and `REF_CACHE` are
  `/bundle/reference/ref_cache/%2s/%2s/%s`. htslib then finds every CRAM
  reference sequence by its `@SQ M5` in the bundle's MD5 cache, never through
  a header `UR` path or a download. SVIM 2.0.0 opens CRAMs without
  `--reference` and depends on this. The caller rules export the same values
  (`HTS_ENV`, from `ref_cache`); the panel build's `cxsv_container` sets them
  from `CXSV_REF_CACHE`.
- **SVIM output check** (`workflow/rules/caller_checks.smk`,
  `scripts/check_svim_output.py`): after the callers, `qc/svim_check.tsv`
  fails the sample when the SVIM log has a CRAM decode / truncation / reference
  error or lacks SVIM's closing `Done.`, or when SVIM has no calls on a contig
  where Sniffles2 has at least 20 (`svim_check.min_reference_calls`). This is
  the failure of the panel-v1.0 build when a chr7 cache file went missing
  (~269 samples; all re-called before the release). The check is a target of
  the panel, panel_build and patient workflows, and `score_cxsv` waits for it.
  It is outside the calling fingerprint, and no caller rule or output changed.

`tests/lab/check_container_env.sh BUNDLE WORKDIR` checks both settings under
Singularity and Docker.

**After every patient run (interim, panel-v1.2 bundles):** read `cn/stats.json`
(`sex_call`, `x_ratio`, `y_ratio`). If `sex_call` is `ambiguous`, flag the sample to the
analyst as a possible sex-chromosome aneuploidy or mosaicism, and treat the run's chrX
and chrY copy-number calls as unreliable. The patient report does not show this yet
(docs/CN_LINKAGE.md, section "After every patient run"; ROADMAP item 16 (c)).

**Open question for the restricted environment's IT** (unchanged): which
workflow engine and container runtime are allowed. The bundle assumes
Snakemake (inside the image) plus Singularity/Apptainer, and also ships a
Docker archive.

## 6. Validation record

### 6-v1.2. Bundle with panel-v1.2-r2, 2026-11-09

Record: [panel_v1.2/README.md](panel_v1.2/README.md) (bundle
`cxsv-bundle-20261109-84d380e08cc4`, release `cxsv-panel-20261108-fdac9ff80401`, 9 checks
including check 9 "literature set shipped").

**Known defect in earlier bundles (fixed in PR #37, main 761f75d).** `make_bundle.sh`
chose the literature set with `ls -d <resources>/literature/lit-v* | sort -V | tail -1`,
and that picked the file `lit-v1.0.build.log`, not the directory `lit-v1.0`. The literature
set was therefore silently left out of bundles `cxsv-bundle-20261103-33c075a22071`,
`cxsv-bundle-20261106-169503679a60` (the deployed patient bundle) and
`cxsv-bundle-20261108-c40e26528d5f`. Patient runs from them carry no `literature_match`
annotation (`private/cxsv_literature_hits.tsv`). Tiers and scores are unaffected, because
the literature set is annotation only. The script now selects directories only, stops when
no set is found unless `--no-literature` is given, and validation check 9 fails a bundle
without the set.

### 6a. Bundle with panel-v1.0 (the first real release), 2026-11-03

Built with `scripts/make_bundle.sh bundle` and validated with
`tests/lab/validate_bundle.sh` on freya (Singularity CE 4.4.1, Docker 29.9.2),
networking off. The run directory was freya local disk; the bundle was on the
NAS. Records: `validation.log`, `compare_*.tsv`
(saved with the validation deliverables).

| | |
|---|---|
| bundle | `/mnt/sargasso1/CxSV/bundles/cxsv-bundle-20261103-33c075a22071` |
| commit | `33c075a22071dd229d268427600847ec01ea020e` (branch `validation-panel-v1.0`) = dc7f582 (`item10-graded-privacy`, merged as PR #17; contains `panel-v1.0` = da6cb90, item 9 and item 10) + the `locus_depth` fix below. Calling code and `config/calling.yaml` unchanged since dc7f582 (`git diff` empty) |
| panel release | `cxsv-panel-20261102-cd7e1d194cd4`, **not a test release**, tag `panel-v1.0`. Fingerprint `cd7e1d194cd4…`, release source commit da6cb90, 1,634 files, `release_manifest.tsv` sha256 `928d721ac91bd179be195ab8af5cf7324a4b0db09e043dd1aeb799f2fa835326`; identical before and after the bundle builds and the validation |
| image | reused `cxsv-pipeline:ba820f7b2bd3`, id `sha256:b6c1bded…`, SIF sha256 `9f4e049fe5c5eb89717e1f531d284836efc7c306e0103374872c8165c01ff967` (the SIF the panel was built in, checked on every node at preflight) |
| recorded in the bundle | `bundle.json` / `bundle.env`: commit, release id, fingerprint, release source commit, image id, SIF and Docker-archive sha256, env-lock sha256. `MANIFEST.sha256` (1,900 files; its own sha256 `c6dbb153…`) covers the release's `release_manifest.tsv` (`928d721a…`), `release.json` (`fede1a48…`), `resources/resources_manifest.tsv` (`2713b466…`) and the 4 catalogue manifests. `release.json` itself has `runtime.image_id`/`sif_sha256` null (open defect, item 12), so `bundle.json` `built_in_sif` is null too |
| size | 10,524,465,214 bytes (9.9 GiB). By component: reference 6,296,932,150; panel release 2,432,071,602; container 1,511,291,981; resources 229,422,362; catalogues 46,466,039; source 4,691,591; repo 2,928,889 |
| build time | 83 s (61 s copying, the rest checksums); image reused |
| stand-in patient | **HG01868**: 1kGP-LRSC, KHV/EAS, female, ONT R9 LSK110 (dorado 0.5.0 sup), **whole genome**, ~28x (median flank depth at MAPQ ≥ 20). Public, and **not** in the panel (absent from `registry/individuals.tsv`); from the 1kGP ONT bucket's `ALIGNED_TO_CHM13/MINIMAP_ALIGNED_BAMS_2026/` (90,944,951,543-byte BAM, sha256 `61181675…`; minimap2 2.28 to CHM13v2, contigs identical to `hs1.fa.fai`). Copy kept at `datasets/standins/HG01868/` with a README. HG002 was considered, but no HG002 long-read data were on the NAS, and it would have needed downloading and aligning. The earlier stand-in HG00130 is a panel individual now |

| Check | Result |
|---|---|
| 1. bundle integrity | PASS: all 1,900 files (107 s) |
| 2. panel release integrity | PASS: 1,634 files verified, manifest sha256 `928d721a…` |
| 3. network off (Singularity, Docker) + controls | PASS: loopback only; TCP/DNS/HTTPS fail; `/mnt/sargasso1`, conda and the repo clone absent; `/bundle` read-only. The `--no-netns` controls connect |
| 4. patient run, Singularity, bundle only | PASS (2,894 s, 12 cores) |
| 5. patient run, Docker, image loaded from the bundle archive | PASS (2,907 s) |
| 6. patient run outside the bundle (clean worktree at the commit, conda env, NAS release, resources and catalogues) | PASS (2,652 s) |
| 7. records, 4 vs 6 and 5 vs 6 | **FAIL as scripted; explained, no data difference.** Of 57 files: 51 identical; 3 identical after path normalisation; 2 expected run-record differences; **1 differs: `plots/complexity_class_heatmap.pdf`**. Its panel B rows are ordered by base R `sort()`, which depends on the locale: the image sets `LC_ALL=C.UTF-8`; the host job had no `LANG` (POSIX). `CN_plus_orientation_change` and `Chromothripsis_like` swap places (the slice stand-ins had no `CN_plus…` locus). Re-rendering only that plot in run 6 with `LC_ALL=C.UTF-8` gives 52 identical, 3 path-normalised, 2 expected, **0 unexpected**, for both 4 vs 6 and 5 vs 6. Runs 4 and 5 are identical to each other. *Rerun 2026-11-04 with the fixed script (small fixes 2026-11): PASS, see below* |
| 8a/8b. tampered release (one digit / one extra file) | PASS: refused before any job ran |
| total | 8,897 s (2 h 28 min) |

**Check 7 rerun, 2026-11-04 (branch `small-fixes-2026-11`, fix 1).** From
hardlinked copies of the three result trees above (originals untouched), only
`complexity_class_heatmap.pdf` was re-rendered with the fixed
`plot_complexity_class.R`, each in the environment of its run: Singularity
(the bundle SIF), Docker (the bundle image), and the host conda env (POSIX
locale). Comparing them with the same `compare_patient_runs.py` call as check
7 gives, for 4 vs 6 and 5 vs 6 alike, 57 files: 52 identical, 3 identical
after path normalisation, 2 expected run-record differences, **0 unexpected:
PASS**. The rest of the bundle validation was not re-run; the bundle itself
still carries the old script until it is rebuilt from a commit with the fix.

Stand-in result (identical in runs 4, 5 and 6):

- 974 CxSV loci: 0 private, 62 `shared_junctions`, 912 public;
- 0 Tier 1, 5 Tier 2 (4 public, 1 shared_junctions), 22 Tier 3,
  947 Benign_likely.

That is 0 Tier-1 calls for a healthy genome that is not in the panel. Healthy
panel genomes average 0.36 at 20-30x (item-10 calibration, leave-one-family-out).

**Defect found and fixed (PI-approved).** A first bundle,
`cxsv-bundle-20261103-dc7f582c01cc` (from dc7f582, 86 s), passed checks 1-5
on the stand-in, but the whole-genome positive control failed at rule
`locus_depth`. Sniffles2 2.8.1 gave a PAR1 BND the ALT mate contig `chrN`
(`CHR2=chrY`). `scripts/locus_segments.py` emitted depth intervals on it, and
`samtools bedcov` exited 2. Commit `33c075a`:

- skips breakends on contigs absent from the reference `.fai` and logs the count;
- adds a unit test;
- leaves the fingerprint unchanged (`locus_segments.py` is not calling code).

The dc7f582 bundle's validation was stopped after check 5 and the bundle renamed
`*.SUPERSEDED`; it was deleted on 2026-11-04 (PI-approved). The same commit adds
`tests/lab/pos_control_full_genome.sh` and the `STANDIN_RUN` label to
`validate_bundle.sh`. Positive control against this bundle:
docs/POSITIVE_CONTROL.md, "Full-panel result".

Findings, not fixed (outside this validation):

- The heatmap row order depends on the locale. A fix is explicit factor
  levels in `plot_complexity_class.R`, or `LC_ALL=C.UTF-8` for check 6.
- `run_cxsv.sh` gives the container a 64 MB tmpfs `/tmp` (`--containall`,
  `TMPDIR` unset). Whole-genome runs did not need more; the panel build binds a
  disk `/tmp`. *Fixed in small fixes 2026-11 (section 5).*
- Whole-genome CxSV detection loses events inside oversized clusters
  (docs/ROADMAP.md item 13).

### 6b. Bundle with scoring (roadmap items 4–5), 2026-10-29

Built and validated on freya with the same scripts (`make_bundle.sh image|bundle`,
`tests/lab/validate_bundle.sh`). The bundle now also ships `resources/scoring/`
(ClinGen dosage, gnomAD v4.1 constraint, HPO v2026-09-01, CHM13v2 PARs) and
`resources/resources_manifest.tsv`, verified against the manifest before
copying. `bundle_site.yaml` points the scorer at `/bundle/resources/scoring`
and at the catalogues in `/bundle/catalogues`. Full log: `validation.log`
(saved with the item-5 deliverables).

| | |
|---|---|
| bundle | `/mnt/sargasso1/CxSV/bundles/cxsv-bundle-20261029-1dcabd9770f7` |
| commit | `1dcabd9770f70e2e7f5a194ad3b1d370b23bc15e` (branch `item5-scoring`; this docs record was committed after it) |
| size | 8,089,067,107 bytes (7.6 GiB), 238 files in `MANIFEST.sha256`. By component: reference 6,296,932,150; container 1,511,291,981; resources 229,273,577; catalogues 46,466,039; source 2,767,920; repo 1,694,273; panel release 570,431 |
| image | reused, `cxsv-pipeline:ba820f7b2bd3` (Dockerfile and lock unchanged since 144a4ff), SIF sha256 `9f4e049f…` |
| build time | test release 431 s (panel workflow in the SIF, network off), bundle assembly 52 s (35 s copying + 16 s checksums) |
| panel release | `TEST-posctrl10slice-20261029-cd7e1d194cd4`, **test release**, the same 10 sliced genomes as before (slices from the item-5 regression workdir). Fingerprint `cd7e1d194cd4` (was `ad454c604bdd`: the classifier is fingerprinted). It equals the fingerprint of the conda-built regression panel |
| stand-in patient | HG00130 slice, as before |

| Check | Result |
|---|---|
| 1. bundle integrity | PASS: all 238 files (14 s) |
| 2. panel release integrity | PASS: 25 files verified (manifest sha256 `e5f14d9e…`) |
| 3. network off (Singularity, Docker) + controls | PASS: loopback only, TCP/DNS/HTTPS fail, host paths absent, `/bundle` read-only; the controls with `--no-netns` connect |
| 4. patient run, Singularity, bundle only | PASS (85 s) |
| 5. patient run, Docker, image loaded from the bundle archive | PASS (87 s) |
| 6. patient run outside the bundle (clean worktree at the commit, conda env, NAS resources and catalogues) | PASS (85 s) |
| 7. records, 4 vs 6 and 5 vs 6 | PASS, both. Of 58 output files: 52 identical; 4 identical after path normalisation (bundle vs NAS paths of resources and catalogues: `catalogue_check.tsv`, `scoring_params.json`, plus the two from item 2); 2 expected run-record differences; 0 unexpected |
| 8a/8b. tampered release (one digit / one extra file) | PASS: refused before any job ran |

Stand-in result:

- 13 CxSV loci: 11 public and 2 shared_junctions.
- Scores: 12 Benign_likely and 1 Tier 3 (the chr9 *EHMT1* shared_junctions
  hotspot locus).
- The same in all three runs.

No bundle failed or was superseded in this round. At the time,
`/mnt/sargasso1/CxSV/bundles/` held only this bundle and the previous validated one
(`cxsv-bundle-20261028-26085cdd9977`); nothing was removed.

### 6c. First bundle (roadmap item 2), 2026-10-28

Run on freya (Singularity CE 4.4.1, Docker 29.9.2), 2026-10-28, with
`tests/lab/validate_bundle.sh`. Full log: `validation.log` (saved with the
item-2 deliverables).

| | |
|---|---|
| bundle | `/mnt/sargasso1/CxSV/bundles/cxsv-bundle-20261028-26085cdd9977` |
| commit | `26085cdd99777c9d9b1b6e168e92731b332c7cb7` (branch `item2-panel-release`; no tag, see section 1) |
| size | 8,028,013,781 bytes (7.5 GiB), 216 files + 2 manifests. By component: reference 6,296,932,150; container 1,511,291,981; resources 169,239,697; catalogues 46,466,039; source 2,203,948; repo 1,243,273; panel release 571,206 |
| build time | image 89 s (reused: same Dockerfile + lock), test release 389 s (panel workflow in the SIF, network off), bundle assembly 50 s (36 s copying + 14 s checksums): **about 9 min from a commit** |
| panel release | `TEST-posctrl10slice-20261028-ad454c604bdd`, **test release**. 10 ONT-1019 genomes (the first 10 with ≥17x, HG00096 … HG00129), each sliced to 16 regions (34.37 Mb) around the positive-control events, with the E13 spike-in in three of them. Built from the bundle's commit inside the bundle's SIF. 902 junction sites and 31 CxSV loci; 701 sites callable in all 10 individuals, 873 in at least half, 15 in none |
| stand-in patient | HG00130 (the 11th ≥17x ONT-1019 genome, 28.3x, not in the panel), sliced over the same 16 regions (59.6 MB CRAM) |

| Check | Result |
|---|---|
| 1. bundle integrity | PASS: all 216 files match `MANIFEST.sha256` (13 s) |
| 2. panel release integrity | PASS: 25 files verified (manifest sha256 `c3df7c58…`) |
| 3a. network off (Singularity `--containall --net --network none`) | PASS: interfaces `lo` only; TCP to 1.1.1.1:443 and 140.82.112.3:443 "Network is unreachable"; DNS fails; `curl https://github.com` exit 6. `/mnt/sargasso1`, `~/miniforge3` and the repo clone are absent; `/bundle` is read-only |
| 3b. network off (Docker `--network none`) | PASS: the same probe, same result |
| 3c. controls | PASS: the same probe with `--no-netns` sees enp5s0/docker0 (Singularity) or eth0 (Docker), connects and resolves; the probe can detect a network |
| 4. patient run, Singularity, bundle only | PASS (76 s) |
| 5. patient run, Docker, image loaded from `container/cxsv-pipeline.docker.tar.gz` (local tag removed first) | PASS (76 s) |
| 6. patient run outside the bundle: worktree at the same commit (clean), conda env `cxsv-pipeline`, lab paths, the original release directory | PASS (72 s) |
| 7. records, 4 vs 6 and 5 vs 6 | PASS, both. Of 49 output files, 45 are identical (VCF records, Sniffles2 `.snf` decompressed content, PDFs minus timestamps, all tables byte for byte). 2 are identical after path normalisation (`merged_vcf_list.txt` lists output paths; `refcheck.ok` names the MD5 cache path). 2 are run records that differ as expected and agree on status (`panel_release_check.txt`: release path and runtime; `panel_compatibility.txt`: release path; fingerprint COMPATIBLE `ad454c604bdd…` in both). 0 unexpected differences |
| 8a. tampered release: one digit of `panel_junction_carriers.tsv` changed (same size) | PASS: refused when the workflow is parsed, before any job ran ("sha256 differs: panel_junction_carriers.tsv") |
| 8b. tampered release: one unlisted file added | PASS: refused before any job ran ("not in manifest: extra_sites.tsv") |

Stand-in result (identical in all three runs):

- 13 CxSV loci: 11 public and 2 `shared_junctions` (both chr9 hotspot loci,
  `Panel_Hotspot` yes), none private.
- `Panel_Callable` per locus is 3–10 of 10. The chr4:866,851–882,936 locus has
  6 carriers but only 3 individuals callable at all three of its junctions:
  carriers can sit below `min_depth`, which is why both numbers are reported.
- 325 junctions: 304 with panel carriers, 21 without; `panel_callable` is 10
  at 267 of them.

The validation found three faults, all fixed before this record:

- the builder's checksum list included its own temp file;
- the image smoke test left a root-owned `/tmp/.cache`, so Snakemake failed
  under Docker as a normal user;
- Snakemake's `getpass.getuser()` failed for a uid without a passwd entry
  (Docker; `USER`/`LOGNAME` are now set).

The bundles built along the way are kept, renamed `*.DEFECTIVE-manifest` /
`*.SUPERSEDED`, under `/mnt/sargasso1/CxSV/bundles/`. They are not used and can
be deleted.

The lab positive-control regression (`tests/lab/pos_control_regression.sh`,
10 backgrounds, conda env, commit 26085cd) gives the same result as before
item 2:

- 40 of 41 junctions exact, 1 partial, 0 missed;
- 12 of 13 events as expected; E12 is the listed known failure (roadmap
  item 5);
- 0 new failures;
- 668 s run time.

So the patient-side changes (release verification, callability columns)
leave the calls and labels unchanged.

Limits of this validation:

- It ran on freya with the container isolating the run, not on a separate
  clean machine.
- The panel is a test panel: its labels mean nothing outside the 16 sliced
  regions.
- No release tag was made: release tags are the PI's decision. The commit
  hash is recorded instead.
