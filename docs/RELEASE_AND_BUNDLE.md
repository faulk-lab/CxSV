# Panel releases and the offline bundle (roadmap item 2)

The patient scorer must run in a restricted environment, on a machine with no
network, from material whose origin can be checked. This document defines:

1. the **panel release**: what the panel builder hands over;
2. **callability**: how a release says where the panel could have seen a junction;
3. the **container image**: one pinned environment for both workflows;
4. the **offline bundle**: everything needed to score a patient with no network;
5. how to **build** and **run** each of them;
6. the **validation record** of the first bundle.

## 1. Panel release

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

**Release id:** `<release.name>-<YYYYMMDD>-<first 12 hex digits of the
calling fingerprint>`, e.g. `TEST-posctrl10slice-20261028-ad454c604bdd`. The
fingerprint part says which calling a patient must match. The name and date
tell apart releases that share a calling but differ in samples. The id is fixed
the first time it is made (`results/panel/panel/release_id.txt`).
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
- `runtime`: the container image tag, **image id**, **SIF sha256** and
  environment-lock sha256 when built in the container; otherwise
  `CONDA_PREFIX`;
- `files`: the file list.

Release tags are **not** created by the workflow. The commit hash is recorded
instead; tagging a release is the PI's decision.

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
| `resources/` | `genes.final.bed`, `exons.filtered.bed`, `repeats.final.bed`, `segdup.final.bed` |
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

**Open question for the restricted environment's IT** (unchanged): which
workflow engine and container runtime are allowed. The bundle assumes
Snakemake (inside the image) plus Singularity/Apptainer, and also ships a
Docker archive.

## 6. Validation record

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
