# panel-v1.3: build record, bundle validation and deployment attempt

Status (2026-11-13):

- **Tagged `panel-v1.3`** (tag object `06a6504c7fc9…` on main `718df799a826`); the release
  `cxsv-panel-20261111-aa64baf08272` is **finalised** (source commit 1c9c80d, calling code
  identical at the tag commit; PR #55).
- **Patients switched** to the tagged bundle `cxsv-bundle-20261113-718df799a826` (section 6).
  Rollback: `cxsv-bundle-20261109-84d380e08cc4` (panel-v1.2).
- The first deployment attempt (2026-11-12, bundle `1c9c80d0e2c6`) stopped on two findings
  (section 4); they were resolved by PI decisions 21 (scorer quarantine of non-reference contigs,
  PR #53) and 22 (locus end = contig length + 1 accepted, PR #54).

Plan: [PLAN_panel_v1.3.md](../PLAN_panel_v1.3.md). Design check: [design_check_v13/README.md](../design_check_v13/README.md). ROADMAP item 16, decisions 19 and 21.

## 1. Release build (Phase 2)

| | |
|---|---|
| build root | `/mnt/sargasso1/CxSV/build/panel-v1.3` |
| source | main `1c9c80d0e2c608c6709eb8aaea19b2c61fa34298` (git archive; `code/COMMIT`; release.json `source.dirty` false) |
| calling fingerprint | `aa64baf08272fc45a730a54ac19e706cf804abbe8d884199e5d97dc0dac60533` (code digest `fa204bc02219…`) |
| manifest | `release_manifest.tsv` sha256 `25dbee9ee5effa626c8b2f0ced8025ac2072d16896a520c15f14c0434254cb42`; 1,640 files, 2,465,592,378 bytes; VERIFIED |
| image used for the build | `cxsv-pipeline:ba820f7b2bd3`, SIF `9f4e049fe5c5…` (as for panel-v1.2-r2) |
| OLD (reused upstream) | `build/panel-v1.2-r2` (cxsv-panel-20261108-fdac9ff80401, tag panel-v1.2) |
| depth | copied (not hard-linked) from `build/cn_depth_v1` (item a); 1,619/1,619 sha256 equal, nlink 1 |

Procedure (docs/PANEL_BUILD.md, "Downstream-only rebuild"):

1. `git archive 1c9c80d` → `code/`. The fingerprint check matched.
2. The held caller code was unchanged since 84d380e: `git diff` was empty for `calling.smk`, `per_sample.smk`, `caller_checks.smk`, `panel_callability.smk`, `panel_aggregate.smk`, `callable_mask.sh`, `check_svim_output.py` and the conda lock. Only the CN-track files and `cn_params` changed.
3. Guard snapshots were taken before the build (panel-v1.1-r2, panel-v1.2-r2, cn_depth_v1).
4. `make_aggregate_yaml.py`, then `rebuild_downstream.sh --plan-old` (319 s), then `link_upstream.py --apply` (153 s).
5. **Tooling gap:**
   - `link_upstream.py` hard-linked 7,190 files that no dry run plans but that the rebuild rewrites:
     - `cn_params.json`, a params-only output;
     - 3,951 `cn/work/*` files, which `cn_linkage.py` writes in place;
     - 2 × 1,619 cn_depth side files (`.csi`, mosdepth summary).
   - They were unlinked by hand (`fix_cn_sidefiles.sh`, after checking that every one was a hard link to OLD). The inventory went to 58,646 files (1,410.5 GB).
   - **Fixed in PR #50:** `link_upstream.py` excludes params-only outputs and every file of a planned `cn/` directory. `rebuild_downstream.sh` stops on a hard-linked side file. There are 4 new unit tests.
   - A read-only check of the fixed tool on the real trees reproduces the manual fix exactly: 58,646 files to link; 1 params-only + 7,189 cn side files left out.
6. `rebuild_downstream.sh OLD NEW 20`: 18,041 s, 17,835 jobs, one run. The inventory verify found 0 changed files.
7. Guard snapshots after the build: **0 differences** for panel-v1.1-r2, panel-v1.2-r2 and cn_depth_v1.

| | panel-v1.2-r2 | panel-v1.3 |
|---|---|---|
| junction sites | 1,304,862 | 1,304,862 |
| CxSV calls | 1,733,997 | 1,733,865 |
| panel CxSV loci | 13,053 | 11,899 |
| CN segments | 37,562 | 37,562 (panel_cn_segments.tsv byte-identical) |
| chrX/chrY public counting excluded (D2 discordance) | – | 4 libraries (NA19238, HG03499; LRSC + ONT-1019 each) |

### Healthy-genome LOFO (full panel, 1,619 genomes; the binding gate)

| stratum | n | v1.2-r2 | **v1.3** |
|---|---|---|---|
| all | 1,619 | 0.198 | **0.2075** (336 Tier 1 loci) |
| ONT LRSC R10 | 247 | 0.247 | 0.2551 |
| ONT LRSC R9 | 220 | 0.159 | 0.1682 |
| ONT 1kGP-1019 | 1,017 | 0.204 | 0.2124 |
| PacBio HiFi | 135 | 0.133 | 0.1481 |
| male / female | 795 / 806 | 0.224 / 0.171 | 0.2327 / 0.1811 |

Stop rules:

- overall 0.2075 ≤ 0.248: pass;
- public Tier 1 loci = 0: pass.

The R10 stratum (0.255) is above 0.248, but the PI's gate applies to the overall rate. Tier 1 changes were +18 / −3:

- (h) type-agnostic caller agreement: +17;
- (g) CN locus span: +1 / −2 (v1.2 malformed CN ids);
- (D1) public cap: −1 (HG00369);
- (D2) and (D3): 0.

Literature test A through the release:

- junction level: 8/63 called, 7 Tier 1;
- CN level: 40/63 called, 34 Tier 1;
- DUP-TRP/INV-DUP: 26/26 called, 23 Tier 1;
- per event, identical to the panel-v1.2-r2 branch run.

### Design decisions (Phase 1; design_check_v13)

- **D1 (public tier ceiling):** public loci are capped at Tier 2 (BOC_2057 *BRSK2*, HG00369).
- **D2 (sex-chromosome public rule):** a chrX/chrY segment is public only when ≥ 2 sex-matched individuals carry it. Libraries of a sex-discordant individual count for no one. L01 is private at female bounds 1.12, 1.18 and 1.25.
- **D3 (depth-noise gate):** sigma_max = median + 3.5 × scaled MAD per stratum (R10 0.233, R9 0.220, ONT-1019 0.401, PacBio 0.198). It flags BOC_5883 and BOC_7308.
- **D4:** cuteSV maximum size is kept at 100 kb (an upstream re-run is a separate item).
- **D5:** `nested_min_z` is kept at 6.0.
- **D1/h re-scoring:** BOC_4433 *NF1* agreement went from 0/2 to 2/2 (score 1.66).
- **PI decision (2026-11-12, subset LOFO):** PROCEED. The subset LOFO 0.2566 → 0.2832 on the enriched 113-genome subset was accepted. The binding gate is the full-panel LOFO ≤ 0.248, which passed (0.2075).

## 2. Bundle (Phase 3)

`/mnt/sargasso1/CxSV/bundles/cxsv-bundle-20261112-1c9c80d0e2c6`:

- built from a clean clone of main 1c9c80d with `--literature lit-v1.0`;
- 2,296 files, 10,600,341,417 bytes;
- `MANIFEST.sha256` sha256 `8dc8dfda518518044a408b07f34926aece78194760aa97da86e0c823114549f1`;
- bundle.json `panel_release.release_manifest_sha256` = `25dbee9e…cb42` (item n); `git_tag` null.

**New image** `cxsv-pipeline:ddfe0617e229` (image id `sha256:b3b02c2c…`, SIF `b1aed2bb…bf6e`):

- conda explicit lock (unchanged) plus the hash-pinned PyPI wheels of the clinician report (`envs/cxsv-pipeline.pip.txt`: python-docx 1.1.2, reportlab 4.2.5, lxml 6.1.3, chardet 7.7.0).
- **Defect:** `make_bundle.sh` failed on 1c9c80d. `container/Dockerfile.dockerignore` (`*` plus a re-admit list) did not re-admit the pip lock, so the docker build's `COPY` failed. CI does not build the image.
  - Fixed in **PR #49** (merged, c2205aa), with a test that every `COPY` source is re-admitted.
  - The image key covers only the Dockerfile, the conda lock and the pip lock, so the fix does not change it. The image was built from the fix commit 2d7d6f9, and the bundle from 1c9c80d reuses it.
- The lab conda env `cxsv-pipeline` on freya gained the same four wheels (`pip install --no-deps --require-hashes -r envs/cxsv-pipeline.pip.txt`; nothing was replaced), so that validation check 6 (repository run) can render the report.

## 3. Validation and acceptance tests

### Bundle validation (freya; `tests/lab/validate_bundle.sh` + `validate_v13_extras.sh`)

- **PASS:**
  - checks 1–6 and 8–9: integrity, release, network off, offline Singularity, offline Docker, repository run, tampered release refused, lit-v1.0;
  - check 10, clinician report: docx + pdf + graphic; 2 pages; sex call printed; Singularity = Docker;
  - check 11, intake offline on a 20,000-read modBAM subset: verdict PASS; 20,000 = 19,957 mapped + 43 unmapped; one `@RG` with SM = sample_id; minimap2 2.24-r1122 `@PG`; MM/ML kept; a PacBio row refused;
  - check 12, **fingerprint recomputed in the bundle = `aa64baf08272…`, COMPATIBLE**;
  - check 13, bundle.json manifest sha256;
  - check 14, **HG01868 0 Tier 1** (1,077 loci, 4 Tier 2).
- **Check 7 FAIL (validation tool at 1c9c80d):**
  - The bundle runs and the repository run differ only in the clinician report's run-identity row (pipeline bundle id / "not run from an offline bundle"; "code 1c9c80d0e2c6").
  - All 82 other outputs are identical or differ only as expected (76 identical, 4 after path normalisation, 2 expected run records).
  - **PR #51** classifies this row as an expected run-record difference (report text compared with the identity masked). With it: 0 unexpected differences.

### Full-genome positive control (sample 3559 + E13, L01–L04; through the bundle)

- 6,258 s.
- **All 13 planted Tier 1 events stay Tier 1, with scores identical to v1.2.** Verdicts are 15/17 (E13 and L03 are the known failures).
- Junctions: 86 exact, 1 partial, 0 missed of 87.
- 16 Tier 1 loci, all inside planted events.
- L01 locus id `chrX:101752400-102493000` (v1.2 `…-102492226`; item g). It is **private at female chrX bounds 1.12, 1.18 and 1.25**:
  - v1.3 rule: 0 matches at every bound;
  - v1.2 rule: 1 match at 1.18 and at 1.25, NA19238, which is now excluded as sex-discordant.

### Nine test genomes (raw modBAM → intake → patient workflow + clinician report)

NAS `/mnt/sargasso1/CxSV/patient_test_runs/monkeymeth_v13_20261112/`; Drive mirror `CxSV-v2/patient_test_runs/monkeymeth_v13_20261112/` (992 files, MD5SUMS verified on Drive).

| criterion | value | result |
|---|---|---|
| 9 genomes, exit 0, network off, COMPATIBLE | 9/9 | PASS |
| BOC_4433 *NF1* chr17:31,902,782-32,067,453 Tier 1, caller agreement ≥ 1/2 | Tier 1, 1.66, 2/2, private | PASS |
| BOC_2057 *BRSK2* chr11:1,465,980-1,509,727 ≤ Tier 2 | Tier 2 (Tier 1 by points, 1.06; public 65/1,582, capped) | PASS |
| sex call in every report | 9/9 | PASS |
| intake read counts (reads in = primary mapped + unmapped) | 9/9 | PASS |
| 2-page clinician report, docx + pdf | 9/9 | PASS |
| **no malformed locus IDs anywhere** | **BOC_5883 `chrN:62452837-62454007`** (BND_cluster, Benign_likely) | **FAIL** |

Tier 1: **1 locus in 9 genomes (0.111 per genome)** at v1.3, versus 2 at v1.2. Healthy comparators: R10 0.2551, all 0.2075. Tier 2 rose by one in BOC_2057 (*BRSK2*, capped), BOC_3537 and BOC_3559 (two public Chromothripsis_like loci now at 0.50 and 0.59). D3 flags BOC_5883 (sigma 0.322) and BOC_7308 (0.237).

`chrN` comes from the Sniffles2 2.8.1 PAR artefact (known since panel-v1.0; ROADMAP item 18):

- For BNDs near the chrY PAR, Sniffles2 writes the ALT mate contig as `chrN` while CHR2=chrY, and SURVIVOR keeps `chrN`.
- The v1.0 fix only drops depth intervals on contigs that are not in the reference.
- Spread: `junctions.tsv` in 5 of 9 test genomes; one scored locus in BOC_5883 (two at v1.2); 266 junction sites and 18 panel loci in the v1.3 release (v1.2: 266 / 21).

## 4. First deployment attempt (2026-11-12): stopped

**STOPPED before tagging** (resolved in section 6). Failed criteria:

1. Malformed locus id `chrN:…` (BOC_5883).
2. Bundle validation check 7 with the 1c9c80d tool (an expected report difference; PR #51).

Actions not taken: no tag `panel-v1.3`; no `write_release_metadata.py --finalise`; no tagged-bundle rebuild; no patient switch.

The patient configuration is unchanged: `cxsv-bundle-20261109-84d380e08cc4` (tag panel-v1.2), with `cxsv-bundle-20261106-169503679a60` (panel-v1.1) as the older rollback.

Bundles on the NAS (none deleted; deletion needs PI approval):

| bundle | role | bytes |
|---|---|---|
| cxsv-bundle-20261109-84d380e08cc4 | **current patient bundle** (tag panel-v1.2) | 10,577,419,887 |
| cxsv-bundle-20261112-1c9c80d0e2c6 | v1.3 candidate, untagged, validated (this record) | 10,600,341,417 |
| cxsv-bundle-20261106-169503679a60 | superseded (panel-v1.1), older rollback | 10,562,414,920 |
| `_images/img-ddfe0617e229` / `img-ba820f7b2bd3` | v1.3 image / v1.1–v1.2 image | 1,529,204,324 / 1,511,405,022 |

PI options (decision 21):

- **(a)** Accept both findings: the `chrN` locus is Benign and a known caller artefact, and check 7 is a tool false positive. Then tag panel-v1.3 on 1c9c80d, finalise, rebuild the bundle from the tag (byte identity), and switch patients.
- **(b)** Fix `chrN` first: normalise the BND mate contig to CHR2, or drop such records, in fingerprinted calling code. This makes a new release, panel-v1.4.
- **(c)** Report-level fix only: quarantine loci on contigs absent from the reference in the scorer. This keeps the release but needs a new bundle and validation.

## 5. Artifacts and scripts

- Scripts on the NAS, `patient_test_runs/monkeymeth_v13_20261112/code/`:
  - `run_regression_v13.sh`: claim-locked intake + patient per genome;
  - `validate_v13_extras.sh`: bundle checks 10–14;
  - `summarise_regression_v13.py`: criteria and tables;
  - `collect_derived_v13.sh`: Drive mirror.
- Scripts on freya: `/home/minknow/scratch/v13/bundle/` (validation `bundle_validation/`, PC `pc_full_L/`, L01 bounds `l01_bounds/`).

## 6. Rounds 2–3, tag, finalised release, tagged bundle and patient switch (2026-11-13)

**Round 2: scorer quarantine (decision 21, option c; PR #53, main 5a50dde).** Loci on contigs
absent from the reference are removed, and non-reference junctions are dropped from scored loci
(SCORING.md 6f; `scoring/quarantined_loci.tsv`; one QC line in the score and clinician reports).
Validation check 10 (locus ids well-formed) was added. Healthy LOFO: 0 of 336 Tier 1 loci touch a
non-reference contig; re-scored with the quarantine, Tier 1 stays 118 → 118. Bundle
`cxsv-bundle-20261112-5a50ddef8710` passed everything except the strict check 10: 9 locus ids
in 7 of the nine genomes, 2 in HG01868 and 1 in the positive control end at contig length + 1
(all Benign_likely).

**Contig-end ids (decision 22; PR #54, main 718df79).** `vcf_to_bed.py` writes a BND breakend as
`[POS - bnd_window, POS + bnd_window)` without clamping to the contig, so a breakend on the last
base gives end = length + 1 (24 panel loci in panel-v1.2-r2 and in panel-v1.3; 603 healthy-LOFO
loci, 0 Tier 1). Check 10 and the regression criterion now accept `end <= length + 1` and still
reject chrN / other non-reference contigs, zero length, a start on or past the contig end and
anything further past it. The clamp in calling code is a v1.4 item (ROADMAP item 19).

**Round 3: validated bundle `cxsv-bundle-20261112-718df799a826`** (main 718df79, release
unchanged, lit-v1.0, image `cxsv-pipeline:ddfe0617e229` / SIF `b1aed2bb…`; 2,302 files,
10,600,572,975 bytes; MANIFEST.sha256 sha256 `e58fbdd6…0967`):

- Bundle validation (`validate_bundle.sh`, 8,450 s): checks 1–10 PASS, failures 0. Check 7: no
  unexpected difference. Check 10: 2,155 locus ids, 0 malformed, in each of the three runs.
- Extras 10–14 PASS: 2-page clinician report with the sex call (Singularity = Docker), offline
  intake, fingerprint `aa64baf08272…` COMPATIBLE, `bundle.json` release-manifest sha256,
  **HG01868 0 Tier 1** (1,047 Benign_likely, 26 Tier 3, 4 Tier 2). HG01868 `cxsv_scores.tsv`,
  `cxsv_score_features.tsv`, `quarantined_loci.tsv` and `cn_loci.tsv` are byte-identical to the
  5a50dde bundle's runs.
- Positive control (full genome, callers re-run, 6,122 s): `score.events.tsv` (md5 `03e6d798…`),
  `score.junctions.tsv`, `score.offtarget.tsv`, `cn/segments.tsv`, `cn/cn_loci.tsv`, `cn/stats.json`
  and the scoring tables are byte-identical to the 5a50dde run; CN segments are identical to the
  1c9c80d candidate. 15/17 events pass (E13, L03 known), 16 Tier 1 all in planted events,
  86/87 junctions exact; **L01 Tier 1, `CN_linked_gain`, 1.53, private** (bounds 1.12/1.18/1.25
  unchanged, CN track identical).
- Nine genomes re-scored through the bundle (scoring + report stages only): 9/9 exit 0,
  COMPATIBLE, network off; every criterion PASS (NF1 Tier 1 1.66, callers 2/2; BRSK2 Tier 2;
  sex call, intake, 2-page report and quarantine line 9/9; 0 malformed ids). 839 of 866 output
  files are byte-identical to the 5a50dde runs; the other 27 (clinician reports,
  `panel_release_check.txt`) differ only in the bundle id and code commit.

**Tag (PI).** Annotated tag `panel-v1.3`, object `06a6504c7fc9eefec521d07a740fb6935c117bec`, on
`718df799a826e14a19982a11eda47be2d5489a98` (tagger Chris Faulk). The message records the release
id, source commit 1c9c80d, the calling fingerprint and decisions 21/22.

**Finalise.** `write_release_metadata.py --finalise … --git-tag panel-v1.3
--allow-descendant-tag` (tool from main a53a6fc, PR #55). The tag commit descends from the release
source 1c9c80d and the 17 calling-code files and `config/calling.yaml` are identical.
`source.commit` stays 1c9c80d; `source.git_tag_commit` 718df79 and `source.git_tag_relation` are
recorded; `finalised_utc` 2026-11-13T14:41:13Z. Only the `release.json` row of the manifest changed:
manifest sha256 `25dbee9e…` → `6ae39d767bc28fcf11d54e5286a55bb89fccedc6ac8b5deb75676fcab6cc13a4`;
`release.json` sha256 `043b37c1c8cb…`; 1,640 files VERIFIED. Pre-finalise copies:
`build/panel-v1.3/provenance/pre_finalise/`.

**Tagged bundle.** `make_bundle.sh bundle --commit panel-v1.3` from a node-local clone at the tag
(81 s). The UTC date changed, so the bundle id is new; the untagged bundle was renamed (not
deleted) to `cxsv-bundle-20261112-718df799a826.pretag`.

| | |
|---|---|
| bundle | `/mnt/sargasso1/CxSV/bundles/cxsv-bundle-20261113-718df799a826` |
| size | 10,600,573,455 bytes; 2,302 files in `MANIFEST.sha256` (sha256 `8cf5d5f58bb3f5185cd8f1a316d5c0bf7c7e4e8b526a9d849e4a4112a12d29ae`) |
| `bundle.json` | `source.git_tag` `panel-v1.3` (object `06a6504c…`); `panel_release.git_tag` `panel-v1.3`, `finalised_utc` 2026-11-13T14:41:13Z, `release_manifest_sha256` `6ae39d76…`, `release_json_sha256` `043b37c1…` |
| checks | `run_cxsv.sh verify-bundle`: all 2,302 files verified; `verify-release`: VERIFIED |
| vs the validated `.pretag` bundle | 2,295 of 2,302 manifest rows identical, the same file set. Differences: `bundle.json` (id, created time, tag fields, release-manifest/json sha256, sizes), `bundle.env` and `README.md` (bundle id), the embedded `release.json` and `release_manifest.tsv` (byte-identical to the finalised release; all other manifest rows identical), `source/cxsv-718df799a826.bundle` (same single head 718df79; pack bytes only) and `repo/workflow/lib/__pycache__/cxsv_bundle_meta.cpython-310.pyc` (bytecode cache written during the build; header mtime, `.py` source identical). Every other `repo/` file, the image, reference, catalogues and literature set are byte-identical. |

**Patient switch.** The current patient configuration is now:

| | |
|---|---|
| bundle | `/mnt/sargasso1/CxSV/bundles/cxsv-bundle-20261113-718df799a826` (tag `panel-v1.3`, commit 718df79) |
| panel release | `cxsv-panel-20261111-aa64baf08272` (panel-v1.3, finalised), embedded in the bundle |
| replaces | `cxsv-bundle-20261109-84d380e08cc4` (panel-v1.2), kept as the rollback |
| new for patients | public loci capped at Tier 2; type-agnostic caller agreement; depth-noise gate; sex call in the report; clinician report rule (docx + pdf + Tier 1 graphic); intake module; non-reference-contig quarantine |

Patient-run instructions are unchanged (RELEASE_AND_BUNDLE.md section 5): `run_cxsv.sh
verify-bundle`, `run_cxsv.sh verify-release`, then `run_cxsv.sh patient …`; read the sex call and
the quarantine line in the report.

Bundles on the NAS (none deleted; deletion needs PI approval):

| bundle | role | bytes |
|---|---|---|
| cxsv-bundle-20261113-718df799a826 | **current patient bundle** (tag panel-v1.3) | 10,600,573,455 |
| cxsv-bundle-20261109-84d380e08cc4 | rollback (tag panel-v1.2) | 10,577,419,887 |
| cxsv-bundle-20261112-718df799a826.pretag | superseded: validated untagged build of 718df79 | 10,600,572,975 |
| cxsv-bundle-20261112-5a50ddef8710 | superseded: round 2 (strict check 10 failed) | 10,600,554,420 |
| cxsv-bundle-20261112-1c9c80d0e2c6 | superseded: round 1 | 10,600,341,417 |
| cxsv-bundle-20261106-169503679a60 | superseded (panel-v1.1) | 10,562,414,920 |

Records on the NAS: `patient_test_runs/monkeymeth_v13_20261112/final_bundle_5a50ddef8710/` (round 2)
and `final_bundle_718df799a826/` (round 3: re-scored runs, `summary/`, `identity_vs_5a50dde.tsv`);
freya `/home/minknow/scratch/v13/r3/` (validation, PC).
