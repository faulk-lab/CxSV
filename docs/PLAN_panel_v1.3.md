# Build plan: panel-v1.3

*Drafted 2026-11-11 from roadmap item 16, the panel-v1.2 build record
([panel_v1.2/README.md](panel_v1.2/README.md)) and the nine-genome test run (roadmap item 17).
Nothing here starts until the PI approves the plan and the decisions in section 2.*

## 1. What v1.3 must fix

| # | Issue | Evidence | Component | Fingerprinted? |
|---|---|---|---|---|
| a | Precomputed depth is hard-linked, so a rebuild changes `cn_depth_v1` mtimes and can stop the run | v1.2 rebuild stopped once (HG01137); 1,619 mtimes changed, content intact | `scripts/cn_depth.sh` | yes |
| b | One sex-matched panel individual makes a chrX/chrY copy-number segment public | PC L01 is private only at female chrX bound 1.12; one library (NA19238 LRSC, chrX 1.17) makes it public at 1.18/1.25 | `cn_params.public_min_genomes` (1), `cn_linkage.py` | yes |
| c | Sex call and possible sex-chromosome aneuploidy are not in the patient report; two sex inferences disagree | report prints only ploidy sex and chrX ratio; `score_cxsv.py` runs its own flank-depth `infer_sex` (0.75/0.85) | `score_cxsv.py`, report | no |
| d | CN literature test A runs through the Phase 1 prototype, not the release | `calibration/cn_linkage/literature_cn.py` uses `phase1/seg` and the 110-genome PoN | calibration | no |
| e | Public loci can reach Tier 1 | BOC_2057 BRSK2 locus: 65/1,582 panel carriers, HPRC AF 0.52, score 1.06 (dosage +0.90, chromothripsis +0.85, public −1.00); 1 public Tier 1 locus in the healthy LOFO | `scoring.yaml` tiers / `score_cxsv.py` | no |
| f | No depth-noise QC: noisy libraries produce long CN-only candidate lists and unreliable CN calls | BOC_5883 CN sigma 0.32 (341 CN-only candidates), BOC_7308 0.24 (110); others 0.16–0.21; library QC PASS for both | CN stats → QC/report, conservative mode | no (gate); sigma already computed |
| g | Malformed copy-number locus id | BOC_5883 Tier 2 CN_linked_loss over 17p13.3 reported as `chr17:0-2` (Start 0, End 2) | `cxsv_cnlink.py` / locus construction | yes |
| h | Caller agreement counts a junction as unsupported when callers label it with different SV types | BOC_4433 NF1 deletion: Sniffles2 DEL, SVIM INV, cuteSV none → 0/2 junctions agreed (−0.10) | `score_cxsv.py` caller-agreement feature | no |
| i | cuteSV runs with its default maximum SV size (100 kb), so it cannot call larger events | BOC_4433 NF1 deletion (165 kb) has no cuteSV call | `rules/calling.smk` (caller parameters) | yes, and upstream |
| j | Nested gain-in-gain detection is sensitive to the threshold on the female chrX reference | L01 nested triplication z 5.73 < `nested_min_z` 6.0 (it was 6.14 on the pooled reference); NA19024 lost a nested DUP-TRP | `cn_params.nested_min_z` | yes |
| k | No standard patient intake: the test genomes needed ad-hoc realignment scripts | 9 modbams realigned by hand to hs1 with minimap2 2.24; BOC_3559 CRAM carries SM `pos_control_3559` | new `workflow/intake` (or `scripts/intake/`) | no (alignment is upstream of calling; recipe must match the panel) |
| l | The clinician report is generated outside the pipeline | `make_clinician_reports.py` (test run) reads run outputs after the fact | new reporting rule in the patient workflow | no |
| m | SVIM can run out of memory | BOC_3559 attempt 1 on athena: SVIM `std::bad_alloc`; rerun on freya passed | rule resources / node qualification | no |
| n | `bundle.json` does not record the git tag or the release-manifest sha256 | tag only in the embedded `release.json` | `scripts/make_bundle.sh` | no |
| o | `docs/PANEL_BUILD.md` lacks the panel-v1.2 measurements | edit blocked by the safety check 2026-11-08; content is in `panel_v1.2/README.md` | docs | no |

Not in scope for v1.3 (decision in section 2): OMIM/HPO features (licence and phenotype
input; both are 0 today), exclusion of the residual noisy panel libraries (HG01879, HG00268,
HG01867; library QC flags them and conservative mode handles them), the read-level axes of item 8
(on hold), the known positive-control failures E13 and L03.

## 2. Decisions for the PI

| # | Decision | Options | Recommendation |
|---|---|---|---|
| D1 | Public loci and Tier 1 (e) | (1) tier ceiling: a locus whose panel status is `public` cannot be above Tier 2; (2) larger public penalty; (3) leave as is | (1): simple, explainable in the report ("seen in healthy genomes"), removes the BRSK2-type call without retuning every weight |
| D2 | Sex-chromosome CN public rule (b) | (1) `public_min_genomes` 2 on chrX/chrY only; (2) drop libraries whose chrX/chrY ratio disagrees with another library of the same individual; (3) both | (3): (1) makes one artefactual library insufficient, (2) removes the known cause; both are measured in the design check |
| D3 | Depth-noise gate (f) | threshold on CN segmentation sigma, calibrated per platform on the panel (e.g. robust z > 4), applied as conservative mode to CN-only and CN-linked classes and printed in the report | approve; threshold set from the panel distribution in the design check, not from the test genomes |
| D4 | cuteSV maximum size (i) | (1) keep 100 kb (no upstream rerun); (2) raise to 1 Mb (or the long-span limit) and re-run cuteSV + merge for all 1,619 panel genomes | (1) for v1.3, with (h) fixing the agreement penalty; revisit (2) as its own item: it re-runs an upstream caller on 80.8 TB (estimate to be measured on 10 genomes first) |
| D5 | `nested_min_z` (j) | keep 6.0, or a value chosen from a sweep (5.0–6.0) on the healthy panel and the positive control | run the sweep in the design check; change only if healthy nested calls do not rise |
| D6 | Out-of-scope items (section 1) | defer | defer |

## 3. Work plan

Phase 1 — code (three branches, one PR each, unit tests and synthetic CI green before merge):

1. **cn-fixes** (fingerprinted): (a) copy precomputed depth (`cp --reflink=auto`, verify sha256);
   (b) per D2; (g) locus construction for CN_linked loci, with a unit test that rejects
   zero-length or out-of-segment coordinates; (j) only if D5 changes the value.
2. **scoring-report** (not fingerprinted): (c) report section "Sample QC and sex check"
   (`sex_call`, `x_ratio`, `y_ratio`, ambiguous-sex flag), remove `infer_sex`; (e) per D1;
   (f) per D3, using the sigma the CN track already writes; (h) junction-level caller agreement
   independent of SV type label (breakend positions within the matching tolerance);
   (l) clinician-report rule (docx + pdf) built from the test-run generator, with the tier
   explanation box and the research-use footer.
3. **intake-infra** (not fingerprinted): (k) intake rule or script: uBAM/modBAM/FASTQ/CRAM
   from any aligner → hs1 CRAM with the panel recipe (minimap2 2.24-r1122
   `-ax map-ont --eqx -k 17`, `-y` with MM/ML/MN/qs, read-count check, SM tag from the sample
   sheet), shipped in the bundle; (m) explicit memory request for SVIM and one retry on
   another node; (n) git tag and release-manifest sha256 in `bundle.json`; (o) re-raise the
   PANEL_BUILD.md edit with the PI's approval.

Phase 2 — design check (before any panel rebuild; Phase 1 subset of 110 genomes + the 14
L01 matchers + the 9 test genomes + the positive control):

- D2 variants: chrX/chrY public segments by sex; L01 private at female bounds 1.12, 1.18, 1.25.
- D3: CN sigma distribution per platform over all 1,619 panel libraries (from the existing
  CN stats); proposed threshold; which panel and test genomes it flags.
- D5: nested calls per genome and PC L01/E-events across `nested_min_z` 5.0–6.0.
- D1 and (h): re-score the 9 test genomes and the subset; BRSK2 ≤ Tier 2; BOC_4433 NF1
  caller agreement ≥ 1/2 and still Tier 1.

Checkpoint with the PI before Phase 3.

Phase 3 — release:

1. Downstream-only rebuild from main into `build/panel-v1.3` with the committed tooling
   (upstream caller outputs hard-linked from `panel-v1.2-r2` unless D4 = (2); depth copied from
   `cn_depth_v1`; guard snapshots of `panel-v1.1-r2` and `panel-v1.2-r2`). Measured cost of the
   v1.2 rebuild: 17,143 s (4.8 h) on freya.
2. Release VERIFIED; full-panel LOFO (9,737 s for v1.2) by platform, sex and library QC.
3. CN literature test A through the release pipeline (d), reported next to the prototype numbers.
4. Checkpoint with the PI.

Phase 4 — bundle and deployment:

1. Bundle from main with lit-v1.0; recompute the fingerprint in the bundle.
2. Validation (9 checks + new report and intake checks; 8,413 s for v1.2); HG01868.
3. Full-genome positive control (6,184 s for v1.2).
4. Regression set: the nine test genomes (≈ 4,100–10,200 s each, run in parallel on the fleet).
5. Checkpoint; the PI tags `panel-v1.3`; finalise; tagged bundle rebuild and identity check;
   patient switch; superseded bundles listed for the PI; docs and Drive mirror;
   manuscript numbers updated to v1.3.

Estimated wall time without D4 (2): about 2 days of compute across the fleet, dominated by the
rebuild, LOFO, validation and regression runs; code and design checks come first.

## 4. Acceptance criteria and stop rules

| Check | Pass | Stop |
|---|---|---|
| Healthy LOFO Tier 1 per genome (all) | ≤ 0.198 (v1.2) | > 0.248 (v1.2 + 0.05) |
| Healthy LOFO public Tier 1 loci | 0 (D1) | > 0 |
| PC planted events | E01–E09, L01, L02, L04 Tier 1; E10 below Tier 2; E11, E12 not CxSV; L01 private at female bounds 1.12–1.25 | any planted Tier 1 event lost |
| HG01868 | 0 Tier 1 | any Tier 1 |
| Nine test genomes | BOC_2057 BRSK2 ≤ Tier 2; BOC_4433 NF1 Tier 1 and caller agreement ≥ 1/2; no malformed locus id; BOC_5883 and BOC_7308 flagged by the depth-noise gate if D3 sets it there; sex call printed for all nine | NF1 lost, or any new Tier 1 without review |
| Bundle | all validation checks pass offline; fingerprint in bundle = release | any failure |
| Literature test A | junction and CN counts reported (no threshold; v1.2: 8/63 junction, 41/63 CN with the prototype) | — |
