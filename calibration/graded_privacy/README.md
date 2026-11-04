# Graded privacy calibration (roadmap item 10)

Tunes the graded privacy penalty of `config/scoring.yaml`
(`weights.panel_privacy.graded`) on the **background panel only**: every panel
genome's own CxSV loci are scored with the unchanged patient scorer against the
panel with that genome's individual (LOIO) or whole family (LOFO) removed, and
the rule that gives healthy genomes the fewest Tier-1 calls is chosen under a
chance-match constraint. The positive control is only a check, never an input.
Results and the chosen parameters: `docs/SCORING.md` (section "Graded privacy")
and the calibration report.

Nothing here changes the panel release or the calling fingerprint. Paths are
set in `config.sh` (override by environment).

| Step | Script | Where | Output under `$CAL_ROOT` |
|---|---|---|---|
| 1 families | `families.py` (in `stage_a.sh`) | 1 node | `families.tsv` |
| 2 per-locus leave-out callability | `locus_callability.py` (in `stage_a.sh`) | 1 node, 16 threads, ~13 min | `locus_callable.tsv.gz` |
| 3 leave-out statuses | `leaveout_status.py` (in `stage_a.sh`) | 1 node, 8 threads, ~16 min, ~8 GB | `leaveout/<cohort>/<sample>/status_{LOIO,LOFO}.tsv`, `jstatus_*.tsv` |
| 4 chance-match null | `null_shift.py` | 1 node, 8 threads | `null_shift.tsv.gz` |
| 5 locus depth + catalogue hits | `launch_depth.sh` -> `run_node_depth.sh` -> `sample_depth_catalogues.sh` | fleet (16 streams; the NAS is the limit) | `persample/<c>/<s>/{depth.tsv,cats.tsv}` |
| 6 scoring | `run_scores.sh` -> `sample_score.sh` (+ `make_params.py`) | fleet, `SHARD=i/k` | `persample/<c>/<s>/scores_{LOIO,LOFO}.tsv` |
| 7 collect | `collect_scores.py` | 1 node | `loci_scores.tsv.gz` |
| 8 sweep + selection | `sweep.py` | anywhere (pandas) | `sweep/` (curve, per-individual table, chosen parameters, stability) |
| 9 figures | `figures.py` | anywhere | `sweep/fig_*.png` |
| 10 shortcut check | `validate_rebuild.sh` (+ `compare_rebuild.py`) | 1 node per family, ~30 GB | `rebuild/<family>/` |

How leave-out is done without rebuilding the panel: `gp_common.ExcludingIndex`
wraps the released `panel_junction_carriers.tsv` index and drops the left-out
individuals from every matched site's carriers; `cxsv_matching.classify_locus`
then runs unchanged. `validate_rebuild.sh` compares this shortcut with a true
`build_panel_sites.py` rebuild without the family.

`leaveout_status.py --modes NONE` reproduces `classify_patient_cxsv.py`
(locus tables identical column for column on 4 panel samples; the two
junction-level columns the scorer does not read, `panel_callable` and
`panel_af`, are written as NA / '.').

Depth: by default (`DEPTH_LOCI=nonpublic`) read depth is measured only for the
loci that are not public in either leave-out mode, the only loci whose points a
privacy rule can change; public loci then score without copy-number features.
Measuring every locus reads nearly every CRAM container of a long-read
alignment (~80 TB over the panel).

## Where the results are

The summary tables of the item-10 calibration (`results/`: chosen_parameters.json,
strata_before_after.tsv, sweep_curve.tsv.gz, stability.tsv, …) are ignored by
git (`.gitignore`: `results/`). A copy with SHA256SUMS is on the NAS at
`/mnt/sargasso1/CxSV/calibration/item10/graded_privacy_results/` (copied
2026-11-04 from the item-10 worktree); the bulk data are in
`/mnt/sargasso1/CxSV/calibration/item10/`.
