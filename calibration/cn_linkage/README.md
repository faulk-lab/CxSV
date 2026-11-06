# calibration/cn_linkage — copy-number-anchored linkage, Phase 1 prototype (roadmap item 15)

Report: [docs/CN_LINKAGE_PHASE1.md](../../docs/CN_LINKAGE_PHASE1.md). Plan: [docs/DUP_TRP_PLAN.md](../../docs/DUP_TRP_PLAN.md).
Nothing here is used by the workflows; it runs in the pipeline image beside a panel release.

| Step | Script | What it does |
|---|---|---|
| depth | `depth_pass.sh`, `depth_worker.sh` | mosdepth 1 kb bins (MAPQ ≥ 20, fast mode; CRAM 3.1 transcoded), timing per genome; fleet claims |
| genomes | `make_genomes.py`, `subset100.tsv` | genomes.tsv for the subset + patients; sex checked on chrX depth |
| segmentation | `cnlib.py`, `segment_all.py` | PoN per stratum (unflagged genomes), ploidy-aware HMM, 1 kb boundary refinement, nested-level scan |
| linkage | `cnlink.py`, `augment.py`, `link_all.py` | candidates, distance-free record association, structure classes, panel CN frequency, loci added to the discovery tree |
| scoring | `patient_rescore.py` | privacy / catalogues / depth for the new loci, deployed scorer |
| chain | `run_phase1.sh` | make_genomes → segment_all → link_all → subset manifest → `cn_lofo.sh` (test D) |
| test A | `literature_cn.py` | literature events with copy number modelled on a healthy donor's depth |
| test B | `grid.py`, `grid_worker.sh`, `train_model.sh` | read-level simulation grid (design, simulate, run, collect) on the fleet |
| test C/E | `smoke_pc.py`, `patient_rescore.py` | positive control (truth-CN smoke test; real depth), HG01868 |
| test D | `cn_lofo.sh`, `testD_compare.py` | LOFO on the augmented subset vs panel-v1.1-r2 |
| comparator | `spectre_compare.sh` | Spectre on the same mosdepth tracks |
| figures | `make_figures.py` | the three report figures |

Unit tests: `tests/unit/test_cn_linkage.py` (numpy only).
