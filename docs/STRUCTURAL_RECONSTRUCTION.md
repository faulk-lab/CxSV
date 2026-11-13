# Structural reconstruction of Tier 1 loci (clinician report)

PI request 2026-11-13: for every Tier 1 hit the clinician report shows a graphical
reconstruction of the variant against the reference, in the style of the positive-control
event figure (`docs/figures/pos_control_events.png`).

Code (reporting only; **not** in the calling fingerprint, `cxsv_fingerprint.CALLING_CODE`):

| file | role |
|---|---|
| `workflow/lib/cxsv_recon.py` | segments, junctions, copy number, derivative walk, ambiguity |
| `workflow/lib/cxsv_recon_figure.py` | drawing primitives shared with `tests/positive_control/plot_events.py` (the PC figure renders pixel-identically from them) and the report panel |
| `workflow/scripts/clinician_report.py` | section "Structural reconstruction of the Tier 1 loci" after the Tier 1 points graphic |
| `workflow/rules/clinician_report.smk` | new inputs (below) and outputs `report/<sample>.reconstruction/NN_<locus>.png`, `report/<sample>.tier1_reconstruction.tsv` |
| `tests/positive_control/recon_vs_truth.py` | validation against the positive-control truth |
| `tests/unit/test_recon.py` | synthetic and truth-derived unit tests |

## Inputs (this run's outputs only; works offline in the bundle, no CRAM access)

* `junctions/junctions.tsv`: the locus's junctions; distinct junctions as in scoring
  (`cxsv_classes.distinct_junctions`). Composite records (SVIM's averaged inversion) are not
  breakpoints; one is used only if both its ends coincide with breakpoints of kept junctions.
* `callers/{sniffles,cutesv,svim}_filtered.vcf`: the join orientation of every caller record
  at the junction (BND ALT brackets; DEL `+-`, DUP `-+`; cuteSV INV `STRAND`; a symbolic INV
  without strand stands for both inversion junctions) and its read support (SUPPORT / RE).
  The merged call set is a fallback only: SURVIVOR keeps one record per merged group and types
  a BND pair on one chromosome as INV whatever its brackets say.
* `scoring/locus_depth.tsv`: per-segment copy number from the scorer's own locus depth
  (`cxsv_scoring.cn_segments`); `cn/segments.tsv` for segments that the locus depth does not
  measure; a CN boundary more than 20 kb from any junction is reported as a note.
* The genes BED of the bundle. It merges overlapping genes into one interval (for example
  `EVI2A,EVI2B,NF1,OMG`); such an interval is reported as one unit
  (`NF1/EVI2A/EVI2B/OMG`, part: x of y kb) and never split into invented extents.

## Method

1. Breakpoints = junction ends clustered within 500 bp (the caller-agreement tolerance).
   Segments are the intervals between consecutive breakpoints; a 20 kb flank on each side
   (drawn shortened). Several chromosomes: flanks `p_<chrom>` / `q_<chrom>`.
2. Heterozygous model: the derivative carries `CN - (ploidy - 1)` copies of a segment.
3. The derivative is a walk through the adjacency graph (reference adjacencies plus the
   novel junctions in their called orientation) from a flank to a flank, every flank used
   once. Cost: 10 per junction left unused, 5 per open end (the walk had to stop: a junction
   is probably missing; drawn joined at a hollow triangle), 5 x (|depth CN - path CN| - 0.25)
   per segment, 1 per copy away from one for a segment too short for depth.
4. Among equal-cost walks: more read support of the joins used; fewer joins that only a
   symbolic INV implies; no derivative joining two left (or two right) chromosome ends
   (dicentric/acentric when the breaks lie on the same side of both centromeres; the bundle
   has no centromere positions, so this is noted, not enforced); fewer inverted bp.
5. **Ambiguous** (stamped on the panel, reasons and alternatives in the text and TSV) when
   walks tie on cost and support (the same joins taken in another order: an adjacency-only
   reconstruction cannot order them; spanning reads would), a junction stays unused, an end
   stays open, depth disagrees with the path by more than 0.6 copies, callers disagree on a
   join (deletion/duplication vs inversion) where depth cannot decide, or the search limit
   is hit. Nothing is drawn that the junctions and depth do not show.
6. Caller type labels are not trusted over orientation: a junction called DEL by one caller
   and INV by another is resolved by the record orientation and depth (a loss supports a
   deletion; an inversion keeps the copy number), and the reasoning is printed.

## Report

One panel per Tier 1 locus, highest score first, at most `reconstruction_max_panels` (6;
the rest are in the TSV), each followed by a table: inferred structure and derivative,
segments (coordinates, CN and its source, fate on the derivative, genes), junctions (both
ends, join orientation, reads per caller, the callers' type labels), ambiguity and
alternatives, notes. Caption: *Inferred from split-read junctions and read depth in this
sample; not orthogonally confirmed.* No Tier 1 locus: no section.

Page rule: a report may take the pages of the same report without reconstructions (2 for
the accepted layout) plus one page per started group of `reconstruction_panels_per_page`
(2) panels: `clinician_report.page_limit()`. The script measures both and warns in the log
when over; the synthetic end-to-end check enforces it.

## Validation (positive control, panel-v1.3 bundle 718df79 run outputs)

`tests/positive_control/recon_vs_truth.py --run <pos_control run> --truth events.yaml --truth addons/L0*/addon.yaml`:
the 16 Tier 1 loci of the PC run, 14 of them on planted events with a truth derivative
(E01-E09, L01, L02, L04; E07 and L04 each give two loci; L03 has no truth derivative,
windowed chromothripsis, 2 loci skipped).

| outcome | loci | events |
|---|---|---|
| drawn derivative = truth (segment order and orientation) | 13 | E01, E02, E04, E05, E06, E07 (2), E08, E09, L01, L02, L04 (2) |
| truth among the equally consistent alternatives listed | 1 | E03 (dupINVdup: 3 orders use the same joins) |
| mismatch | 0 | |

Of the 13 exact loci, 4 are flagged ambiguous: E01 and L01 (DUP-TRP/INV-DUP, 4 orders
each; L01 also depth CN 3.9/4.6 vs 3/4 on the path), E06 (2 orders: U2 orientation) and
E08 (the third chromoplexy junction is in another locus: drawn across a gap). The
tie-breaking rules were written while developing on these events; there is no held-out set.
Truth-derived unit tests (junctions and depth generated from each events.yaml truth path):
12/13 exact, E03 in alternatives.
