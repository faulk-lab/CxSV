# Positive control: sample 3559 with a catalogue of planted CxSVs

The lab positive control is a real long-read genome (sample 3559, ONT R10.4.1,
~39x) into which twelve heterozygous structural rearrangements have been
planted at known coordinates: ten complex SVs of different mechanisms and two
simple SVs that must **not** be called complex. It tests the whole patient path
on real background sequence and real host variation: calling, CxSV clustering
and classification, panel matching (a planted event must come out *private*)
and the private VCF.

> **Data handling.** Sample 3559 is human-subject data (lab-only). The control
> CRAM, the host CRAM and anything else containing its reads stay on the lab
> NAS. Only the generator code, `events.yaml`, the truth set and summary
> counts are in the repository. The control is never used in CI; CI runs the
> fully synthetic test in `tests/synthetic/`.

![Planted events](figures/pos_control_events.png)

*Figure: structure of each planted event, drawn from `events.yaml` by
`tests/positive_control/plot_events.py`. Top track per panel: reference,
segments to scale, deleted segments hatched. Lower track(s): derivative
haplotype(s), arrows pointing left are reverse-complemented. A segment keeps its
letter and colour in every track; red triangles mark novel junctions (the truth
set). Flanks are shortened; the two inter-chromosomal events are drawn
schematically.*

## The catalogue

Coordinates are T2T-CHM13v2 (UCSC hs1 names), 0-based half-open. All events
are heterozygous (host is female, so the chrX event at *MECP2* is
heterozygous as well). Path notation: segments 5'->3' on the derivative
haplotype, `+` reference orientation, `-` reverse complement, `p`/`q` flanks.

| ID | Event | Gene(s) | Span (hs1) | Derivative path | Expected class |
|---|---|---|---|---|---|
| E01 | DUP-TRP/INV-DUP | *MECP2* | chrX:152,320,014-152,349,014 | `der: p+ B+ C+ D+ C- B+ C+ D+ q+` | DUP+INV |
| E02 | delINVdel | *NF1* | chr17:32,243,849-32,271,849 | `der: p+ I- q+` | DEL+INV |
| E03 | dupINVdup | *PMP22* | chr17:15,140,559-15,166,559 | `der: p+ A+ C- B- A- C+ q+` | DUP+INV |
| E04 | DEL+INV+DUP | *SHANK3* | chr22:51,207,890-51,230,915 | `der: p+ A+ C- D+ E+ E+ q+` | DEL+DUP+INV |
| E05 | nested inversion with inner deletion | *NRXN1* | chr2:50,453,329-50,488,904 | `der: p+ E- C+ B- A- q+` | DEL+INV |
| E06 | templated insertions at a deletion (FoSTeS/MMBIR) | *AUTS2* | chr7:71,340,258-71,406,158 | `der: p+ T1+ U1+ T1+ T2- U2+ T2+ q+` | INS+DEL |
| E07 | reciprocal translocation with junction deletions | *SATB2*/*FOXP2* | chr2:199,854,293-199,856,293; chr7:115,868,955-115,869,905 | `der(2): p_chr2+ q_chr7+ / der(7): p_chr7+ q_chr2+` | BND_cluster |
| E08 | chromoplexy (balanced 3-chromosome chain) | *NSD2*/*EHMT1*/*ZIC2* | chr4:1,931,409-1,931,784; chr9:149,997,337-149,997,837; chr13:99,195,184-99,195,684 | `der(4): p_chr4+ q_chr9+ / der(9): p_chr9+ q_chr13+ / der(13): p_chr13+ q_chr4+` | BND_cluster |
| E09 | localized chromothripsis | *CHD2* | chr15:90,515,727-90,925,728 | `der: p+ S1+ S5- S3+ S9- S7+ S2+ S11- S6+ S10- q+` | Chromothripsis_like |
| E10 | multi-breakpoint cluster | intergenic | chr5:25,003,550-25,041,800 | `der: p+ A+ C- D+ E- G+ q+` | Multi_breakpoint_cluster |
| E11 | simple deletion (specificity control) | intergenic | chr10:29,098,000-29,115,650 | `der: p+ q+` | not a CxSV |
| E12 | simple inversion (specificity control) | intergenic | chr12:74,754,175-74,814,325 | `der: p+ I- q+` | not a CxSV |

Placement rules (`design_events.py`, seed 3559): each catalogue template is
centred on its gene (gene events must disrupt an exon); every breakpoint is
moved by up to +/-5 kb to sequence with no segmental duplication within 10 kb
and <=30% soft-masked bases within +/-250 bp, so that the control measures the
pipeline rather than mappability. 37 novel junctions in total.

## How it is built

All stages are scripts in `tests/positive_control/`; outputs live under
`/mnt/sargasso1/CxSV/datasets/pos_control/`.

1. **Host realignment.** The original unaligned modBAM
   (`/mnt/pacific/monkey_meth/Human_raw_data/modbams/BloodofChildren_3559.allC.m6A.mod.bam`;
   dorado 0.8.3, `dna_r10.4.1_e8.2_400bps_sup@v5.0.0`) is realigned exactly like
   the background panel: `samtools fastq -T MM,ML,MN,qs | minimap2 2.24-r1122
   -y -x map-ont -a --eqx -K 4g` (k17 hs1 index), RG `ID:3559 SM:pos_control_3559`.
   13,928,513 reads in = 13,928,513 primary alignments out. -> `host/`
2. **Event design** - `design_events.py` -> `events.yaml`.
3. **Haplotypes and truth** - `build_haplotypes.py events.yaml hs1.fa build --check`
   writes, per event, `alt.fa` (derivative haplotype(s) with 150 kb flanks),
   `ref.fa` (the matching reference windows) and `replace.bed` (event span
   +/- 50 kb), plus `truth/` (`truth_junctions.tsv`, `truth.vcf` as VCF 4.2
   breakends with `GT 0/1`, `truth_events.bed`, `copy_number.tsv`). -> `build/`
4. **Read model** - Badread 0.4.2 error and q-score models trained on 20,000
   host reads aligned with `minimap2 -c`; read lengths and identities fitted to
   the host (200,000 reads: mean 7,110 bp, SD 8,683, N50 12,836; identity mean
   98.09%, SD 2.96 over 17,236 MAPQ>=20 alignments). The beta maximum is capped
   at 99.8% because the host 99.5th-percentile identity is 100%. Chimeras, junk
   and random reads, glitches and adapters are off. Check run: simulated reads
   had mean length 7,112 bp, identity 97.98% (SD 2.97) and mean read Q 39.5
   (host 39.1). -> `sim_model/`
5. **Simulate and splice** - `simulate_and_splice.py`: per event, local depth
   = median host depth over the replace windows; the derivative and the
   reference windows are each simulated at half that depth (one copy each).
   Simulated reads are aligned like the host. Every host read with any alignment
   in a replace window is removed (all its records); every simulated read with
   any alignment in a replace window is added (all its records). Simulated
   reads carry `RG ID:sim` with the host's `SM`, are named
   `sim_<event>_<alt|ref>_<n>` (origins in `sim_read_origins.tsv`) and have
   **no MM/ML methylation tags**. `check_splice.py` compares host and spliced
   depth in 2 kb bins across each window. -> `release/`

## Running the regression

```bash
conda activate cxsv-pipeline
bash tests/lab/pos_control_regression.sh ~/scratch/cxsv_posctrl_regression 10
RERUN=1 bash tests/lab/pos_control_regression.sh ...   # plus determinism check
```

The script slices every replace window +/- 1 Mb from the control and from 10
background genomes (registry sheet, >=17x), builds a panel from the
background slices, runs the patient workflow on the control slice and scores
it with `score_pos_control.py` (junctions matched within +/-1 kb; per-event
verdict against `expect` in `events.yaml`). Exit 0 means every event meets its
expectation.

## Splice checks (2026-10-28)

- Host: 13,928,513 reads in the original modBAM = 13,928,513 primary
  alignments after realignment; ~39x on autosomes, chrX 39.7x, chrY 0.92x
  (female; chrX/autosome mean-depth ratio 1.03).
- Simulated: local depth 32-43x per event, split 1:1 between the derivative
  and the reference copy (`release/truth/simulation_plan.tsv`).
- Spliced: 11,176 host reads removed, 11,855 simulated reads added; final
  primaries 13,929,192 (= 13,928,513 - 11,176 + 11,855).
- Depth (2 kb bins, windows +/- 25 kb, `check_splice.py`): spliced/host ratio
  median 0.99 over 1,501 bins, 0.96 over the 29 window-edge bins; 5th-95th
  percentile 0.62-1.39, the spread expected for two independent samplings at
  ~38x with 7 kb reads. Simulated depth is flat at the window median, so the
  host's own local depth variation does not continue inside a window.
- Release (`/mnt/sargasso1/CxSV/datasets/pos_control/release/`):
  `pos_control_3559.hs1.minimap2.cram` (md5 77d1b04f9c27fcd9d6b4b25a469743eb),
  `truth/`, `MD5SUMS`.

## Latest regression result (2026-10-28, 10 background slices)

Patient workflow on the control slice (15 merged regions, 32.2 Mb) against a
panel of 10 ONT background genomes (>=17x) sliced the same way. A second run
gave identical `merged_sv.vcf`, `cxsv_panel_status.tsv` and `private_cxsv.vcf`
records.

| ID | Junctions exact/partial/missed | Pipeline class | Expected class | Panel status (carriers) | Verdict |
|---|---|---|---|---|---|
| E01 | 2/0/0 | Chromothripsis_like | DUP+INV | private (0) | PASS |
| E02 | 2/0/0 | Chromothripsis_like | DEL+INV | private (0) | PASS |
| E03 | 2/0/0 | Chromothripsis_like | DUP+INV | private (0) | PASS |
| E04 | 3/0/0 | Chromothripsis_like | DEL+DUP+INV | **public (4)** | known failure |
| E05 | 4/0/0 | Chromothripsis_like | DEL+INV | private (0) | PASS |
| E06 | 3/0/0 | DUP+INV | INS+DEL | private (0) | PASS |
| E07 | 2/0/0 | Complex_other (2 loci) | BND_cluster | private (0) | PASS |
| E08 | 2/1/0 | Complex_other | BND_cluster | private (0) | PASS |
| E09 | 9/0/0 | Chromothripsis_like | Chromothripsis_like | private (0) | PASS |
| E10 | 4/0/0 | Chromothripsis_like | Multi_breakpoint_cluster | private (0) | PASS |
| E11 | 1/0/0 | not a CxSV | not a CxSV | - | PASS |
| E12 | 2/0/0 | **Nested_SV (CxSV)** | not a CxSV | private (0) | known failure |

Summary: 36 of 37 junctions recovered exactly, 1 partially, none missed. The
partial one is E08's chr4-chr9 junction: no chr4-chr9 call was made; it scores
partial only because the chr13 junction calls end within the 1 kb slop of it
(the chromoplexy junction deletions are 375-500 bp).
Ten of twelve events meet their expectation. Outside the events, the slices
hold 17 CxSV loci of the host's own, 16 of them public in the 10-genome panel.

Findings (listed in `tests/positive_control/known_failures.tsv`; the
regression fails on any new failure, or if a listed one starts passing):

1. **E04 is labelled public.** The planted event is absent from the panel, but
   its locus overlaps a different panel CxSV at the *SHANK3* subtelomeric
   hotspot (chr22:51,150,520-51,229,698, 5 of 10 background genomes). This is
   the masking failure that roadmap item 3 (breakpoint-level matching) targets.
2. **A simple inversion is called a CxSV (E12).** Sniffles2 reports both
   reciprocal breakends of one inversion; SURVIVOR keeps them as two INV
   records with identical coordinates, and C2 (nested SV) fires on a record
   "nested" in its own mate. Every balanced inversion is exposed to this.
   Fix: collapse records with the same breakpoint pair before C1-C5.
3. **Class labels.** `annotate_cxsv.py` tests C5 (>=3 SVs in 10-50 kb) before the
   type-based labels, so compact events with 2-4 junctions are labelled
   `Chromothripsis_like`; the label matched the expected class for 1 of 10
   complex events (E09). The README class table also omits the `Nested_SV`,
   `Overlapping_types` and `CN_plus_orientation_change` labels the code emits.
   Belongs with roadmap item 5.
4. **Caller support.** Of 37 junctions, Sniffles2 recovers 32, cuteSV 19 and
   SVIM 15; 11 junctions are carried by a single caller. Raising
   `survivor_params.min_callers` above 1 would lose them. The inter-chromosomal
   junctions (E07, E08) come mostly from SVIM/cuteSV records.

Score tables: `score.events.tsv`, `score.junctions.tsv` in the regression
work directory; the off-target table (`score.offtarget.tsv`) lists host loci
and stays on the NAS.

## Replaces

The previous control (`/mnt/sargasso2/CxSV/synthetic_pos_control/BloodofChildren_3559.T2T.nested-inv.bam`,
RefSeq contig names) was meant to carry a nested inversion but holds a single
~100 kb homozygous inversion at chr1:10,750,000-10,850,000: VISOR/LASeR
simulated only one haplotype, replaced the real reads of chr1:10-11 Mb, and
gave ~2x depth there. It is superseded by this control.
