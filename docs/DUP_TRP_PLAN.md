# Plan: sensitive detection of large duplication–triplication events (DUP-TRP/INV-DUP)

Status: **proposed** (2026-11-05). Nothing is implemented yet. This plan needs PI approval before Phase 1 starts.

## 1. The target

The literature set lit-v1.0 contains 26 published DUP-TRP/INV-DUP events (23 MECP2, 3 PLP1; Grochowski 2024, Bahrambeigi 2019).

| Property | Value |
|---|---|
| Location and carriers | All on chrX, all in male carriers, so CN on X goes 1 → 2 → 3 → 2 → 1. In a female or on an autosome it would go 2 → 3 → 4 → 3 → 2. |
| Total span on hs1 | 0.42–7.4 Mb, median 0.60 Mb |
| Junctions | 2 per event (4 breakends) |
| Junction 1 | Inverted junction between inverted low-copy repeats. 19/26 share the breakend chrX:152,631,178 (Xq28 LCR pair). |
| Junction 2 | Non-recurrent (template switch) |
| Triplicated segment | 131–317 kb where curated (n = 2); 63 kb–1.6 Mb for the duplicated flanks |

CxSV calls **none** of the 26 under any long-span rule. The positive-control add-on L01 (PLP1 DUP-TRP/INV-DUP) is also not called.

## 2. Why the current pipeline misses them

1. **Junction span.** The junction records span the duplicated segment. Most are > 500 kb, so the item-13 long-span rule excludes them from discovery (L01's 532.8 kb junction).
2. **Clustering distance.** The two junctions are usually > 50 kb apart. Distance clustering therefore splits them into two clusters, each holding a single INV, and neither cluster is "complex".
3. **No copy number in discovery.** The defining feature of the event, a nested copy-number gain, is never used to link records. Depth enters only at scoring, per locus.
4. **Repeat context.** Junction 1 sits in inverted segmental duplications. Callers may miss it or call it only once, and it may match panel junction sites in the segdup.

Mechanisms 1–3 also cause most of the 55 literature events lost under the exclude rule: 28 because their long junctions are excluded, 26 because they are split into non-complex clusters. So the fix proposed here is a general **copy-number-anchored** discovery path, not a DUP-TRP special case.

## 3. Approach: copy-number-anchored linkage (a second discovery track)

The new track runs beside the existing junction-clustering track and adds loci. It removes nothing.

1. **Binned depth per genome.** Run mosdepth `--by` 10 kb bins at MAPQ ≥ 20, plus 1 kb bins inside candidate regions for refinement, on every panel CRAM and on the patient. The CRAM 3.1 → 3.0 transcode path from `callable_mask.sh` is reused. The panel currently keeps only the mosdepth summary and the callable mask, so this is a new pass over the 1,619 alignments.
2. **Panel-of-normals normalisation.** Build it per platform stratum (ONT-1019, LRSC R9, LRSC R10, PacBio) and per sex for X and Y. For each bin, store the median and MAD of normalised depth across the healthy panel. This removes GC, mappability and platform bias empirically. Bins with high panel variance (segdups, satellites, common CNV) are down-weighted. The 1,619-genome panel is the main asset here.
3. **Segmentation.** A ploidy-aware HMM over CN states 0–6 (or CBS), with segments of at least 3 bins. Outputs are CN segments with a confidence each.
4. **Candidate structures.**
   - A gain containing a higher gain (CN +1 with a nested +2): the DUP-TRP signature.
   - Adjacent gains or losses at different levels.
   - Any gain or loss whose boundaries coincide with SV breakends.
5. **Distance-free junction association.** Attach SV records of any length to a candidate when their breakends fall at a segment boundary (tolerance set by bin size, then refined with split reads) or inside the segment. This includes the records now set aside in `long_span_records.tsv` and `single_caller_filtered.tsv`.
   - A long record is readmitted only when copy number confirms it.
   - Item 13 found that 73–77% of long DEL/DUP records in healthy genomes have no depth support, so those stay excluded.
6. **Structure model.** Classify as DUP-TRP/INV-DUP when two conditions hold:
   - the nested-gain pattern is present;
   - an inverted junction (INV/BND ++ or --) is consistent with the Carvalho model, using segdup orientation for junction 1.

   Otherwise the locus becomes a "CN-linked complex gain/loss" (for example DUP-NML-DUP, DEL-NML-DEL or DUP-INV-DUP).
7. **Scoring.** Use the existing scorer, with no weight changes in Phase 1. Dosage comes from the CN segments:
   - Triplosensitivity uses ClinGen genes and regions (item 9).
   - Privacy uses the existing graded junction-level f, plus a segment-level CN frequency in the panel: a gain seen in healthy genomes is public.
8. **Read-level confirmation (Phase 2, optional).** Look for reads spanning the inverted junction, and ultra-long reads that cross the LCR.

## 4. Validation and success criteria

| Test | Data | Pass criterion (proposed) |
|---|---|---|
| A. Literature, synthetic | 26 DUP-TRP events with CN modelled from the published structure; also the other 37 non-DUP-TRP junction-level events | ≥ 80% of DUP-TRP called; report the gain on the other 37 |
| B. Read-level simulation grid | Planted DUP-TRP/INV-DUP at Xq28 (MECP2) and PLP1 in male hosts, plus autosomal cases in a female host. Total span 0.2–2 Mb, triplicated segment 50–500 kb. Coverage 10, 20 and 30x on ONT R10, R9 and PacBio. Uses the existing positive-control simulation model. | Sensitivity curve by size × coverage × platform; ≥ 90% at ≥ 20x for triplicated segments ≥ 100 kb |
| C. Full-genome positive control | E01–E13 + L01–L04 from the bundle | L01 called at Tier 1; every other event unchanged |
| D. Healthy specificity | Full-panel LOFO, 1,619 genomes | Healthy Tier 1 rises by no more than 0.05 per genome (from 0.251); stop rule set by the PI |
| E. Stand-in and host | HG01868, host 3559 | 0 Tier 1 |

## 5. Phasing

1. **Phase 1 — prototype outside the pipeline** (`calibration/cn_linkage/`).
   - Depth on a 100-genome stratified subset (the item-13 evaluation set), the host and HG01868.
   - Panel-of-normals from those 100 genomes.
   - Run tests A, B and C, plus D on the subset.
   - Measure depth-pass cost per genome before estimating the full-panel pass. Reading every CRAM again over the NAS is the main cost, and athena's 1 GbE link is the bottleneck.
   - Decision point for the PI.
2. **Phase 2 — integrate.**
   - Add a per-sample `cn_depth` rule and a `cn_linkage` discovery rule.
   - Ship the panel-of-normals as a release component.
   - Run the full-panel depth pass on the fleet (freya, shadowfax, keurig, athena).
   - Then a downstream rebuild as panel-v1.2, full LOFO, bundle, validation, PC and tag.
   - This changes the calling fingerprint. Nothing else upstream is recomputed: callers, merge and callability are reused.

## 6. Decisions for the PI before Phase 1

1. **Scope.** DUP-TRP only, or all CN-linked complex events from the same machinery? The same machinery would also address most of the 55 missed literature events. Recommended: all, with DUP-TRP as the primary acceptance target.
2. **Specificity stop rule.** Proposed: a healthy Tier 1 increase of at most +0.05 per genome.
3. **Tooling.** In-house (mosdepth bins + panel-of-normals + HMM, no new container dependency), or adding a long-read CNV caller such as Spectre (ONT) or HiFiCNV (PacBio) to the image. That would change the image and need validation. Recommended: in-house for Phase 1, with one external caller run outside the pipeline as a comparator.

## 7. Risks

- **Depth noise.** Small triplicated segments at < 15x will be hard to call; the coverage floor decision means such genomes exist in the panel.
- **Segdup repeats.** Mapping in inverted LCRs is poor, so junction 1 may be unseen. The structure model must allow a CN-only call, which would then score lower.
- **Chimeric libraries.** They add false inverted junctions, so the chimerism QC (in progress) should be in place first.
- **Compute.** A full depth pass over about 1,600 CRAMs is the largest new compute since the panel build.
