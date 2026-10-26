# Four-mark factorial replacement: FDR participation and stability

This directory now contains only the computation and reports needed for
same_stage_bulk_replacement_all_marks and the requested 16-way stability analysis.
Other previous local analyses were cleared. Upstream source projects are retained.
Source replacement assesses sensitivity to track provenance; it does not perform
model-based batch-effect removal or identify a causal batch effect.

## Design and computation

2cell, 8cell and ICM; 16 subsets of four replacement marks per stage (48 runs).
Bit order: bit0 H3K36me3, bit1 H3K9me3, bit2 H3K27me3, bit3 H3K4me3.
The 16 include the unchanged original (mask 0) and 15 nonempty replacements.
Replacement sources follow project/bulk_factorial_triplets/manifest.json:
H3K36me3 GSE112834; H3K9me3 GSE97778 replicate mean;
H3K27me3 and H3K4me3 GSE73952 replicate means. All other tracks are fixed.
Three missing historical H3K36me3 aliases are rebuilt as local symlinks to the
original candidate data files recorded in project/batch_effect/metadata/candidates.json.

All runs use the same 8 signals (exclude H2A.Z), 5 bins, 40 mm10 region strata,
max 100000 loci per stratum, sampling seed 42, permutation seed 1234,
and 10000 independent within-block permutations. Null orders are capped at 3.
Lower-tail P = (1 + number(null <= observed)) / 10001.
BH correction is across all 56 triplets within each variant x region, before
filtering by sign. Significant synergy: O-information < 0 and BH q <= 0.05.
No adjustment is performed across variants or stages; the BH family is explicit.

## Participation and late reference

For each mark and region, participation = significant triplets containing that
mark / 21 possible triplets. The early bars average the 40 region fractions.
All 16 variants and the late reference appear in each of the 3 stage panels.
Late reference: current 10000-permutation shared-eight-signal results from
project/mouse_after_E11.5/synergistic_triplets; fresh BH56 is applied identically.
Reference tissues: Forebrain, Midbrain, Hindbrain, Heart and Liver, present at
all 7 later stages (E11.5-P0). First average 40 regions within a tissue-stage,
then stages within a tissue, then average the 5 tissues equally. Error bars
are sample SD across the five tissue means, not confidence intervals.

## Stability

S_variant is the set of significant (genomic region, canonical triplet) keys.
For each early stage, strict stability = |intersection of all 16 S_variant| /
|union of all 16 S_variant|. Original mask 0 is included. The 40 regions retain
their identities. An empty union is undefined (blank/NA), not perfect stability.
Also reported per region and per participating mark. This is the multi-set
intersection/union ratio, not the average pairwise Jaccard index.
BH thresholds can change between variants: unchanged raw triplet statistics
need not retain the same corrected significance when other triplets change.
An unaffected-triplet check confirms observed O-information is unchanged for
triplets not containing any replaced mark (tolerance 1e-9).

## Programs, outputs and execution

scripts/compute_bulk16.slurm: 48 tasks x 16 CPUs, 24 h, no concurrency cap.
scripts/plot_batch_effect.py: BH56 tables, numeric stability summaries and the single primary participation figure.
scripts/plot.slurm: scheduled after successful completion of all 48 runs.
Full input provenance and parameters: manifest_16_fdr.json; submission IDs: submission.json.

plots/same_stage_bulk_replacement_all_marks.{png,pdf}: the only retained display figure.
Width 18 cm; each of the three plotting rectangles is exactly 3 cm high.
Total height 15 cm includes stage titles, rotated signal labels and bottom legend.
All text is Arial 8 pt. Four-column legend and brief method caption appear below the panels.
The replacement bars use light, moderately saturated colors; Original is gray and Late black.
Other batch figure products have been removed; their numeric analyses are retained.
Redraw validated summary tables without recomputing FDR: python scripts/plot_batch_effect.py --plots-only.
stability_by_stage.tsv: 3 primary ratios and intersection/union sizes.
stability_by_stage_region.tsv and stability_by_stage_mark.tsv: stratified ratios.
same_stage_mark_summary.tsv: all16 early participation fractions.
data/: full BH56 tests, per-region participation, 16-way support frequencies,
late-reference tests/tissue means and unaffected-triplet numerical validation.
validation.json is written only after all calculations and figures complete.

Historical support-frequency and sensitivity tables remain in data/ and the threshold-specific directories.

Original developmental heatmap and occurrence scripts, data and figures have moved to
`../synergistic_triplets/`. See that directory for the current original-only developmental figures.
