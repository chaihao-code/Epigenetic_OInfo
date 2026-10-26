# Original developmental synergistic triplets

Main figure: `plots/mark_occurrence_per_tissue_by_stage_combined.pdf` (also PNG).
Preimplantation and E11.5–P0 are shown side by side. Arial 8 pt, titles 10 pt, 18 × 6 cm. Separate period figures are 10 × 6 cm. Palette: DNA methylation/H3K9me3/H3K27me3 cool colors; ATAC/H3K4me1/H3K4me3/H3K27ac warm colors; H3K36me3 gray.

Eight shared signals; five bins; 10,000 permutations; BH correction over all 56 triplets per sample × annotation block; synergistic calls require observed O-information < 0 and FDR <= 0.05.
Count each mark's occurrence across 40 annotation blocks within a sample. Early stages each have one sample. Late stages average the same five common tissues (Forebrain, Midbrain, Hindbrain, Heart, Liver). Supporting heatmap values are occurrence / 840.

Reproduce from preserved triplet-test inputs, then redraw the occurrence figures:

```bash
/lustre/home/2101110354/anaconda3/envs/hoi_env/bin/python scripts/plot_original_development_heatmap.py
/lustre/home/2101110354/anaconda3/envs/hoi_env/bin/python scripts/plot_mark_occurrence_per_tissue_by_stage.py
```

Run from this directory, or invoke the scripts by absolute path. Scripts resolve data relative to their own location.
`inputs/source_manifest.json` records provenance. The additional four early stages have a preserved BH56 test snapshot in `inputs/`. Original 2cell/8cell/ICM and late tests are read from the retained `../batch_effect_figures/data/`, ensuring exact agreement with the previous original bars. Upstream permutation results remain in mouse_preimp and mouse_after_E11.5.
`data/` contains plotting tables and validation records. No dependency on the deleted ABC or four/triple replacement directories remains.

Stacked A/B/C count heatmaps for all 40 annotation blocks:
`plots/triplet_ABC_stacked_post_E11p5_n5.{png,pdf}`. Reproduce with
`python scripts/plot_triplet_ABC_stacked_post_E11p5.py`. These use the same
eight signals, E11.5–P0 stages, five tissues, BH56 calls and A/B/C definitions
as the final trajectories below, including H3K36me3 in the H3K9me3/H3K27me3
repressive group. Cells show mean counts, with independent Reds scales per
class to match `mouse_development/synergistic_triplets/triplet_ABC_stacked.png`.
The canvas is exactly 18 × 12.3 cm. All text is Arial 8 pt except the cell
numbers, which are 4 pt to fit the 40 columns. Fixed canvas export preserves
the requested width; validation checks label bounds, stage-label overlap
and cell-number fit.
The matching `data/triplet_ABC_stacked_post_E11p5_n5_*` files preserve
per-sample counts, stage means and validation against the trajectory and
mark-occurrence tables.

Preimplantation companion heatmaps are saved beside the late figure as
`plots/triplet_ABC_stacked_preimplantation.{png,pdf}`. Reproduce with
`python scripts/plot_triplet_ABC_stacked_preimplantation.py`. The seven rows
are zygote, 2-cell, 4-cell, 8-cell, morula, ICM and TE, with one original sample
per stage and ICM/TE shown separately. Signals, A/B/C definitions, 40 annotation
blocks and BH56 selection are identical to the late figure. The retained mask-0
tests supply 2cell/8cell/ICM; the preserved BH56 snapshot supplies the other four
stages, matching the existing developmental mark-occurrence analysis.
The figure shares the late heatmap's 18 × 12.3 cm layout, Arial 8 pt labels and
4 pt cell numbers through a shared renderer. Cells show sample counts;
each class and period uses its own color range. Supporting sample counts,
summary tables and source/validation metadata are in `data/` with the same
prefix. Validation independently recomputes all BH56 values, checks the eight
signals and H3K36me3 classification, and matches every per-mark occurrence
count in all seven samples to the existing analysis.

Selected-region overlaid companion to the seven-signal developmental figure:
`plots/triplet_ABC_selected_region_trajectories_8signals_post_E11p5_n5.{png,pdf}`.
Run `python scripts/plot_triplet_ABC_selected_regions_post_E11p5.py`.
Six panels overlay A/B/C in two rows of three. The top row is exon/SINE,
distal enhancer CTCF−/SINE and proximal enhancer CTCF+/CpG island; the bottom
row is intron/SINE, distal enhancer CTCF−/LINE and None/None (`None_AND_None`).
The middle column compares SINE/LINE within distal CTCF− enhancers. None/None
has zero significant A/B/C triplets throughout E11.5–P0 and is labeled
explicitly to distinguish coincident zero curves from missing data.
The region order, colors, y limits and six-panel layout
match the seven-signal selected-region figure. Both figures are 18 cm wide
and approximately 9.67 cm high, with all text in Arial 8 pt.
Both share x axes between rows and show a single lane of stage labels only
on the bottom row; the seven-signal figure retains all 17 stage labels.
Both the canvas and each panel are compressed to 3/4 of the preceding version's
height, with a panel height/width ratio of approximately 0.5606 and matching
physical axes sizes. All region titles are on one line, and the two rows share
one y-axis title. Data retain the eight-signal definitions (including
H3K36me3 alongside H3K9me3/H3K27me3), A/B/C denominators 15/31/10,
BH56 q <= 0.05, O-information < 0, and the same five common tissues at
every E11.5–P0 stage. A contains DNA methylation and at least one of the
three repressive marks; B contains a repressive mark without DNA methylation;
C contains none of those repressive marks. SEM uses population SD / sqrt(5), as in the original
trajectory plot. The seven-signal companion uses all available tissues;
these figures therefore do not isolate the effect of adding H3K36me3.
Supporting `_by_sample.tsv`, `_summary.tsv` and `_metadata.json` under
`data/` share the figure prefix. Validation checks complete BH56 families,
H3K36me3 classification, agreement with all-region counts and label bounds.

Superseded 3 × 4 ABC trajectory figures and their derived tables/validation
records have been removed at the user's request, along with exploratory
region-comparison outputs. Final selected-region figures and their supporting
tables, heatmaps, input data and source scripts are retained. The exact deleted
paths, sizes and SHA256 hashes are recorded in `data/trajectory_cleanup_manifest.json`.
