# Shared-eight-signal triplet nulls

66 later stage x tissue samples from mouse_after_E11.5.
Inputs are the original ten-signal lists, excluding H3K4me2 and H3K9ac.
The eight signals are ATAC, DNA methylation, H3K4me1, H3K4me3, H3K27ac,
H3K9me3, H3K27me3 and H3K36me3.
5 bins, max 100000 loci/block, sampling seed 42; 10000 within-block independent
permutations, permutation seed 1234, null orders <= 3, 16 CPUs/sample.
The existing ten-signal full-system outputs remain intact.
The downstream analysis applies BH56 per sample x block.
Current developmental ABC plots/statistics are in project/H3K36me3_developmental_switch/synergistic_triplets/;
replacement-source robustness summaries are in its batch_effect_figures/ directory.
