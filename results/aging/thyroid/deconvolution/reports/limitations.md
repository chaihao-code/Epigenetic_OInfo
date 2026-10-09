# Limitations

1. **n = 6 donors.** All age-association analyses have six independent
   samples. Results emphasize effect sizes and confidence intervals; p-values
   are descriptive. Sex (5F/1M), disease status and batch cannot be adjusted
   simultaneously with age without overfitting and were not.

2. **No methylation coverage.** The WGBS methylation bigWigs are 0–100
   percent values without per-CpG coverage or read counts. UXM_deconv
   (fragment-level) cannot run. Deconvolution uses constrained NNLS on beta
   values (linear-mixture approximation). Consequently:
   - beta-space mixing is an approximation of read-level mixing;
   - low immune fractions (<5%) are detected in only 27–72% of simulations;
   - immune fractions are estimates with chromosome-block bootstrap CIs, not
     exact measurements.

3. **Reference–target comparability.** The reference (Loyfer et al. 2023,
   GSE186458) is WGBS with coverage; the targets are coverage-free bigWig
   percentages. Marker missingness in targets ≤ 0.4%. Endothelial and
   fibroblast references come from non-thyroid tissues (aorta/saphenous vein,
   colon/heart) — stromal fractions are the least interpretable components.

4. **Histone-track input loss for O-information recomputation.** The
   original O-information used divInput per-bp signals (lost). The recomputed
   estimator (5-state beds) matches stored values at r = 0.848 across-region
   but does NOT reproduce stored age slopes (r = 0.29; no significant age
   effect under recomputation). Therefore the Step 9 immune-exclusion
   sensitivity analysis is inconclusive about the original finding
   (see `results/qc/step9_recomputed_slope_divergence.md`).

5. **Methylation-based deconvolution ≠ correction of all seven tracks.**
   Estimated immune fractions reflect DNA methylation only. They cannot
   correct histone-mark or RNA-level composition effects; the O-information
   adjustment (Step 8) is a one-covariate sensitivity analysis, not a
   complete cell-composition correction of the 7-track system.

6. **Cross-sample ChIP-seq quantification.** Histone bigWig absolute values
   are input-normalized but not guaranteed comparable across samples; Step 10b
   therefore reports rank and locus-relative (sample self-calibrated) metrics
   as primary. Two input files were partially unusable:
   `CEMT0044.H3K27me3` bigWig is corrupt (sample excluded, pairwise-complete
   statistics); `CEMT0040.H3K9me3` base file is truncated (the complete
   `.bigWig.1` copy was used instead).

7. **Immune-marker exclusion ≠ purified epithelium.** Removing
   immune-specific methylation bins from bulk data does not remove the immune
   cells' contribution to the remaining bins. Step 9 is a marker-exclusion
   sensitivity analysis, not a reconstruction of epithelium-intrinsic
   O-information.

8. **Sum-to-one constraint.** Deconvolution forces proportions to sum to 1
   over the modeled cell types. Unmodeled cell types (e.g. adipocytes,
   pericytes) load onto the modeled components; fractions are relative to the
   modeled set.

9. **Correlation ≠ causation.** Even where immune fraction correlates with
   age or O-information, no causal claim is made in either direction.
