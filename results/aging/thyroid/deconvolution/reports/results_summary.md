# Results summary — cell-composition confounding assessment (thyroid aging)

All numbers below are from actual runs in this project directory
(`result/hoi_analysis/projects/aging/thyroid/deconvolution`).
Primary tables are in `results/tables/`, figures in `results/figures/`.
n = 6 donors in all age models.

## Q1. Do thyroid samples differ in immune composition?

Yes, modestly. Total immune fraction (Panel A/Model A, methylation
deconvolution, constrained NNLS, fit R² 0.91–0.96, marker missingness ≤0.4%):

| age | sample | total immune | dominant immune component |
|-----|--------|-------------|---------------------------|
| 30  | CEMT0098 | 4.8% | NK 2.0%, Monocyte 1.5% |
| 37  | CEMT0130 | 1.2% | Granulocyte 1.1% |
| 37  | CEMT0132 | 1.8% | Fibroblast* 4.5% (stromal), Granulocyte 1.0% |
| 46  | CEMT0042 | 8.7% | B_cell 2.8%, T_CD8 2.7% |
| 55  | CEMT0044 | 10.8% | T_CD4 3.5%, B_cell 2.3% |
| 67  | CEMT0040 | 5.4% | Monocyte 2.1%, T_CD4 2.0% |

Thyroid epithelium: 77–88%. Stromal (endothelial+fibroblast): 10–16%.
Simulation-validated precision: bias ≈0 at 0% immune; 5% immune detected in
27–72% of simulations — treat individual small fractions with caution.

RNA-seq marker check (gene-body RPKM, strand-aware): immune markers are
expressed at low level (PTPRC mean 1.3 RPKM, TYROBP 9.3, LST1 2.1, CD3D 0.9);
thyroid markers strongly expressed (TG 423, TPO 89, PAX8 86).

## Q2. Does estimated immune fraction change with age?

Weak positive trend, not significant. Pearson r = 0.45–0.50 (p = 0.31–0.38
across marker-panel/model combinations); Spearman ρ = 0.64 (p = 0.17).
OLS slope ≈ +0.10–0.14 percentage-points/year, 95% CI spans zero
(e.g. Panel A/Model A: +0.139%/y, CI [−0.19, +0.47], p = 0.31).
RNA-seq immune markers corroborate a mild trend: PTPRC ρ_age = 0.81 (p = 0.050),
LST1 ρ = 0.84 (p = 0.036), CD3D ρ = 0.75 (p = 0.084); immune-marker mean
correlates with estimated immune fraction (PTPRC ρ = 0.83, p = 0.042).
**Conclusion: suggestive evidence of a mild age-associated increase in immune
content, not statistically definitive at n = 6.**

## Q3. Is O-information correlated with estimated immune fraction?

No robust association. In Model 1 (O-info ~ age + immune fraction) across all
40 genomic contexts, beta_immune is non-significant everywhere (smallest
p = 0.20). For the focus context (CTCF-bound proximal enhancers without CpG
islands), beta_immune p = 0.50. Age–immune correlation r = 0.50 (VIF 1.33),
so the two covariates are separable.

## Q4. Does the age effect on O-information survive immune-fraction adjustment?

Yes. Four contexts show FDR-significant O-information decline with age
(BH-FDR < 0.05), all CTCF-bound proximal enhancer subtypes (CpG island, LINE,
LTR, SINE). After adding immune fraction:

| context | β_age unadj (p) | β_age adj (p) | attenuation |
|---------|----------------|---------------|-------------|
| CTCF+ prox. enh. × CpG island | −0.0076 (6.7e-5) | −0.0073 (5.7e-4) | 4% |
| CTCF+ prox. enh. × LINE | −0.0380 (0.002) | −0.0378 (0.013) | 0.5% |
| CTCF+ prox. enh. × LTR | −0.0458 (0.004) | −0.0466 (0.021) | −2% (sign flip of Δ) |
| CTCF+ prox. enh. × SINE | −0.0307 (0.002) | −0.0290 (0.012) | 5% |
| CTCF+ prox. enh. × no island (focus) | −0.0396 (0.013) | −0.0352 (0.054) | 11% |

Attenuation is minimal (0.5–11%); bootstrap CIs and leave-one-out analyses
confirm stability; marker-panel sensitivity r = 0.993 between adjusted slope
sets. The age-associated O-information decline is essentially unchanged by
statistical control of estimated immune fraction.

## Q5. Is the O-information change robust to immune-region exclusion?

Inconclusive — with an important caveat. The fully reproducible recomputed
estimator (5-state discretization) does not reproduce the stored age slopes
(r = 0.29, sign agreement 40%; see `results/qc/step9_recomputed_slope_divergence.md`),
so the exclusion analysis cannot confirm/refute the original finding. Within
the recomputed estimator, excluding immune-specific methylation bins (1,662
bins; strongly specific: 339) changes age slopes no more than matched
bin-count downsampling (100 random subsets per region; arm B ≈ arm M
throughout the CTCF+ proximal enhancer family).

## Q6. Do GO enrichments indicate immune composition confounding?

- The immune-related GO terms cited (T-cell activation, leukocyte adhesion,
  MHC class II antigen presentation) appear in the original GSEA with
  FDR q = 0.12–0.53 — nominally trending but not genome-wide significant;
  the only FDR-significant terms (q < 0.05) are DNA-repair/organelle
  organization terms, none immune.
- Losing-redundancy bins (per-bin p < 0.05, slope < 0) are NOT enriched for
  immune-specific methylation bins in any region (empirical permutation
  p = 0.66–1.00; `immune_region_enrichment.tsv`).
- Re-running GSEA after excluding immune-specific bins changes immune-term
  NES/FDR negligibly (`go_before_after_terms.tsv`).
- MHC-class-II/antigen-presentation terms trend toward the
  losing-redundancy direction (NES ≈ −1.7, FDR ≈ 0.11–0.16) in the LINE
  subcontext, before AND after exclusion equally.

**Conclusion: no evidence that the GO trends are driven by immune-specific
methylation regions.**

## Histone-mark and methylation changes vs. immune fraction

Focus region (CTCF-bound proximal enhancers without CpG islands),
locus-relative metric (log2 region median / sample genome-wide median;
sample-self-calibrated; `histone_immune_association.tsv`, fig5):

| mark | slope/y (p) | ρ_age (p) | ρ_immune (p) |
|------|-------------|-----------|--------------|
| H3K4me1 | −0.038 (0.013) | −0.99 (3e-4) | −0.66 (0.16) |
| H3K4me3 | −0.049 (0.14) | −0.75 (0.08) | −0.77 (0.07) |
| H3K27ac | −0.018 (0.28) | −0.35 (0.50) | −0.49 (0.33) |
| H3K27me3 (n=5) | +0.021 (0.32) | +0.67 (0.22) | **+0.90 (0.037)** |
| H3K9me3 (n=5) | +0.016 (0.29) | +0.56 (0.32) | **+0.90 (0.037)** |
| DNA methylation | +0.0009 (0.17) | +0.64 (0.17) | +0.77 (0.07) |

Directions match the reviewer's description (H3K4me1/4me3/K27ac down;
H3K27me3/K9me3/DNA methylation up). The two repressive marks track the
estimated immune fraction (ρ = 0.90, p = 0.037), and the active marks trend
negatively with it — the repressive-mark increases are therefore the
component most plausibly linked to composition change, though n = 6 and
two corrupted input files (CEMT0044 H3K27me3, truncated CEMT0040 H3K9me3
superseded by the complete `.bigWig.1`) limit precision.

Cross-metric caveat: the rank metric (monotonic-transform-robust) confirms
the H3K4me1 age decline (ρ = −0.81, p = 0.050) but disagrees for H3K27ac
(rank ρ_age = +0.93, p = 0.008, vs locus-relative −0.35): the oldest sample
has the highest genome-wide H3K27ac level, so global-vs-locus normalization
changes the sign. We therefore do not draw conclusions about H3K27ac.

## Overall conclusion

At n = 6, bulk thyroid shows a mild, non-significant age-associated increase
in estimated immune content that is supported (also non-significantly) by
RNA-seq immune markers. The age-associated decline of O-information in
CTCF-bound proximal enhancers is essentially unattenuated by statistical
control of immune fraction, and immune-specific methylation regions do not
overlap the declining loci. However, the deconvolution is a methylation-only,
coverage-free approximation with limited power below 5% immune content, and
the marker-exclusion sensitivity analysis is inconclusive because the
recomputed O-information estimator does not reproduce the original age
slopes. The data are therefore **consistent with, but do not prove, a
thyroid-epithelium-intrinsic component of the O-information aging signal**;
claims of cell-intrinsic epigenetic aging should be stated as suggestive and
qualified by these limitations.
