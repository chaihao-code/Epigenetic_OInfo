# Response to reviewer — immune cell composition confounding (thyroid aging O-information)

We thank the reviewer for raising the important possibility that age-associated
immune infiltration contributes to the bulk-tissue O-information changes we
report. We performed a dedicated cell-composition confounding assessment
(Loyfer et al., Nature 2023 reference atlas; constrained NNLS methylation
deconvolution; RNA-seq marker validation; immune-region exclusion sensitivity
analysis). All analyses were run on the six donor samples of the original
study (ages 30–67 y). Key numbers: `results/tables/`, figures:
`results/figures/`, methods: `reports/methods.md`.

## 1. Is there age-associated immune infiltration?

We observe a **mild, non-significant increase** in estimated immune fraction
with age. Total immune fraction (methylation deconvolution, simulation-validated,
chromosome-block bootstrap CIs) ranges 1.2–10.8% across the six donors and
correlates positively with age (Pearson r = 0.50, p = 0.31; Spearman ρ = 0.64,
p = 0.17; slope ≈ +0.14 percentage-points/year, 95% CI [−0.19, +0.47]).
RNA-seq marker genes corroborate the trend: PTPRC (ρ = 0.81, p = 0.050),
LST1 (ρ = 0.84, p = 0.036), CD3D (ρ = 0.75, p = 0.084), and the estimated
immune fraction itself is supported by expression of immune markers
(PTPRC ρ = 0.83, p = 0.042). At n = 6 these trends are suggestive, not
definitive; we state them as such.

## 2. Is the deconvolution reliable?

Within its stated limits. Because the WGBS methylation files carry no
per-CpG coverage, fragment-level methods (UXM) are not applicable; we used
coverage-free constrained NNLS on methylation beta values, validated on
simulated mixtures: zero bias at 0% immune content, MAE 0.002–0.012, and
27–72% detection of a 5% immune fraction. We therefore treat fractions below
~5% as low-precision and report bootstrap confidence intervals rather than
point estimates. The thyroid-epithelial fraction (77–88%) and fit quality
(R² 0.91–0.96) are stable across two independent marker panels (panel
agreement r = 0.993).

## 3. Does the O-information age effect survive adjustment for immune fraction?

**Yes, essentially unchanged.** In the CTCF-bound proximal enhancer contexts
where O-information declines significantly with age (BH-FDR < 0.05 in four
sub-contexts), adding estimated immune fraction as a covariate attenuates the
age slope by only 0.5–11% (e.g. focus context −0.0396 → −0.0352 per year),
with no significant immune-fraction coefficient anywhere (all p > 0.20).
Age and immune fraction are only moderately correlated (r = 0.50, VIF 1.33),
so the two effects are statistically separable; results are stable to
leave-one-donor-out, donor bootstrap, and marker-panel changes.

## 4. Is the result robust to exclusion of immune-specific methylation regions?

The immune-marker-exclusion sensitivity analysis is **inconclusive** and we
report it as such. The original O-information values were computed from
input-normalized per-bp signals that are no longer available; our exact
reimplementation from the archived 5-state discretized data reproduces the
stored region-level values well (r = 0.848) but does not reproduce the stored
age slopes (r = 0.29). Within the reproducible estimator, excluding
immune-specific methylation bins (1,662 of 2,875,012 autosomal bins; 339
strongly specific) changes age slopes no more than matched bin-count
downsampling (100 random subsets per region). We therefore cannot claim the
exclusion analysis validates the original finding, and we have added this as
an explicit limitation.

## 5. Do the GO enrichments indicate composition confounding?

No evidence for it. First, the immune-related GO terms (T-cell activation,
leukocyte adhesion, MHC class II antigen presentation) were nominally trending
but **not genome-wide significant** in the original GSEA (FDR q = 0.12–0.53);
the only FDR-significant terms are DNA-repair/organelle-organization terms.
Second, bins with significant age-declining local O-information are **not
enriched** for immune-specific methylation bins in any region (permutation
p = 0.66–1.00). Third, re-running the GSEA after excluding immune-specific
bins leaves the immune-term enrichment scores and FDRs essentially unchanged
(MHC-class-II terms trend toward declining O-information with NES ≈ −1.7,
FDR ≈ 0.11–0.16, both before and after exclusion). The immune-marker
expression trend (point 1) and the O-information signal are therefore not
currently linked through shared genomic loci.

## 6. Do the histone-mark/methylation changes track immune fraction?

In the focus context (CTCF-bound proximal enhancers without CpG islands),
sample-self-calibrated (locus-relative) continuous signals change with age in
the directions the reviewer describes: H3K4me1 declines strongly
(slope −0.038/y, p = 0.013; Spearman ρ_age = −0.99), H3K4me3 declines
(ρ = −0.75), and H3K27me3, H3K9me3 and DNA methylation increase
(ρ_age = 0.56–0.67). Notably, the two repressive marks correlate with the
estimated immune fraction (ρ = 0.90, p = 0.037), and the active marks trend
negatively with it — the repressive-mark increases are the component most
consistent with a composition contribution. However, (i) n = 6; (ii) the
rank-based metric (robust to cross-sample normalization) confirms only the
H3K4me1 decline (ρ = −0.81, p = 0.050) and disagrees on H3K27ac; and
(iii) two input files were partially corrupted and handled with pairwise
complete data. We therefore present these as exploratory associations, not
evidence that composition changes drive the histone-mark shifts.

## 7. Revisions to our interpretation

We have revised the manuscript to (i) describe the immune-fraction trends as
suggestive rather than established; (ii) state that controlling estimated
immune fraction does not explain the O-information age decline; (iii) add the
deconvolution and recomputation limitations; and (iv) qualify the
cell-intrinsic interpretation: our data are **consistent with a
thyroid-epithelium-intrinsic component** of the O-information aging signal,
but with n = 6, coverage-free methylation input, and an inconclusive
marker-exclusion test, we do not claim proof of cell-intrinsic epigenetic
aging. We also clarify that methylation-based deconvolution does not correct
the other six epigenetic tracks, so bulk-level vs. cell-intrinsic
contributions to the 7-track O-information cannot be fully separated with
the current data.
