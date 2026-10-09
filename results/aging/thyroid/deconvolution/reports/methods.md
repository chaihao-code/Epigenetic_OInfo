# Methods — Cell-composition confounding assessment of age-associated O-information changes in human thyroid

Date of analysis: 2026-10-09
Project directory: `result/hoi_analysis/projects/aging/thyroid/deconvolution`

## 1. Study design and samples

Six donor thyroid samples from the IHEC/CEEHRC reference epigenome collection
were used in the original O-information analysis (all healthy adults; 5 female,
1 male; ages 30–67 y):

| sample_id | age | sex |
|-----------|-----|-----|
| CEMT0098  | 30  | F   |
| CEMT0130  | 37  | F   |
| CEMT0132  | 37  | F   |
| CEMT0042  | 46  | F   |
| CEMT0044  | 55  | F   |
| CEMT0040  | 67  | M   |

CEMT0086/0088/0090 were excluded by the original analysis and are not used here.
The statistical unit throughout is the donor/sample (n = 6); genomic bins are
never treated as replicates.

Seven epigenetic tracks: DNA methylation (WGBS), H3K4me1, H3K4me3, H3K9me3,
H3K27me3, H3K27ac, H3K36me3. The genome is partitioned into 1-kb bins
(2,875,012 autosomal bins, hg38; chromosome naming `chrN` consistent across all
samples; verified by chromosome-length checks, e.g. chr21 = 46,709,983).

## 2. Input data inventory (Step 1–2)

- DNA methylation: `*.WGB-Seq.methylation_profile.bigWig`, **0–100 percent
  scale, no per-CpG coverage**. Consequence: UXM_deconv (fragment/PAT-level)
  cannot be run; deconvolution uses constrained NNLS on beta values as a
  linear-mixture approximation. This limitation is reported separately.
- Histone marks: `*.signal_unstranded.bigWig` (input-normalized).
  Note: `96920.CEEHRC.CEMT0040.H3K9me3.signal_unstranded.bigWig` is truncated
  (36 MB vs. full `.bigWig.1` ~1 GB); recorded, not used for quantitative work.
- RNA-seq: `*.mRNA-Seq.rpkm_{forward,reverse}.bigWig` for all 6 samples.
- Discrete states: 5-state per-sample discretized beds
  (`data/aging_thyroid_input_normed_binned_discreted/`, row order = bins_1kb
  autosomal order).
- Reference: Loyfer et al., Nature 613:355–364 (2023), GEO GSE186458.
  27 WGBS `.hg38.beta` files (wgbstools format, #meth/#cov uint8 pairs,
  29,401,795 CpGs each), 9 cell types × 3 donors: Thyroid_epithelial
  (GSM5652264–266), T_CD4, T_CD8, B_cell, Monocyte, NK, Granulocyte,
  Endothelial (aorta/saphenous vein), Fibroblast (colon/heart). Plus
  `Markers.U250.hg38.tsv` (UXM marker set) and hg38.2bit.

Full inventory: `results/tables/input_inventory.tsv`,
`results/tables/reference_download_manifest.tsv`.

## 3. Validation of the O-information recomputation (Step 9 prerequisite)

The original per-region O-information values
(`results_oinfo_n_bins_5/`) were compared against an exact reimplementation
(`scripts/hoi_recompute.py`) of the published pipeline
(`compute_hoi_entropy.py` + `compute_oinfo_from_entropy.py`):
5-state discretization, region definitions as the Cartesian product of
{None, promoter, proximal CTCF± enhancer, distal CTCF± enhancer, exon, intron}
× {CpG island, LINE, LTR, SINE, None}, joint entropies from empirical
frequencies, O-info = TC − DTC, n_max = 100,000 bins subsampled with
`pandas.sample(random_state=42)`.

Result: DNA-methylation-track entropies match the original computation
exactly; the methylation O-information values used in Step 8 are the stored
original values. Histone-track entropy inputs (divInput per-bp files) resided
on an unavailable storage volume and could not be recovered; recomputed vs.
stored O-information across all 40 regions × 6 samples correlate at
r = 0.848 (n = 240; `results/qc/hoi_recompute_validation.md`).
**Decision (recorded):** Step 9 three-arm comparisons (all / immune-excluded /
matched downsampling) use the recomputed estimator uniformly, applied to the
stored 5-state matrices (discretization thresholds are NOT re-derived after
filtering, so discretization cannot drive arm differences). Step 8 uses the
stored original O-information values.

## 4. Reference matrix construction (Step 4)

Per-CpG methylation/coverage extracted from the 27 beta files at all 1-kb
reference bins; per-cell-type means over 3 donors. Replicate concordance
(Pearson r of per-bin values between donors of the same cell type) median
0.66 (fibroblasts, cross-tissue) – 0.93. Validation: beta-derived
thyroid-epithelial methylation vs. the GEO-supplied bigWig for GSM5652264
agrees per-CpG at r = 0.9999998; vs. the Loyfer `Markers.U250` target means
at r = 0.67–0.90 (endothelial/fibroblast lower due to cross-tissue mixing —
recorded).

Two marker panels:
- Panel B: UXM U250 markers restricted to the 8 modeled cell types
  (2,149 markers).
- Panel A: independently selected 1-kb-bin markers: |Δ immune vs. other| in
  top range, replicate sd ≤ 0.08, mean coverage ≥ 10; balanced
  1,250 hyper- + 1,250 hypomethylated immune markers.
Immune-specific regions (|d| ≥ 0.4): 666 hyper + 997 hypo bins
(`data/processed/immune_specific_regions.bed`); thyroid-specific
(|d_thyroid| ≥ 0.4): 13,277 bins.

## 5. Deconvolution validation and estimation (Steps 5–6)

Constrained NNLS: min ||W(y − Rp)||² s.t. p ≥ 0, Σp = 1 (augmented-ones
formulation, scipy NNLS iteration; `scripts/deconv_utils.py`).
Two granularity models: A (thyroid/lymphoid/myeloid/stromal), B (thyroid
epithelium, T_CD4, T_CD8, B, monocyte, NK, granulocyte, endothelial,
fibroblast).

Simulation: beta-space mixtures of reference profiles at immune fractions
0–50% (σ = 0.01/0.05, 50 replicates). Zero-immune false-positive bias
≈ 0.0003–0.0019; MAE 0.002–0.012; 5% immune detected in 27–72% of
replicates. → low-level immune fractions have limited precision; results are
reported with chromosome-block bootstrap CIs (200 resamples).

Bulk thyroid estimation: marker missingness ≤ 0.4%, fit R² 0.91–0.96.

## 6. O-information age models (Steps 7–8)

Per genomic context:
- Model 0: O-info ~ age
- Model 1: O-info ~ age + total immune fraction (Panel A/Model A primary;
  Panel B/Model B sensitivity; panel stability r = 0.993).
Attenuation = 1 − β_age,adj/β_age,unadj; donor bootstrap CIs (1,000);
leave-one-donor-out; VIF; BH-FDR across 40 contexts.

## 7. Immune-marker-exclusion sensitivity analysis (Step 9)

Recomputed O-information on (A) all bins, (B) bins excluding immune-specific
regions, (C) bins excluding strongly immune-specific (|d| ≥ 0.5) regions, and
(M) 100 matched random subsets of the same size per region (region composition
matched by construction). This is an immune-marker-exclusion sensitivity
analysis; it does NOT reconstruct epithelium-pure O-information.

## 8. Marker / GO / histone analyses (Step 10)

- RNA-seq markers: gene-body mean RPKM (forward + reverse) from IHEC
  mRNA-Seq bigWigs, gencode v47 gene coordinates.
- Immune-region overlap: losing-redundancy bins (FDR < 0.10, slope < 0) vs.
  immune-specific bins; region background + 1,000 region-matched random
  permutations.
- GO: exact reproduction of the original GSEA prerank protocol
  (ranking = −log10(p)·sign(r), min-p per gene, GO_Biological_Process_2021,
  1,000 permutations, seed 42) on the original and immune-bin-excluded
  TSS-annotated locus tables.
- Histone associations: per-1kb bin means extracted once per file with
  `pyBigWig.values(..., numpy=True)` (cached in `data/processed/bin_means_cache/`);
  region medians computed per sample. Three metrics — raw, locus-relative
  (log2 region median / sample genome-wide median), and cross-sample rank —
  Spearman vs. age and vs. immune fraction with pairwise-complete samples;
  OLS slopes reported as effect sizes (n = 6; cross-sample ChIP-seq absolute
  quantification treated as unreliable, rank/locus-relative preferred).
  Data integrity notes: `97007.CEEHRC.CEMT0044.H3K27me3.signal_unstranded.bigWig`
  is corrupt (fails to open); CEMT0040 H3K9me3 uses the complete
  `.bigWig.1` copy (the base file is a truncated 36 MB stub). Both handled
  with pairwise-complete observations.

## 9. Software

Python 3.8.5 (Ubuntu, kernel 5.15.0-179): numpy 1.24.4, pandas 1.5.3,
scipy 1.10.1, statsmodels 0.13.5, scikit-learn 1.3.2, matplotlib 3.7.1,
seaborn 0.12.2, pyBigWig 0.3.18, py2bit 0.3.0, gseapy 1.0.3, PyYAML 6.0.
Fixed seeds: 42 (entropy subsampling, overlap permutation, GSEA),
20261009 (matched downsampling), 2026 (deconvolution bootstrap).

## 10. Reproducibility

`config.yaml` holds all paths and parameters; scripts `01–10` run
independently in order; see `README.md`.
