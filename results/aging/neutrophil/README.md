# Neutrophil Aging — O-information Analysis

O-information analysis of **mature neutrophils (hg38)** from BLUEPRINT/IHEC
across 12 donors (`donor_metadata.tsv`: 7 healthy adult peripheral-blood donors,
2 cord-blood children, 3 G-CSF/Dex-treated donors), including a locus-permutation
significance test and regression of O-information on donor age.

## Data (7 features per donor)

- Histone ChIP-seq (6 marks): H3K4me1, H3K4me3, H3K9me3, H3K27ac, H3K27me3, H3K36me3
  (BLUEPRINT `pval_signal` bigWig, normalized by `signal/(input+0.01)`)
- WGBS methylation (coverage-weighted)
- 1 kb genomic bins, discretized to **5 bins** (global quantiles; zero signal merged
  into the lowest bin), autosomes only

## Pipeline

### 1. Preprocessing (`scripts/`)

`bin_histone.sh` (bigWig → 1 kb bins) → `weighted_methylation.sh` (coverage-weighted
WGBS methylation) → `divinput.sh` (input normalization) → `make_datainfo.sh`
(per-donor signal list under `data_info/`)

### 2. Entropy + O-information + permutation test (recomputed 2026-10-11)

- Program: `script/hoi_code/compute_hoi_oinfo.py` (SLURM array submission:
  `scripts/submit_oinfo_hoi.slurm`, 12 tasks)
- Genomic blocks: updated *purified* annotations (promoter / proximal|distal ×
  ctcf_bound|not enhancer / exon / intron — pairwise overlap verified to be zero)
  × repeat classes (CpG_island / LINE / SINE / LTR), AND-combined into **40 blocks**
- Per donor × block: joint entropies of all 127 non-empty feature subsets
  (exact frequency counting, Shannon log2)
- O-information: Ω = TC − DTC (Rosas et al. 2019; Ω > 0 redundancy, Ω < 0 synergy),
  for the full system and leave-k-out subsets (n_missing 1–5)
- **Significance**: 1000 locus permutations (each track permuted independently
  within each block's sampled loci, preserving block-level marginals); per-block
  `p_dev` (two-sided deviation from the permutation null), `p_up` / `p_down` (direction)
- Outputs: `results_oinfo_n_bins_5/<IHECRE*>/{*_entropy_all_subsets.csv, *_oinfo.csv,
  *_permutation_oinfo.csv, *_permutation_summary.csv, *_metadata.json}`
- Note: the permutation null is an **independence null**. Observed O-information is
  significant in every donor × block (p_dev ≈ 0.002), i.e. genomic signals are
  generically redundant; these p-values do not vary with age and must not be used
  as age-association tests

### 3. Age regression (`scripts/oinfo_age_regression.py`)

- Per block: full-system O-info / TC / DTC regressed on donor age (midpoint of the
  reported range); G-CSF/Dex-treated donors excluded; OLS (slope / R² / p) + Spearman
- Outputs:
  - `oinfo_age_regression.tsv`: all untreated donors (n=9, incl. children)
  - `oinfo_age_regression_adults.tsv`: adults only (n=7)
  - `oinfo_permutation_block_summary.tsv`: per-block permutation significance
    aggregated across donors
- BH-FDR q-values added (`hoi_q`, `hoi_spearman_q`)

### 4. Figures (`plots/`)

- `oinfo_vs_age_{all_blocks,cpg_blocks}[_adults].png`: O-information vs age scatter
  grid per block
- `volcano_oinfo_slope_vs_p_adults.png`: slope vs -log10(p) volcano plot
- `heatmap_oinfo_age.png`: block × donor heatmap
- `scatter_selected_blocks.png`: fitted scatters for the 6 blocks with the
  largest |slope|

## Main results

- **Direction (adults, n=7):** O-information decreases with age in *all* 40 blocks
  (slopes all negative); 11 blocks nominal p < 0.05, strongest in proximal
  CTCF-bound enhancer blocks (LTR/LINE/SINE) and `None_AND_None`
  (≈ 0.2–0.6 bits/decade)
- **After BH-FDR across the 40 blocks: no block survives q < 0.05** (smallest
  q = 0.113, in the `None_AND_*` blocks, R² ≈ 0.47–0.64). With n = 7 the analysis
  is under-powered; the consistent negative direction suggests a real effect that
  more donors could confirm
- Including children (n=9) flips slopes to positive — child donors dominate the
  fit; the adult-only version is the reliable one

## Legacy results

The first version of this analysis (2026-09; older overlapping annotations, no
permutation test, O-information derived in a second pass from entropy tables)
was deleted. This directory is the 2026-10-11 recomputation.
