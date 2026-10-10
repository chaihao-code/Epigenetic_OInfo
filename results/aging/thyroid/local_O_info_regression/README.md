# Thyroid local O-information — immune-adjusted age regression + GSEA

Assesses whether the age-associated change of local O-information in thyroid
genomic regions survives adjustment for the estimated immune-cell fraction
(cell-composition confounding check).

## Pipeline

| Step | Script | What it does |
|------|--------|--------------|
| 01 | `scripts/01_regress_oinfo_immune.py` | Per 1-kb bin: `O-info ~ age + immune_fraction` (OLS, 6 donors), plus the unadjusted `O-info ~ age` model for comparison. Fully vectorized least squares. |
| 02 | `scripts/02_gsea_pathway.py` | TSS-based gene mapping (Gencode v47, ±2 kb), gene ranking = `-log10(p_age) * sign(beta_age)` from the **adjusted** model, preranked GSEA against GO_Biological_Process_2021 and MSigDB_Hallmark_2020. |
| 03 | `scripts/03_gsea_marginal.py` | Signal-change GSEA: bins classified per mark by signal-vs-age trend (up/flat/down, same thresholds as the original analysis), ranking = adjusted metric, custom mark x trend gene sets. |
| 04 | `scripts/04_summarize.py` | Collects significant terms (NOM p < 0.05 AND FDR q < 0.25) across regions. |

Run everything: `bash scripts/run_all.sh` (from any directory).

## Inputs (all read-only, see `scripts/config.py`)

- `../local_oinfo_age_correlation/<region>/oinfo_per_sample.bed` — per-bin
  local O-information for the 6 donors, from `run_analysis_generic_excluded.py`.
- `../local_oinfo_age_correlation/<region>/extracted_signals.tsv` — per-bin
  discrete states (0-4) of the 7 marks x 6 donors.
- `../deconvolution/.../thyroid_cell_proportions.tsv` — cell fractions
  (Panel A / Model A). Total immune fraction = sum of T_CD4, T_CD8, B_cell,
  NK, Monocyte, Granulocyte per donor:
  CEMT0098 4.8%, CEMT0130 1.2%, CEMT0132 1.8%, CEMT0042 8.7%,
  CEMT0044 10.9%, CEMT0040 5.4%.

## Outputs

- `results/regression/<region>.tsv` — one row per bin: original 6 sample
  O-info columns; `beta_age, se_age, t_age, p_age, beta_immune, se_immune,
  t_immune, p_immune, r2_model` (adjusted model); `slope_age_unadj,
  se_age_unadj, p_age_unadj`; `age_slope_attenuation`.
- `results/tables/gsea_{pathway,marginal}_significant.tsv` — summary of
  significant terms across all 37 regions.
- Per-region GSEA reports (`*.report.csv`), ranked gene lists and gseapy plot
  directories are written next to these tables when the pipeline runs but are
  not tracked in git (bulk, regenerable from scripts + 01 outputs).
- `logs/` — one log per step (versions, parameters, per-region summaries).

## Key decisions

- The regression reproduces the original unadjusted model's structure (one
  regression per bin, OLS, no covariates beyond age/immune) so results are
  directly comparable with the published analysis.
- GSEA parameters mirror the original pipeline: gene ranking
  `-log10(p) * sign(beta)`, min set size 5, 1000 (pathway) / 10000 (marginal)
  permutations, seed 42, significance = NOM p < 0.05 AND FDR q < 0.25.
  For marginal GSEA on regions with > 50,000 ranked bins, permutations are
  reduced to 1,000 to keep gseapy prerank within RAM (logged per region).
- With n = 6 donors and 2 predictors, residual degrees of freedom = 3;
  p-values are exact t-test p-values but necessarily coarse — treat bin-level
  significance as descriptive.

## Dependencies

Python 3.8+, numpy, pandas, scipy, gseapy 1.0.3 (downloads gene-set libraries
from Enrichr at run time), bedtools (TSS intersect).
