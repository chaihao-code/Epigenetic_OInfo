# QC note — recomputed vs. stored O-information age slopes (Step 9)

Date: 2026-10-09

## Background
The stored per-region O-information values (`results_oinfo_n_bins_5/`) were
computed by the original pipeline from input-normalized per-bp signals. The
divInput per-bp files used for the histone tracks are no longer available
(storage volume offline, confirmed unrecoverable). An exact reimplementation
(`scripts/hoi_recompute.py`) reproduces the pipeline from the 5-state
discretized beds (2026-06-01 matrices).

## Validation status
- Region definitions, bin universe, methylation-track entropies: exact match.
- Stored vs. recomputed O-information values across 40 regions × 6 samples:
  Pearson r = 0.848 (n = 240). Most of this correlation is between-region
  (region means differ strongly); see `hoi_recompute_validation.tsv`.

## New finding (Step 9)
Age slopes (O-information ~ age, n = 6 donors) from the recomputed estimator
(`arm A_all`, `results/tables/oinfo_region_filtering.tsv`) do NOT reproduce
the stored age slopes (`results/tables/oinfo_age_adjustment.tsv`):

- Pearson r = 0.29 between the two sets of 40 slopes; sign agreement 40%.
- The recomputed estimator shows no significant age effect in any of the 40
  contexts (all p > 0.6), including the CTCF-bound proximal enhancer contexts
  where the stored analysis found FDR-significant declines.

## Consequence for interpretation
The immune-marker-exclusion sensitivity analysis (Step 9) is therefore only
valid *within* the recomputed estimator: under re-derived 5-state
discretization, excluding immune-specific methylation bins does not change
age slopes beyond what matched bin-count downsampling produces (arm B ≈ arm M
for all focus regions; `oinfo_filtering_summary.tsv`).

It CANNOT be used to confirm or refute the original stored age effect, because
that effect is not reproduced by the recomputed estimator in the first place.
The likely cause is the loss of the original histone-track input signals
(divInput per-bp files); the 5-state discretized beds preserve coarser
information and their within-region, across-sample variance appears dominated
by factors that differ from the original computation.

## What remains valid
- Step 8 (immune-fraction adjustment of the stored O-information values) uses
  the original stored values and is unaffected.
- The Step 9 result stands as: "under the only fully reproducible estimator
  available today, the immune-bin exclusion has no effect on age slopes;
  however this estimator does not reproduce the original age decline, so the
  sensitivity analysis is inconclusive regarding the original finding."

## Decision (recorded)
No downstream claim is made from the recomputed age slopes.
