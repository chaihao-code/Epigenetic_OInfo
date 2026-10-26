# Full-signal synergistic triplet and redundancy core occurrence

Four requested plots compare preimplantation (7 samples, 9 signals) and
post-E11.5 development (all 66 available stage/tissue samples, 10 signals).
Each sample contributes the same 40 annotated genomic blocks. ICM and TE are
separate samples. Gastrulation is outside these two periods.

The common eight signals are ATAC, H3K4me1, H3K4me3, H3K27ac, H3K9me3,
H3K27me3, H3K36me3 and DNA methylation. Preimplantation adds H2A.Z;
post-E11.5 adds H3K9ac and H3K4me2. Source metadata and the complete sets of
subsets are validated for every sample. The legacy `system=9sig` column in the
late source CSVs is not used to infer their signal count; they contain 10.

Triplets follow the supplied occurrence plotting script: all order-three
combinations with observed O-information < 0 (`interpretation == Synergy`).
There is no P-value or FDR filter in these plots. They therefore differ from
the existing significant-triplet trajectory figures. Each identical combination
in a different sample or block counts as another occurrence. A, B and C are
pooled, rather than selecting only one triplet class. There are 84 possible
triplets in the early period and 120 in the late period.

Cores are recomputed from the complete nine- or ten-signal subset tables using
`script/hoi_code/compute_redundancy_core.py`; the existing eight-signal cores
are not reused. For each order from 3 to N-1, test its maximum-O-information
subset against all same-order single-signal replacements. Require positive
O-information, modified z > 2.5, and either a one-component Gaussian mixture
preferred by BIC or a two-component preference with incoming-signal ANOVA
eta > 0.5. Promote an accepted N-1 core to N when full-set O-information is
higher, then prune proper subsets of accepted larger cores. GMM seed is 0,
with five starts. This criterion does not assign a permutation P-value or FDR.

All combinations are ranked by pooled occurrence count, with alphabetical
tie-breaking. Show the top 25 (or all if fewer). Triplet bars show raw counts.
Core bars show occurrence percentages, with denominator equal to all core
occurrences in that period, including combinations outside the top 25. Core
plots also show cumulative percentages with the same denominator. Each column
of the lower dot matrix marks every member of that exact combination; internal
DNA methylation aliases are consistently mapped to its display label.
The first three matrix rows, from top to bottom, are DNA methylation, H3K9me3
and H3K27me3. All remaining signals appear below these in a consistent order.

Individual PDFs/PNGs are 9 x 7.5 cm; an additional 2 x 2 layout is 18 x 15 cm.
All text is Arial 8 pt. Redundancy cores occupy the top row of the combined
figure and synergistic triplets the bottom row. Preimplantation is in the left
column and post-E11.5 in the right column (panels a, b, c, d in reading order).
Different periods have different sample counts and signal universes, so their
raw counts should not be interpreted as directly standardized rates.

Reproduction (use the project `hoi_env` Python):

```bash
python ../scripts/plot_full_signal_combination_occurrence.py --prepare
sbatch ../scripts/compute_full_signal_occurrence_cores.slurm
# After every array task finishes:
python ../scripts/plot_full_signal_combination_occurrence.py --plot
```

The preparation manifest records source and engine checksums. Per-sample core
CSVs, diagnostics and metadata are retained in `results/`. Data tables include
all and top-25 combination frequencies, per-sample counts (including zeros),
pooled core records, and a plot validation report. All new outputs are kept
here; the uploaded three scripts and earlier analyses are unchanged.
