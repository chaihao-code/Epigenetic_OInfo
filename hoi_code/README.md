# HOI Analysis — Higher-Order Interaction & O-Information Pipeline

This directory contains the higher-order interaction (HOI) analysis pipeline
of the `sequence_entropy` project. Starting from the binned genome-wide signal
tracks produced by the `preprocessing/` pipelines, it discretizes the signals,
partitions the genome into functional blocks, and computes joint entropies,
O-information (Ω), total correlation (TC), dual total correlation (DTC),
mutual information, and pairwise correlations per block.

## File Overview

| File | Role |
|:-----|:-----|
| `signal_loader.py` | Module: load binned signal files and discretize them |
| `block_processor.py` | Module: genomic block definitions and block combinations |
| `compute_hoi_entropy.py` | Main entry: joint entropy of all feature subsets per block |
| `compute_oinfo_from_entropy.py` | O-info / TC / DTC from per-block entropy files |
| `compute_oinfo_from_combined_entropy.py` | Same, for combined per-cell-type entropy files |
| `convert_entropy_to_mi.py` | Convert entropy subsets to MI tables for plotting |
| `compute_signal_correlation.py` | Pairwise signal correlation (overall or per block) |
| `submit_correlation_jobs.sh` | SLURM batch submission for the correlation jobs |

## Modules

### `signal_loader.py` — signal loading and discretization

- Reads 4-column signal files (`chr start end signal`, the `.bp` output of
  `preprocessing/generate_bin_avg`), keeping autosomes only.
- Discretizes each signal into `n_bins` quantile bins (**bin 0 always
  contains zero-signal positions**; non-zero values are binned separately by
  quantiles).
- Builds a wide-format matrix: rows = genomic bins (`bin_id`), columns =
  features (e.g. ATAC, H3K4me1, methylation, ...).

### `block_processor.py` — genomic block partitioning

- Block files are TSVs with a `bin_id` column and a 0/1 flag column named
  after the block (e.g. promoter, enhancer, exon, intron, CpG island, LINE,
  SINE, LTR — see `data/genome_features/mm10/`).
- Within a block group, blocks are mutually exclusive; across groups, the
  **Cartesian product** of all block combinations is computed (plus a `None`
  class for bins not covered by any block in a group). Each combination is
  intersected to yield the final set of genomic bins for that block.

## Main Pipeline

### Step 1 — `compute_hoi_entropy.py` (main entry)

For every block combination:

1. Extracts the block's bins from the signal matrix (optional downsampling to
   `n_max`, default 100000).
2. Computes the **joint entropy H(S) of all 2^n − 1 non-empty feature
   subsets**. With `method=binning` (default) this is exact frequency
   counting on the discretized data; continuous estimators (`gc`, `knn`,
   `kernel`, `gauss`) are available via the
   [`hoi`](https://github.com/brainets/hoi) package.
3. Writes per-block `*_entropy_all_subsets.csv` plus `*_summary.json` /
   `*_metadata.json`, and combined per-sample outputs.

Parallelized with multiprocessing + shared memory (the signal matrix is
shared across workers); supports `--resume` to skip finished blocks.

```bash
python compute_hoi_entropy.py \
    -i input_files.txt \           # list of signal .bp files (one per line)
    -o results/ \
    --block-groups "group1_a.tsv;group1_b.tsv::BLOCKGROUP::group2_a.tsv;group2_b.tsv" \
    --n-bins 5 --n-max 100000 --method binning
```

### Step 2 — O-information from entropies

`compute_oinfo_from_entropy.py` (per-block files) and
`compute_oinfo_from_combined_entropy.py` (combined per-cell-type files with a
`block` column) turn the entropy subsets into information-theoretic measures
for the full system and every leave-k-out subset:

- **TC** = Σ H(Xᵢ) − H(X₁…Xₙ) (total correlation)
- **DTC** = Σ H(X₋ᵢ) − (n−1)·H(X₁…Xₙ) (dual total correlation)
- **Ω = TC − DTC** (O-information; Ω > 0 → redundancy-dominated,
  Ω < 0 → synergy-dominated)

Output: `*_oinfo_overall.csv` and `*_oinfo_missing_{k}.csv` per block.
Both scripts parallelize over files and skip existing outputs.

### Step 3 — `convert_entropy_to_mi.py`

Converts entropy subsets into MI tables for the plotting scripts
(`plot_by_block_all_anchors.py`, `plot_joint_vs_pairwise_ratio_heatmap_v3.py`):

- `*_mi_self.csv` — per-feature entropy H(X)
- `*_mi_pairwise.csv` — pairwise MI: I(X;Y) = H(X) + H(Y) − H(X,Y)
- `*_mi_joint.csv` — joint MI: I(Xᵢ; X₋ᵢ) = H(Xᵢ) + H(X₋ᵢ) − H(all)

### Pairwise correlation — `compute_signal_correlation.py`

A lightweight companion analysis: pairwise Pearson / Spearman / Kendall
correlation between all signals, either over the whole genome or per block
(`--per-block`), on discretized (default) or continuous (`--continuous`)
signals. Outputs a 3-column TSV (`feature1 feature2 correlation`, plus a
`block` column in per-block mode).

`submit_correlation_jobs.sh` submits these jobs to SLURM for the
`mouse_after_E11.5` and `mouse_preimp` projects, reusing the same block
groups and `n_bins=5` configuration as the entropy runs (default: per-block
mode; `--overall` for whole-genome correlation).

## Dependencies

- Python packages: `numpy`, `pandas`
- [`hoi`](https://github.com/brainets/hoi) (with JAX) — only needed for the
  continuous entropy estimators and the `Oinfo`/`TC`/`DTC` metrics; the
  default `binning` path is self-contained
- SLURM (`sbatch`) for `submit_correlation_jobs.sh`

## Notes

- Several scripts contain hard-coded absolute paths (cluster-specific);
  adjust them before reuse.
- `__pycache__/` is excluded from version control.
