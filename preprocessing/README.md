# Preprocessing Pipelines — Genome-Wide Signal Binning

This repository contains the preprocessing pipelines of the `sequence_entropy`
project. They convert raw epigenomic signal tracks (ATAC-seq, histone
ChIP-seq, and DNA methylation) into uniformly binned, whole-genome signal
profiles, which serve as the input for downstream sequence-entropy and
correlation analyses.

All pipelines produce the same 4-column BED-like format:

```
chr   start   end   mean_signal
```

## Repository Layout

```
preprocessing/
├── generate_bin_avg/               # Signal binning toolkit (bigWig/bedGraph -> binned mean)
│   ├── compute_bin_average.sh          # Core: bin one bigWig/bedGraph file
│   ├── convert_cpg_to_bedgraph.sh      # per-CpG BED (Bismark/NOMe-seq) -> bedGraph
│   ├── compute_weighted_methylation_average.sh  # Coverage-weighted WGBS binning
│   ├── norm_using_fc_to_input.sh       # ChIP / Input fold-change normalization
│   ├── batch_*.sh / slurm_*.sh         # Batch / SLURM wrappers
│   └── README.md                       # Detailed documentation (Chinese)
├── average_replicates_and_liftover/  # Single-cell replicate averaging + mm9->mm10 liftOver
│   ├── average_one_replicate_group.sh  # Average one replicate group + liftOver
│   ├── submit_replicate_avg_jobs.sh    # Batch submission for all groups
│   ├── merge_stage_schema.slurm        # Merge embryos of the same stage
│   ├── submit_stage_merge_jobs.sh      # Batch submission for stage merging
│   ├── mm9ToMm10.over.chain.gz         # UCSC liftOver chain file
│   └── mm9.chrom.sizes
└── README.md
```

## Pipeline 1: `generate_bin_avg` — Signal Binning

**Core script: `compute_bin_average.sh`**

1. Builds sliding windows with `bedtools makewindows -w BIN_SIZE -s STEP_SIZE`
   over a `chrom.sizes` file.
2. Computes the mean signal per window:
   - **bigWig** input: `bigWigAverageOverBed` (the `mean0` column is used).
   - **bedGraph** input: intervals are expanded to 1-bp resolution, then
     averaged with `bedtools map -o mean`.
3. Writes `${BASENAME}.bin${BIN_SIZE}.step${STEP_SIZE}.bp`
   (`chr  start  end  mean_signal`).

Supporting scripts:

- `convert_cpg_to_bedgraph.sh` — converts per-CpG BED files (Bismark /
  NOMe-seq output) to bedGraph, filtering CpGs below a coverage threshold
  (default >= 8).
- `compute_weighted_methylation_average.sh` — pairs WGBS `coverage` and
  `methylation` bigWigs, filters CpGs by coverage, and bins both
  `methylation x coverage` and `coverage` separately so that a
  coverage-weighted mean methylation can be computed downstream.
- `norm_using_fc_to_input.sh` — normalizes binned ChIP signal tracks by
  dividing by the matching Input track (fold change).
- `batch_compute_bin_average.sh` + `slurm_compute_bin_average.sh` — submit one
  SLURM job per file for a whole directory.

## Pipeline 2: `average_replicates_and_liftover` — Replicate Averaging + liftOver

Designed for single-cell methylation bigWigs in **mm9** coordinates, which
must be averaged across replicates and converted to **mm10** before binning.
It proceeds in two stages:

**Stage 1 — per-group averaging with liftOver**
(`average_one_replicate_group.sh`, driven by `submit_replicate_avg_jobs.sh`):

For each `(group, kmer)` combination (e.g. `2-Cell_embryo3.ACG.TCG`):

1. `bigWigToBedGraph` on every replicate bigWig;
2. expand intervals to 1-bp resolution;
3. `liftOver` mm9 -> mm10 (`mm9ToMm10.over.chain.gz`), keeping only uniquely
   mapped 1-bp intervals;
4. merge all replicates with `bedtools unionbedg` and average over the samples
   that actually have a value at each position.

Output: `${GROUP}.${KMER}.mm10.bedGraph`.

**Stage 2 — stage-level merging**
(`merge_stage_schema.slurm`, driven by `submit_stage_merge_jobs.sh`):

Embryos belonging to the same developmental stage and experiment schema are
merged in the same way (`unionbedg` + mean over covered samples).

Output: `${STAGE}.${SCHEMA}.stage_mean.bedGraph`.

## End-to-End Example: Mouse Preimplantation Data

Applied to `data/mouse_preimplementation/`:

- `original/atac/` — ATAC-seq bigWigs per developmental stage (already mm10);
- `original/histone/` — histone ChIP-seq bigWigs (already mm10);
- `original/methylation/` — single-cell NOMe-seq bigWigs (**mm9**), with k-mer
  context in the file names (e.g. `2-Cell_embryo3.ACG.TCG.bw`).

Workflow:

```
methylation (mm9, single cell)
    -> submit_replicate_avg_jobs.sh   # Stage 1: replicate average + liftOver
    -> original/methylation_rep_averaged/
    -> submit_stage_merge_jobs.sh     # Stage 2: merge embryos per stage
    -> original/methylation_stage_averaged/

atac / histone / methylation_stage_averaged
    -> submit_{atac,histone,methylation}.sh
       (all call batch_compute_bin_average.sh)
    -> 100000_bin_1000_step_binned_results/{atac,histone,methylation}/
```

The final binning uses **BIN_SIZE = 100000, STEP_SIZE = 1000**, i.e. each
100-kb window slides 1 kb at a time (neighboring windows overlap by 99 kb),
producing a smooth, continuous genome-wide signal track:

```
chr1    0       100000    0
chr1    1000    101000    0
chr1    2000    102000    0
```

## Dependencies

- `bedtools`
- `deepTools` (provides `bigWigAverageOverBed`)
- UCSC utilities: `bigWigToBedGraph`, `bedGraphToBigWig`, `liftOver`
- SLURM (`sbatch`) for batch submission

## Notes

- Job logs under `logs/` are excluded from version control.
- Several scripts contain hard-coded absolute paths; adjust them before reuse
  on a different system.
- `generate_bin_avg/README.md` provides a more detailed description of the
  binning toolkit (in Chinese).
