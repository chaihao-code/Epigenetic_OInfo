# Synergistic triplets: 5 bins, 10000 permutations, BH84

7 stages: zygote, 2cell, 4cell, 8cell, morula, ICM, TE.
9 signals give 84 triplets per stage x genomic block; 40 blocks per stage.
All 84 lower-tail P-values are BH-adjusted together, before filtering by sign.
Synergistic: observed O-information < 0 and BH q <= 0.05.
Empirical lower-tail P = (1 + number of null values <= observed) / 10001.

The previous results_n_bins_5 run provided a null only for the complete nine-signal
system, which cannot be used to test individual triplets. A supplemental run in
permutation_results/ computes pair/triplet nulls with --permutation-level all
--permutation-max-order 3, while preserving the completed full-system results.
Observed entropy/O-information is recomputed using identical inputs, 5 bins,
40 mm10 blocks, n_max 100000, sampling seed 42 and permutation seed 1234.
The summarizer verifies every observed triplet against the completed run (1e-9 tolerance).

## Pipeline and outputs

scripts/submit_triplet_null.slurm: seven independent tasks, 16 CPUs each, 24 h, no concurrency cap.
scripts/summarize_fdr.slurm: runs after successful completion of every null task.
Job IDs and submission status: submission.json. Logs: logs/.

After the pipeline completes:
- all_triplets_fdr.tsv.gz: all 23520 tests, including raw P, BH q and significance flags.
- synergistic_triplets_fdr05.tsv: significant negative triplets after BH84.
- summary_by_stage.tsv: tested, negative, nominal significant, FDR significant and removed counts.
- summary_by_stage_region.tsv: the same counts for each of 280 stage-region groups.
- fdr_validation.json: parameters, source checksums and observed-value agreement.

These outputs are created only when all seven triplet null runs have completed;
submission alone does not imply that FDR results are available.
