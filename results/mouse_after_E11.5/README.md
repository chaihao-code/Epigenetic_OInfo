# Mouse development after E11.5: O-information rerun

66 stage-by-tissue samples (E11.5–E16.5 and P0), 10 signals each.
Input signal paths are in `data/input_files/*.txt`; source tracks are external and retained.

## Current analysis

Engine: `../../script/hoi_code/compute_hoi_oinfo.py` relative to this project directory.
Actual absolute path is recorded in `oinfo_rerun_manifest.json`.
Parameters: 5 discretization bins, at most 100000 loci per block, sampling seed 42,
40 mm10 genomic block combinations, all 10 signals retained (drop-signals disabled).
10000 independent within-block permutations, permutation seed 1234 (engine default).
Null distributions are calculated for the full 10-signal system (system mode),
while observed entropy/O-information is calculated for all subsets supported by the engine.
The engine retains a legacy `9sig` system label even for 10 features; consult feature_names
and remaining_indices to identify the actual signals. No reduced subsystem is requested.

Submit: `sbatch script/submit_oinfo.slurm` (66 array tasks, 16 CPUs/task, 24 h, no concurrency cap).
Sample order: `script/oinfo_samples.tsv`.
Results and Slurm logs: `results_oinfo_n_bins_5/`.
Each successful task writes entropy, observed O-information, permutation O-information,
permutation summary, and run metadata files. FDR correction is not performed by this job.

## Cleanup

At the user's request, previous entropy, replacement/append, pairwise correlation results
and replacement_summary.tsv were deleted before resubmission. Input lists and historical
submission scripts remain. Parameters, source hash, cleanup status and submission ID are
recorded in `oinfo_rerun_manifest.json`.

## GitHub source export

See [UPLOAD_NOTES.md](UPLOAD_NOTES.md) for the uploaded scope, excluded legacy/TF analyses,
shared calculation snapshot, original HPC paths and external data requirements.
