# Mouse development: seven-signal integration

82 samples cover 17 stages from zygote to P0: 7 preimplantation, 9 gastrulation and 66 late tissue-stage samples.
The shared panel is ATAC, DNA methylation, H3K4me1, H3K4me3, H3K27ac, H3K9me3 and H3K27me3.

- `result_oinfo_n_bins_5_triplet_perm10000/`: input manifest, per-sample jobs, submission code and the validated calculation snapshot; five bins, 10,000 permutations, 40 blocks × 35 triplets per sample.
- `script/`: per-sample × block BH35 correction, ABC/F1 plots and developmental trajectories.
- `synergistic_triplets/`: retained BH35 test snapshot and small supporting plotting/summary tables. The original run has 114,800 tests and 20,505 significant negative triplets.
- `plot_tc_dtc/`: total seven-signal O-information, TC and DTC trajectory plotting source and small validated plotting inputs.
- `redundancy_core/`: retained core CDF figure resizing source; the original core computation is outside this small figure directory.
- `shared/hoi_code/`: required shared calculation modules archived for all four mouse projects.

See [UPLOAD_NOTES.md](UPLOAD_NOTES.md) for source-export scope, original HPC paths and external input requirements.
See subdirectory guides for statistical definitions and retained commands.
