# demo/results/from_bigwig/oinfo

## Purpose

This directory contains inputs, intermediate files or outputs of the complete human thyroid demo or its parent workflow. Read the repository and `demo/README.md` for the exact stage definitions. Original genomic windows and class memberships remain unchanged.

## Inputs and usage

Run from the release root. `python demo/run_demo.py` uses the supplied complete prepared tracks; `python demo/run_demo.py --from-bigwig` bins the complete bigWigs on the unchanged original grid, applies histone/(Input+0.01), and calls the shared O-information engine. Use `python demo/compare_runs.py` after both routes complete. Local input links are not portable and are listed in `config/local_input_links.tsv`.

## Outputs

The demo writes original-format raw/normalized `.bp` tracks, subset/global O-information, 10,000-permutation summaries, exact sampled loci, state codes, local O-information and local/global checks. Reference cohort results elsewhere in the release retain their historical settings. The original archived thyroid values have a documented upstream missing-input limitation, reported separately from current-input pipeline validation.

## Files

- `CEEHRC.CEMT0098_discrete_states.tsv.gz`
- `CEEHRC.CEMT0098_entropy_all_subsets.csv`
- `CEEHRC.CEMT0098_local_global_check.tsv`
- `CEEHRC.CEMT0098_local_oinfo.tsv.gz`
- `CEEHRC.CEMT0098_metadata.json`
- `CEEHRC.CEMT0098_oinfo.csv`
- `CEEHRC.CEMT0098_permutation_oinfo.csv`
- `CEEHRC.CEMT0098_permutation_summary.csv`
- `CEEHRC.CEMT0098_sampled_bins.tsv.gz`
