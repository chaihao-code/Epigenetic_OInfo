# cases/mouse_development/result_oinfo_n_bins_5_triplet_perm10000/code

Retained numerical outputs and/or new explicitly named demo outputs for the parent analysis. The metadata in each run identifies its actual permutation count and signal panel.

82 samples, seven signals and 35 triplets per sample/class. New triplet calls use 10,000 permutations and BH across all 35 lower-tail tests in each sample/class.

## Inputs

The parent case's sample manifests, existing four-column `.bp` tracks (`chr start end signal`) and unchanged genome-specific annotation masks. Comparison/plot programs also read the existing result CSV/TSV files named in their source.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. For the shortest verified workflow, follow `demo/README.md`. For a case command preview, use `python tools/plan_oinfo.py --case <case> --sample <sample>`. Program CLI options are listed below where present.

- `block_processor.py`: `python cases/mouse_development/result_oinfo_n_bins_5_triplet_perm10000/code/block_processor.py --help`; options include `--signal-matrix`, `--block-files`, `--block-groups`, `--output-dir`, `--n-max`, `--seed`, `--quiet`.
- `compute_hoi_oinfo.py` forwards imports and command-line execution to the single shared engine.
- `signal_loader.py`: `python cases/mouse_development/result_oinfo_n_bins_5_triplet_perm10000/code/signal_loader.py --help`; options include `--input`, `--output`, `--n-bins`, `--quiet`.

## Outputs

Each script writes the result/figure paths declared in its source. New O-information runs write all-subset entropy, observed Omega/TC/DTC, a permutation summary and metadata; triplet-null runs also write a BH-corrected significance table. Local/sampled-bin exports are enabled with `--save-local-oinfo --save-sampled-bins`.

## Files

- `block_processor.py`
- `compute_hoi_oinfo.py`
- `correction_utils.py`
- `finish_analysis.py`
- `signal_loader.py`
