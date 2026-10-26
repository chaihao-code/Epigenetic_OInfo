# cases/aging/neutrophil/scripts

Analysis programs or support modules for the parent workflow. O-information arithmetic is provided only by common/02_oinformation/compute_hoi_oinfo.py.

12 retained samples; original O-information and age-regression reference tables. This cleanup does not invent a new neutrophil triplet study.

## Inputs

The parent case's sample manifests, existing four-column `.bp` tracks (`chr start end signal`) and unchanged genome-specific annotation masks. Comparison/plot programs also read the existing result CSV/TSV files named in their source.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. For the shortest verified workflow, follow `demo/README.md`. For a case command preview, use `python tools/plan_oinfo.py --case <case> --sample <sample>`. Program CLI options are listed below where present.

## Outputs

Each script writes the result/figure paths declared in its source. New O-information runs write all-subset entropy, observed Omega/TC/DTC, a permutation summary and metadata; triplet-null runs also write a BH-corrected significance table. Local/sampled-bin exports are enabled with `--save-local-oinfo --save-sampled-bins`.

## Files

- `bin_histone.sh`
- `discretize.sh`
- `divinput.sh`
- `make_datainfo.sh`
- `oinfo_age_regression.py`
- `plot_oinfo_vs_age.py`
- `plot_paper_figures.py`
- `plot_volcano_oinfo_age.py`
- `run_all.slurm`
- `submit_oinfo_hoi.slurm`
- `weighted_methylation.sh`
