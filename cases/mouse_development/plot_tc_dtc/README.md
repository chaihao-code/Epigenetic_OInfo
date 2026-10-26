# cases/mouse_development/plot_tc_dtc

Analysis programs or support modules for the parent workflow. O-information arithmetic is provided only by common/02_oinformation/compute_hoi_oinfo.py.

82 samples, seven signals and 35 triplets per sample/class. New triplet calls use 10,000 permutations and BH across all 35 lower-tail tests in each sample/class.

## Inputs

The parent case's sample manifests, existing four-column `.bp` tracks (`chr start end signal`) and unchanged genome-specific annotation masks. Comparison/plot programs also read the existing result CSV/TSV files named in their source.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. For the shortest verified workflow, follow `demo/README.md`. For a case command preview, use `python tools/plan_oinfo.py --case <case> --sample <sample>`. Program CLI options are listed below where present.

- `plot_combined_hoi_trajectories_2x4.py`: `python cases/mouse_development/plot_tc_dtc/plot_combined_hoi_trajectories_2x4.py --help`; options include `--results-dir`, `--tree-path`, `--output-dir`.

## Outputs

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Subdirectories

- [inputs](inputs/README.md)

## Files

- `README_zh.md`
- `combined_tc_dtc_trajectories_6x2_n5_validation.json`
- `oinfo_overall_values_for_plot.tsv`
- `plot_combined_hoi_trajectories_2x4.py`
- `plot_tc_dtc_stacked_18cm.py`
- `plot_total_oinfo_trajectories_18cm.py`
- `tc_dtc_missing_inputs.csv`
- `tc_dtc_overall_values_for_plot.tsv`
