# cases/mouse_development/script

Analysis programs or support modules for the parent workflow. O-information arithmetic is provided only by common/02_oinformation/compute_hoi_oinfo.py.

82 samples, seven signals and 35 triplets per sample/class. New triplet calls use 10,000 permutations and BH across all 35 lower-tail tests in each sample/class.

## Inputs

The parent case's sample manifests, existing four-column `.bp` tracks (`chr start end signal`) and unchanged genome-specific annotation masks. Comparison/plot programs also read the existing result CSV/TSV files named in their source.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. For the shortest verified workflow, follow `demo/README.md`. For a case command preview, use `python tools/plan_oinfo.py --case <case> --sample <sample>`. Program CLI options are listed below where present.

- `plot_synergistic_triplet_mark_trajectories.py`: `python cases/mouse_development/script/plot_synergistic_triplet_mark_trajectories.py --help`; options include `--plot-only`.
- `plot_triplet_A_B_C_counts_trajectories.py`: `python cases/mouse_development/script/plot_triplet_A_B_C_counts_trajectories.py --help`; options include `--distal-ctcf`, `--prepare-only`.
- `plot_triplet_counts_corrected_by_stage.py`: `python cases/mouse_development/script/plot_triplet_counts_corrected_by_stage.py --help`; options include `--stacked-only`.

## Outputs

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Files

- `analyze_seven_signal_ABC_regions.py`
- `plot_special_triplet_f1.py`
- `plot_synergistic_triplet_mark_trajectories.py`
- `plot_triplet_ABC_and_F1_panels.py`
- `plot_triplet_ABC_selected_regions.py`
- `plot_triplet_A_B_C_counts_trajectories.py`
- `plot_triplet_A_C_D_counts_trajectories.py`
- `plot_triplet_counts_corrected_by_stage.py`
- `recalculate_triplet_fdr.py`
- `submit_entropy_computation.sh`
