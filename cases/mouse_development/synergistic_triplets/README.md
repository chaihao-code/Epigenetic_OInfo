# cases/mouse_development/synergistic_triplets

Directory index and supporting files for the parent workflow.

82 samples, seven signals and 35 triplets per sample/class. New triplet calls use 10,000 permutations and BH across all 35 lower-tail tests in each sample/class.

## Inputs

The parent case's sample manifests, existing four-column `.bp` tracks (`chr start end signal`) and unchanged genome-specific annotation masks. Comparison/plot programs also read the existing result CSV/TSV files named in their source.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. Read `cases/mouse_development/synergistic_triplets/f1_by_sample_block.tsv` with pandas (`sep="\t"`; gzip compression is inferred). Existing result tables do not execute a pipeline.

## Outputs

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Files

- `analysis_metadata.json`
- `deleted_runs.json`
- `f1_by_sample_block.tsv`
- `f1_by_stage_block.tsv`
- `input_snapshot.tsv`
- `plot_validation.json`
- `significant_negative_triplets_10000.tsv`
- `special_triplet_f1_inputs.tsv`
- `special_triplet_f1_metadata.json`
- `stage_coverage.tsv`
- `summary_by_block.tsv`
- `summary_by_sample.tsv`
- `summary_by_sample_block.tsv`
- `summary_by_stage.tsv`
- `summary_by_triplet.tsv`
- `summary_by_triplet_class.tsv`
- `synergistic_triplet_mark_trajectories_7signals_by_sample.tsv`
- `synergistic_triplet_mark_trajectories_7signals_metadata.json`
- `synergistic_triplet_mark_trajectories_7signals_summary.tsv`
- `total_synergistic_count_by_stage_region.tsv`
- `total_synergistic_triplet_count_by_stage.tsv`
- `triplet_ABC_and_F1_panels_ab_layout.json`
- `triplet_ABC_selected_region_trajectories_7signals_by_sample.tsv`
- `triplet_ABC_selected_region_trajectories_7signals_metadata.json`
- `triplet_ABC_selected_region_trajectories_7signals_summary.tsv`
- `triplet_ABC_stacked_layout.json`
- `triplet_A_synergistic_count_by_stage.tsv`
- `triplet_B_synergistic_count_by_stage.tsv`
- `triplet_C_synergistic_count_by_stage.tsv`
- `triplet_counts_plot_metadata.json`
- `validation.json`
