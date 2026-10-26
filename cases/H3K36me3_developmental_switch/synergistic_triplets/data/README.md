# cases/H3K36me3_developmental_switch/synergistic_triplets/data

Directory index and supporting files for the parent workflow.

The eight-signal BH56 comparison reuses preimplantation/late-development observations. Full-panel occurrence uses the retained sign-only rule. Replacement-background redundancy cores have a different null from permutation-tested triplets.

## Inputs

The parent case's sample manifests, existing four-column `.bp` tracks (`chr start end signal`) and unchanged genome-specific annotation masks. Comparison/plot programs also read the existing result CSV/TSV files named in their source.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. Read `cases/H3K36me3_developmental_switch/synergistic_triplets/data/mark_occurrence_per_tissue_by_stage.tsv` with pandas (`sep="\t"`; gzip compression is inferred). Existing result tables do not execute a pipeline.

## Outputs

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Files

- `mark_occurrence_per_tissue_by_stage.tsv`
- `mark_occurrence_per_tissue_by_stage_post_E11p5.tsv`
- `mark_occurrence_per_tissue_by_stage_preimplantation.tsv`
- `mark_occurrence_validation.json`
- `original_development_heatmap_validation.json`
- `original_development_mark_participation.tsv`
- `original_development_mark_participation_by_sample.tsv`
- `trajectory_cleanup_manifest.json`
- `triplet_ABC_selected_region_trajectories_8signals_post_E11p5_n5_by_sample.tsv`
- `triplet_ABC_selected_region_trajectories_8signals_post_E11p5_n5_metadata.json`
- `triplet_ABC_selected_region_trajectories_8signals_post_E11p5_n5_summary.tsv`
- `triplet_ABC_stacked_post_E11p5_n5_by_sample.tsv`
- `triplet_ABC_stacked_post_E11p5_n5_summary.tsv`
- `triplet_ABC_stacked_post_E11p5_n5_validation.json`
- `triplet_ABC_stacked_preimplantation_by_sample.tsv`
- `triplet_ABC_stacked_preimplantation_summary.tsv`
- `triplet_ABC_stacked_preimplantation_validation.json`
