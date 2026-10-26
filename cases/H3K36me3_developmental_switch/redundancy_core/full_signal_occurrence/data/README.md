# cases/H3K36me3_developmental_switch/redundancy_core/full_signal_occurrence/data

Directory index and supporting files for the parent workflow.

The eight-signal BH56 comparison reuses preimplantation/late-development observations. Full-panel occurrence uses the retained sign-only rule. Replacement-background redundancy cores have a different null from permutation-tested triplets.

## Inputs

The parent case's sample manifests, existing four-column `.bp` tracks (`chr start end signal`) and unchanged genome-specific annotation masks. Comparison/plot programs also read the existing result CSV/TSV files named in their source.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. Read `cases/H3K36me3_developmental_switch/redundancy_core/full_signal_occurrence/data/counts_by_sample.tsv` with pandas (`sep="\t"`; gzip compression is inferred). Existing result tables do not execute a pipeline.

## Outputs

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Files

- `counts_by_sample.tsv`
- `occurrence_summary.tsv`
- `plot_validation.json`
- `redundancy_core_combination_occurrence_post_E11p5_all.tsv`
- `redundancy_core_combination_occurrence_post_E11p5_top25.tsv`
- `redundancy_core_combination_occurrence_preimplantation_all.tsv`
- `redundancy_core_combination_occurrence_preimplantation_top25.tsv`
- `redundancy_cores.tsv`
- `synergistic_triplet_combination_occurrence_post_E11p5_all.tsv`
- `synergistic_triplet_combination_occurrence_post_E11p5_top25.tsv`
- `synergistic_triplet_combination_occurrence_preimplantation_all.tsv`
- `synergistic_triplet_combination_occurrence_preimplantation_top25.tsv`
- `triplet_combinations_by_sample.tsv`
