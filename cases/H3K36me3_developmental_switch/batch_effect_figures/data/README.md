# cases/H3K36me3_developmental_switch/batch_effect_figures/data

Directory index and supporting files for the parent workflow.

The eight-signal BH56 comparison reuses preimplantation/late-development observations. Full-panel occurrence uses the retained sign-only rule. Replacement-background redundancy cores have a different null from permutation-tested triplets.

## Inputs

The parent case's sample manifests, existing four-column `.bp` tracks (`chr start end signal`) and unchanged genome-specific annotation masks. Comparison/plot programs also read the existing result CSV/TSV files named in their source.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. Read `cases/H3K36me3_developmental_switch/batch_effect_figures/data/intersection_fraction_each_variant.tsv` with pandas (`sep="\t"`; gzip compression is inferred). Existing result tables do not execute a pipeline.

## Outputs

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Files

- `bulk_replacement_plot_format.json`
- `intersection_fraction_each_variant.tsv`
- `intersection_fraction_summary.tsv`
- `late_common_tissue_means.tsv`
- `mark_region_participation.tsv`
- `stable_intersection_triplet_summary.tsv`
- `stable_intersection_triplets.tsv`
- `triplet_support_across_variants.tsv`
- `triplet_support_distribution.tsv`
- `triplet_support_distribution_validation.json`
- `unaffected_triplet_validation.tsv`
