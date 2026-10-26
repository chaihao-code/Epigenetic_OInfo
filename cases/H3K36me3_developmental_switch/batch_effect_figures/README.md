# cases/H3K36me3_developmental_switch/batch_effect_figures

Directory index and supporting files for the parent workflow.

The eight-signal BH56 comparison reuses preimplantation/late-development observations. Full-panel occurrence uses the retained sign-only rule. Replacement-background redundancy cores have a different null from permutation-tested triplets.

## Inputs

The parent case's sample manifests, existing four-column `.bp` tracks (`chr start end signal`) and unchanged genome-specific annotation masks. Comparison/plot programs also read the existing result CSV/TSV files named in their source.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. Read `cases/H3K36me3_developmental_switch/batch_effect_figures/same_stage_mark_summary.tsv` with pandas (`sep="\t"`; gzip compression is inferred). Existing result tables do not execute a pipeline.

## Outputs

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Subdirectories

- [data](data/README.md)
- [fdr_lt_0p1](fdr_lt_0p1/README.md)
- [inputs](inputs/README.md)
- [scripts](scripts/README.md)

## Files

- `manifest_16_fdr.json`
- `same_stage_mark_summary.tsv`
- `stability_by_stage.tsv`
- `stability_by_stage_mark.tsv`
- `stability_by_stage_region.tsv`
