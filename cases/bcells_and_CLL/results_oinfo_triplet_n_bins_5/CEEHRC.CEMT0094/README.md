# cases/bcells_and_CLL/results_oinfo_triplet_n_bins_5/CEEHRC.CEMT0094

Retained numerical outputs and/or new explicitly named demo outputs for the parent analysis. The metadata in each run identifies its actual permutation count and signal panel.

17 epigenomic donors. RNA male contrast: 4 CLL-R versus 5 CLL-S, adjusting RNA method and IGHV. DMR: 13 CLL donors, M-value limma with age/sex. New per-sample triplet BH covers 1,400 tests; a separate global step covers all 23,800 tests.

## Inputs

The parent case's sample manifests, existing four-column `.bp` tracks (`chr start end signal`) and unchanged genome-specific annotation masks. Comparison/plot programs also read the existing result CSV/TSV files named in their source.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. Read `cases/bcells_and_CLL/results_oinfo_triplet_n_bins_5/CEEHRC.CEMT0094/CEEHRC.CEMT0094_blocks.csv` with pandas (`sep=","`; gzip compression is inferred). Existing result tables do not execute a pipeline.

## Outputs

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Files

- `CEEHRC.CEMT0094_blocks.csv`
- `CEEHRC.CEMT0094_null_summary.csv`
- `CEEHRC.CEMT0094_oinfo.csv`
- `CEEHRC.CEMT0094_subsets.csv`
