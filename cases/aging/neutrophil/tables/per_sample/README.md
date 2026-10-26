# cases/aging/neutrophil/tables/per_sample

Observed O-information tables grouped by biological sample. Each table retains the subsystem orders and exact observed values of its source run.

12 retained samples; original O-information and age-regression reference tables. This cleanup does not invent a new neutrophil triplet study.

## Inputs

The parent case's sample manifests, existing four-column `.bp` tracks (`chr start end signal`) and unchanged genome-specific annotation masks. Comparison/plot programs also read the existing result CSV/TSV files named in their source.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. Read `cases/aging/neutrophil/tables/per_sample/IHECRE00000004.tsv.gz` with pandas (`sep="\t"`; gzip compression is inferred). Existing result tables do not execute a pipeline.

## Outputs

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Files

- `IHECRE00000004.tsv.gz`
- `IHECRE00000094.tsv.gz`
- `IHECRE00000095.tsv.gz`
- `IHECRE00000109.tsv.gz`
- `IHECRE00000124.tsv.gz`
- `IHECRE00000159.tsv.gz`
- `IHECRE00000256.tsv.gz`
- `IHECRE00000265.tsv.gz`
- `IHECRE00000284.tsv.gz`
- `IHECRE00000308.tsv.gz`
- `IHECRE00000346.tsv.gz`
- `IHECRE00000353.tsv.gz`
