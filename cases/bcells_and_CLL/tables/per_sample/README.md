# cases/bcells_and_CLL/tables/per_sample

Observed O-information tables grouped by biological sample. Each table retains the subsystem orders and exact observed values of its source run.

17 epigenomic donors. RNA male contrast: 4 CLL-R versus 5 CLL-S, adjusting RNA method and IGHV. DMR: 13 CLL donors, M-value limma with age/sex. New per-sample triplet BH covers 1,400 tests; a separate global step covers all 23,800 tests.

## Inputs

The parent case's sample manifests, existing four-column `.bp` tracks (`chr start end signal`) and unchanged genome-specific annotation masks. Comparison/plot programs also read the existing result CSV/TSV files named in their source.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. Read `cases/bcells_and_CLL/tables/per_sample/CEEHRC.CEMT0004.tsv.gz` with pandas (`sep="\t"`; gzip compression is inferred). Existing result tables do not execute a pipeline.

## Outputs

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Files

- `CEEHRC.CEMT0004.tsv.gz`
- `CEEHRC.CEMT0005.tsv.gz`
- `CEEHRC.CEMT0006.tsv.gz`
- `CEEHRC.CEMT0025.tsv.gz`
- `CEEHRC.CEMT0026.tsv.gz`
- `CEEHRC.CEMT0027.tsv.gz`
- `CEEHRC.CEMT0028.tsv.gz`
- `CEEHRC.CEMT0029.tsv.gz`
- `CEEHRC.CEMT0030.tsv.gz`
- `CEEHRC.CEMT0094.tsv.gz`
- `CEEHRC.CEMT0095.tsv.gz`
- `CEEHRC.CEMT0096.tsv.gz`
- `CEEHRC.CEMT0097.tsv.gz`
- `CEEHRC.CEMT0129.tsv.gz`
- `CEEHRC.CEMT0134.tsv.gz`
- `CEEHRC.CEMT0135.tsv.gz`
- `CEEHRC.CEMT0138.tsv.gz`
