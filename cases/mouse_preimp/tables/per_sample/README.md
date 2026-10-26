# cases/mouse_preimp/tables/per_sample

Observed O-information tables grouped by biological sample. Each table retains the subsystem orders and exact observed values of its source run.

Seven early stages; nine-track full panels and their retained reduced panels. Historical triplet tables use the documented BH84 universe.

## Inputs

The parent case's sample manifests, existing four-column `.bp` tracks (`chr start end signal`) and unchanged genome-specific annotation masks. Comparison/plot programs also read the existing result CSV/TSV files named in their source.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. Read `cases/mouse_preimp/tables/per_sample/2cell.tsv.gz` with pandas (`sep="\t"`; gzip compression is inferred). Existing result tables do not execute a pipeline.

## Outputs

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Files

- `2cell.tsv.gz`
- `4cell.tsv.gz`
- `8cell.tsv.gz`
- `ICM.tsv.gz`
- `TE.tsv.gz`
- `morula.tsv.gz`
- `zygote.tsv.gz`
