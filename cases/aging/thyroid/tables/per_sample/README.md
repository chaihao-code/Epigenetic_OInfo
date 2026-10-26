# cases/aging/thyroid/tables/per_sample

Observed O-information tables grouped by biological sample. Each table retains the subsystem orders and exact observed values of its source run.

Six healthy donors: CEMT0098, CEMT0130, CEMT0132, CEMT0042, CEMT0044 and CEMT0040. Seven signals, hg38, the original fixed 1-kb grid and 40 annotation classes. The complete CEMT0098 donor is the release demo.

## Inputs

The original hg38 1-kb grid, seven binned tracks per donor, unchanged genomic masks and retained donor ages/cell fractions. The complete input-to-output example is `demo/README.md`; raw/processed local input locations are listed in `config/local_input_links.tsv`.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. Read `cases/aging/thyroid/tables/per_sample/CEMT0040.tsv.gz` with pandas (`sep="\t"`; gzip compression is inferred). Existing result tables do not execute a pipeline.

## Outputs

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Files

- `CEMT0040.tsv.gz`
- `CEMT0042.tsv.gz`
- `CEMT0044.tsv.gz`
- `CEMT0098.tsv.gz`
- `CEMT0130.tsv.gz`
- `CEMT0132.tsv.gz`
