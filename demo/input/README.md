# demo/input

Directory index and supporting files for the parent workflow.

## Inputs

See the files/subdirectories below and the parent README. Coordinate inputs use zero-based, half-open genomic intervals; signal lists use repository-root-relative paths. Read `docs/PATHS.md` for external-input restoration.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. Read `demo/input/bigwig_sources.tsv` with pandas (`sep="\t"`; gzip compression is inferred). Existing result tables do not execute a pipeline.

## Outputs

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Subdirectories

- [bigwig](bigwig/README.md)
- [binned](binned/README.md)

## Files

- `bigwig_sources.tsv`
- `hg38_bins_1kb.bed.gz`
- `signals.txt`
