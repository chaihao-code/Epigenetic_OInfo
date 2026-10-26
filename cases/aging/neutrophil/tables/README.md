# cases/aging/neutrophil/tables

Directory index and supporting files for the parent workflow.

12 retained samples; original O-information and age-regression reference tables. This cleanup does not invent a new neutrophil triplet study.

## Inputs

The parent case's sample manifests, existing four-column `.bp` tracks (`chr start end signal`) and unchanged genome-specific annotation masks. Comparison/plot programs also read the existing result CSV/TSV files named in their source.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. Read `cases/aging/neutrophil/tables/per_class_index.tsv` with pandas (`sep="\t"`; gzip compression is inferred). Existing result tables do not execute a pipeline.

## Outputs

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Subdirectories

- [per_class](per_class/README.md)
- [per_sample](per_sample/README.md)

## Files

- `per_class_index.tsv`
- `per_sample_index.tsv`

## Bin counts

[bins_per_class.tsv](bins_per_class.tsv) gives the 40 genome-specific reference annotation-class membership counts. These are counts before sample-specific filtering/downsampling, not the loci used in every sample. Original class identifiers are retained; source-specific aliases in result files must be interpreted using the original region definitions.
