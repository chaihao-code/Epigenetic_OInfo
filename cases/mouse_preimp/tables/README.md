# cases/mouse_preimp/tables

Directory index and supporting files for the parent workflow.

Seven early stages; nine-track full panels and their retained reduced panels. Historical triplet tables use the documented BH84 universe.

## Inputs

The parent case's sample manifests, existing four-column `.bp` tracks (`chr start end signal`) and unchanged genome-specific annotation masks. Comparison/plot programs also read the existing result CSV/TSV files named in their source.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. Read `cases/mouse_preimp/tables/per_class_index.tsv` with pandas (`sep="\t"`; gzip compression is inferred). Existing result tables do not execute a pipeline.

## Outputs

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Subdirectories

- [per_class](per_class/README.md)
- [per_sample](per_sample/README.md)

## Files

- `per_class_index.tsv`
- `per_sample_index.tsv`
