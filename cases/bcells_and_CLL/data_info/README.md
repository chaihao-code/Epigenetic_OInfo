# cases/bcells_and_CLL/data_info

Input manifests and retained metadata for the parent workflow. Repository paths in text manifests are rooted at .epigenetic_oinfo_root.

17 epigenomic donors. RNA male contrast: 4 CLL-R versus 5 CLL-S, adjusting RNA method and IGHV. DMR: 13 CLL donors, M-value limma with age/sex. New per-sample triplet BH covers 1,400 tests; a separate global step covers all 23,800 tests.

## Inputs

The parent case's sample manifests, existing four-column `.bp` tracks (`chr start end signal`) and unchanged genome-specific annotation masks. Comparison/plot programs also read the existing result CSV/TSV files named in their source.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. Read `cases/bcells_and_CLL/data_info/cll_metadata.tsv` with pandas (`sep="\t"`; gzip compression is inferred). Existing result tables do not execute a pipeline.

## Outputs

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Subdirectories

- [ihec_rna_lookup](ihec_rna_lookup/README.md)

## Files

- `CEEHRC.CEMT0004.txt`
- `CEEHRC.CEMT0005.txt`
- `CEEHRC.CEMT0006.txt`
- `CEEHRC.CEMT0025.txt`
- `CEEHRC.CEMT0026.txt`
- `CEEHRC.CEMT0027.txt`
- `CEEHRC.CEMT0028.txt`
- `CEEHRC.CEMT0029.txt`
- `CEEHRC.CEMT0030.txt`
- `CEEHRC.CEMT0094.txt`
- `CEEHRC.CEMT0095.txt`
- `CEEHRC.CEMT0096.txt`
- `CEEHRC.CEMT0097.txt`
- `CEEHRC.CEMT0129.txt`
- `CEEHRC.CEMT0134.txt`
- `CEEHRC.CEMT0135.txt`
- `CEEHRC.CEMT0138.txt`
- `cll_metadata.tsv`
- `external_input_provenance.json`
