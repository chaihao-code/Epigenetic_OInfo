# cases/H3K36me3_developmental_switch/batch_effect_figures/inputs

Input manifests and retained metadata for the parent workflow. Repository paths in text manifests are rooted at .epigenetic_oinfo_root.

The eight-signal BH56 comparison reuses preimplantation/late-development observations. Full-panel occurrence uses the retained sign-only rule. Replacement-background redundancy cores have a different null from permutation-tested triplets.

## Inputs

The parent case's sample manifests, existing four-column `.bp` tracks (`chr start end signal`) and unchanged genome-specific annotation masks. Comparison/plot programs also read the existing result CSV/TSV files named in their source.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. Read `cases/H3K36me3_developmental_switch/batch_effect_figures/inputs/tasks.tsv` with pandas (`sep="\t"`; gzip compression is inferred). Existing result tables do not execute a pipeline.

## Outputs

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Subdirectories

- [signals](signals/README.md)

## Files

- `2cell_bulk00.txt`
- `2cell_bulk01.txt`
- `2cell_bulk02.txt`
- `2cell_bulk03.txt`
- `2cell_bulk04.txt`
- `2cell_bulk05.txt`
- `2cell_bulk06.txt`
- `2cell_bulk07.txt`
- `2cell_bulk08.txt`
- `2cell_bulk09.txt`
- `2cell_bulk10.txt`
- `2cell_bulk11.txt`
- `2cell_bulk12.txt`
- `2cell_bulk13.txt`
- `2cell_bulk14.txt`
- `2cell_bulk15.txt`
- `8cell_bulk00.txt`
- `8cell_bulk01.txt`
- `8cell_bulk02.txt`
- `8cell_bulk03.txt`
- `8cell_bulk04.txt`
- `8cell_bulk05.txt`
- `8cell_bulk06.txt`
- `8cell_bulk07.txt`
- `8cell_bulk08.txt`
- `8cell_bulk09.txt`
- `8cell_bulk10.txt`
- `8cell_bulk11.txt`
- `8cell_bulk12.txt`
- `8cell_bulk13.txt`
- `8cell_bulk14.txt`
- `8cell_bulk15.txt`
- `ICM_bulk00.txt`
- `ICM_bulk01.txt`
- `ICM_bulk02.txt`
- `ICM_bulk03.txt`
- `ICM_bulk04.txt`
- `ICM_bulk05.txt`
- `ICM_bulk06.txt`
- `ICM_bulk07.txt`
- `ICM_bulk08.txt`
- `ICM_bulk09.txt`
- `ICM_bulk10.txt`
- `ICM_bulk11.txt`
- `ICM_bulk12.txt`

This directory contains 49 files; see the source manifest for the full inventory.
