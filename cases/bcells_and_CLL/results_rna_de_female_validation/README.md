# cases/bcells_and_CLL/results_rna_de_female_validation

Retained numerical outputs and/or new explicitly named demo outputs for the parent analysis. The metadata in each run identifies its actual permutation count and signal panel.

17 epigenomic donors. RNA male contrast: 4 CLL-R versus 5 CLL-S, adjusting RNA method and IGHV. DMR: 13 CLL donors, M-value limma with age/sex. New per-sample triplet BH covers 1,400 tests; a separate global step covers all 23,800 tests.

## Inputs

Use the retained sample metadata, imported counts/TPM in the B-cell case and the frozen GO/GOA annotation snapshot under `data/annotations/go_human/`. O-information programs instead use that case's seven-track signal manifests and original hg38 masks.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. For the shortest verified workflow, follow `demo/README.md`. For a case command preview, use `python tools/plan_oinfo.py --case <case> --sample <sample>`. Program CLI options are listed below where present.

## Outputs

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Files

- `README_zh.md`
