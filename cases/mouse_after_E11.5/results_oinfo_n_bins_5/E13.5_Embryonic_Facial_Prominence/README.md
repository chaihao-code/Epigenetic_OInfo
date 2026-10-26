# cases/mouse_after_E11.5/results_oinfo_n_bins_5/E13.5_Embryonic_Facial_Prominence

Retained numerical outputs and/or new explicitly named demo outputs for the parent analysis. The metadata in each run identifies its actual permutation count and signal panel.

66 tissue/stage samples. Retain the distinction between the full panel and the eight-track comparison panel.

## Inputs

The parent case's sample manifests, existing four-column `.bp` tracks (`chr start end signal`) and unchanged genome-specific annotation masks. Comparison/plot programs also read the existing result CSV/TSV files named in their source.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. Read `cases/mouse_after_E11.5/results_oinfo_n_bins_5/E13.5_Embryonic_Facial_Prominence/E13.5_Embryonic_Facial_Prominence_oinfo.csv` with pandas (`sep=","`; gzip compression is inferred). Existing result tables do not execute a pipeline.

## Outputs

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Files

- `E13.5_Embryonic_Facial_Prominence_oinfo.csv`
- `E13.5_Embryonic_Facial_Prominence_permutation_summary.csv`
