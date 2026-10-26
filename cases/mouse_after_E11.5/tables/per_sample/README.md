# cases/mouse_after_E11.5/tables/per_sample

Observed O-information tables grouped by biological sample. Each table retains the subsystem orders and exact observed values of its source run.

66 tissue/stage samples. Retain the distinction between the full panel and the eight-track comparison panel.

## Inputs

The parent case's sample manifests, existing four-column `.bp` tracks (`chr start end signal`) and unchanged genome-specific annotation masks. Comparison/plot programs also read the existing result CSV/TSV files named in their source.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. Read `cases/mouse_after_E11.5/tables/per_sample/E11.5_Embryonic_Facial_Prominence.tsv.gz` with pandas (`sep="\t"`; gzip compression is inferred). Existing result tables do not execute a pipeline.

## Outputs

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Files

- `E11.5_Embryonic_Facial_Prominence.tsv.gz`
- `E11.5_Forebrain.tsv.gz`
- `E11.5_Heart.tsv.gz`
- `E11.5_Hindbrain.tsv.gz`
- `E11.5_Limb.tsv.gz`
- `E11.5_Liver.tsv.gz`
- `E11.5_Midbrain.tsv.gz`
- `E11.5_Neural_Tube.tsv.gz`
- `E12.5_Embryonic_Facial_Prominence.tsv.gz`
- `E12.5_Forebrain.tsv.gz`
- `E12.5_Heart.tsv.gz`
- `E12.5_Hindbrain.tsv.gz`
- `E12.5_Limb.tsv.gz`
- `E12.5_Liver.tsv.gz`
- `E12.5_Midbrain.tsv.gz`
- `E12.5_Neural_Tube.tsv.gz`
- `E13.5_Embryonic_Facial_Prominence.tsv.gz`
- `E13.5_Forebrain.tsv.gz`
- `E13.5_Heart.tsv.gz`
- `E13.5_Hindbrain.tsv.gz`
- `E13.5_Limb.tsv.gz`
- `E13.5_Liver.tsv.gz`
- `E13.5_Midbrain.tsv.gz`
- `E13.5_Neural_Tube.tsv.gz`
- `E14.5_Embryonic_Facial_Prominence.tsv.gz`
- `E14.5_Forebrain.tsv.gz`
- `E14.5_Heart.tsv.gz`
- `E14.5_Hindbrain.tsv.gz`
- `E14.5_Intestine.tsv.gz`
- `E14.5_Kidney.tsv.gz`
- `E14.5_Limb.tsv.gz`
- `E14.5_Liver.tsv.gz`
- `E14.5_Lung.tsv.gz`
- `E14.5_Midbrain.tsv.gz`
- `E14.5_Neural_Tube.tsv.gz`
- `E14.5_Stomach.tsv.gz`
- `E15.5_Embryonic_Facial_Prominence.tsv.gz`
- `E15.5_Forebrain.tsv.gz`
- `E15.5_Heart.tsv.gz`
- `E15.5_Hindbrain.tsv.gz`
- `E15.5_Intestine.tsv.gz`
- `E15.5_Kidney.tsv.gz`
- `E15.5_Limb.tsv.gz`
- `E15.5_Liver.tsv.gz`
- `E15.5_Lung.tsv.gz`

This directory contains 66 files; see the source manifest for the full inventory.
