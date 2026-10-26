# cases/mouse_development/tables/per_sample

Observed O-information tables grouped by biological sample. Each table retains the subsystem orders and exact observed values of its source run.

82 samples, seven signals and 35 triplets per sample/class. New triplet calls use 10,000 permutations and BH across all 35 lower-tail tests in each sample/class.

## Inputs

The parent case's sample manifests, existing four-column `.bp` tracks (`chr start end signal`) and unchanged genome-specific annotation masks. Comparison/plot programs also read the existing result CSV/TSV files named in their source.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. Read `cases/mouse_development/tables/per_sample/2cell.tsv.gz` with pandas (`sep="\t"`; gzip compression is inferred). Existing result tables do not execute a pipeline.

## Outputs

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Files

- `2cell.tsv.gz`
- `4cell.tsv.gz`
- `8cell.tsv.gz`
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

This directory contains 82 files; see the source manifest for the full inventory.
