# cases/mouse_development/data/input_files

Directory index and supporting files for the parent workflow.

82 samples, seven signals and 35 triplets per sample/class. New triplet calls use 10,000 permutations and BH across all 35 lower-tail tests in each sample/class.

## Inputs

The parent case's sample manifests, existing four-column `.bp` tracks (`chr start end signal`) and unchanged genome-specific annotation masks. Comparison/plot programs also read the existing result CSV/TSV files named in their source.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. For the shortest verified workflow, follow `demo/README.md`. For a case command preview, use `python tools/plan_oinfo.py --case <case> --sample <sample>`. Program CLI options are listed below where present.

## Outputs

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Files

- `2cell.txt`
- `4cell.txt`
- `8cell.txt`
- `E11.5_Embryonic_Facial_Prominence.txt`
- `E11.5_Forebrain.txt`
- `E11.5_Heart.txt`
- `E11.5_Hindbrain.txt`
- `E11.5_Limb.txt`
- `E11.5_Liver.txt`
- `E11.5_Midbrain.txt`
- `E11.5_Neural_Tube.txt`
- `E12.5_Embryonic_Facial_Prominence.txt`
- `E12.5_Forebrain.txt`
- `E12.5_Heart.txt`
- `E12.5_Hindbrain.txt`
- `E12.5_Limb.txt`
- `E12.5_Liver.txt`
- `E12.5_Midbrain.txt`
- `E12.5_Neural_Tube.txt`
- `E13.5_Embryonic_Facial_Prominence.txt`
- `E13.5_Forebrain.txt`
- `E13.5_Heart.txt`
- `E13.5_Hindbrain.txt`
- `E13.5_Limb.txt`
- `E13.5_Liver.txt`
- `E13.5_Midbrain.txt`
- `E13.5_Neural_Tube.txt`
- `E14.5_Embryonic_Facial_Prominence.txt`
- `E14.5_Forebrain.txt`
- `E14.5_Heart.txt`
- `E14.5_Hindbrain.txt`
- `E14.5_Intestine.txt`
- `E14.5_Kidney.txt`
- `E14.5_Limb.txt`
- `E14.5_Liver.txt`
- `E14.5_Lung.txt`
- `E14.5_Midbrain.txt`
- `E14.5_Neural_Tube.txt`
- `E14.5_Stomach.txt`
- `E15.5_Embryonic_Facial_Prominence.txt`
- `E15.5_Forebrain.txt`
- `E15.5_Heart.txt`
- `E15.5_Hindbrain.txt`
- `E15.5_Intestine.txt`
- `E15.5_Kidney.txt`

This directory contains 82 files; see the source manifest for the full inventory.
