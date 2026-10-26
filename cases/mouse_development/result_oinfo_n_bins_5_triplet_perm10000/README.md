# cases/mouse_development/result_oinfo_n_bins_5_triplet_perm10000

Retained numerical outputs and/or new explicitly named demo outputs for the parent analysis. The metadata in each run identifies its actual permutation count and signal panel.

82 samples, seven signals and 35 triplets per sample/class. New triplet calls use 10,000 permutations and BH across all 35 lower-tail tests in each sample/class.

## Inputs

The parent case's sample manifests, existing four-column `.bp` tracks (`chr start end signal`) and unchanged genome-specific annotation masks. Comparison/plot programs also read the existing result CSV/TSV files named in their source.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. For the shortest verified workflow, follow `demo/README.md`. For a case command preview, use `python tools/plan_oinfo.py --case <case> --sample <sample>`. Program CLI options are listed below where present.

## Outputs

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Subdirectories

- [2cell](2cell/README.md)
- [4cell](4cell/README.md)
- [8cell](8cell/README.md)
- [E11.5_Embryonic_Facial_Prominence](E11.5_Embryonic_Facial_Prominence/README.md)
- [E11.5_Forebrain](E11.5_Forebrain/README.md)
- [E11.5_Heart](E11.5_Heart/README.md)
- [E11.5_Hindbrain](E11.5_Hindbrain/README.md)
- [E11.5_Limb](E11.5_Limb/README.md)
- [E11.5_Liver](E11.5_Liver/README.md)
- [E11.5_Midbrain](E11.5_Midbrain/README.md)
- [E11.5_Neural_Tube](E11.5_Neural_Tube/README.md)
- [E12.5_Embryonic_Facial_Prominence](E12.5_Embryonic_Facial_Prominence/README.md)
- [E12.5_Forebrain](E12.5_Forebrain/README.md)
- [E12.5_Heart](E12.5_Heart/README.md)
- [E12.5_Hindbrain](E12.5_Hindbrain/README.md)
- [E12.5_Limb](E12.5_Limb/README.md)
- [E12.5_Liver](E12.5_Liver/README.md)
- [E12.5_Midbrain](E12.5_Midbrain/README.md)
- [E12.5_Neural_Tube](E12.5_Neural_Tube/README.md)
- [E13.5_Embryonic_Facial_Prominence](E13.5_Embryonic_Facial_Prominence/README.md)
- [E13.5_Forebrain](E13.5_Forebrain/README.md)
- [E13.5_Heart](E13.5_Heart/README.md)
- [E13.5_Hindbrain](E13.5_Hindbrain/README.md)
- [E13.5_Limb](E13.5_Limb/README.md)
- [E13.5_Liver](E13.5_Liver/README.md)
- [E13.5_Midbrain](E13.5_Midbrain/README.md)
- [E13.5_Neural_Tube](E13.5_Neural_Tube/README.md)
- [E14.5_Embryonic_Facial_Prominence](E14.5_Embryonic_Facial_Prominence/README.md)
- [E14.5_Forebrain](E14.5_Forebrain/README.md)
- [E14.5_Heart](E14.5_Heart/README.md)
- [E14.5_Hindbrain](E14.5_Hindbrain/README.md)
- [E14.5_Intestine](E14.5_Intestine/README.md)
- [E14.5_Kidney](E14.5_Kidney/README.md)
- [E14.5_Limb](E14.5_Limb/README.md)
- [E14.5_Liver](E14.5_Liver/README.md)
- [E14.5_Lung](E14.5_Lung/README.md)
- [E14.5_Midbrain](E14.5_Midbrain/README.md)
- [E14.5_Neural_Tube](E14.5_Neural_Tube/README.md)
- [E14.5_Stomach](E14.5_Stomach/README.md)
- [E15.5_Embryonic_Facial_Prominence](E15.5_Embryonic_Facial_Prominence/README.md)
- [E15.5_Forebrain](E15.5_Forebrain/README.md)
- [E15.5_Heart](E15.5_Heart/README.md)
- [E15.5_Hindbrain](E15.5_Hindbrain/README.md)
- [E15.5_Intestine](E15.5_Intestine/README.md)
- [E15.5_Kidney](E15.5_Kidney/README.md)
- [E15.5_Limb](E15.5_Limb/README.md)
- [E15.5_Liver](E15.5_Liver/README.md)
- [E15.5_Lung](E15.5_Lung/README.md)
- [E15.5_Midbrain](E15.5_Midbrain/README.md)
- [E15.5_Neural_Tube](E15.5_Neural_Tube/README.md)
- [E15.5_Stomach](E15.5_Stomach/README.md)
- [E16.5_Forebrain](E16.5_Forebrain/README.md)
- [E16.5_Heart](E16.5_Heart/README.md)
- [E16.5_Hindbrain](E16.5_Hindbrain/README.md)
- [E16.5_Intestine](E16.5_Intestine/README.md)
- [E16.5_Kidney](E16.5_Kidney/README.md)
- [E16.5_Liver](E16.5_Liver/README.md)
- [E16.5_Lung](E16.5_Lung/README.md)
- [E16.5_Midbrain](E16.5_Midbrain/README.md)
- [E16.5_Stomach](E16.5_Stomach/README.md)
- [E6.5_Epi](E6.5_Epi/README.md)
- [E6.5_ExE](E6.5_ExE/README.md)
- [E7.0_A](E7.0_A/README.md)
- [E7.0_ExE](E7.0_ExE/README.md)
- [E7.0_P](E7.0_P/README.md)
- [E7.5_A](E7.5_A/README.md)
- [E7.5_AM](E7.5_AM/README.md)
- [E7.5_ExE](E7.5_ExE/README.md)
- [E7.5_P](E7.5_P/README.md)
- [ICM](ICM/README.md)
- [P0_Forebrain](P0_Forebrain/README.md)
- [P0_Heart](P0_Heart/README.md)
- [P0_Hindbrain](P0_Hindbrain/README.md)
- [P0_Intestine](P0_Intestine/README.md)
- [P0_Kidney](P0_Kidney/README.md)
- [P0_Liver](P0_Liver/README.md)
- [P0_Lung](P0_Lung/README.md)
- [P0_Midbrain](P0_Midbrain/README.md)
- [P0_Stomach](P0_Stomach/README.md)
- [TE](TE/README.md)
- [code](code/README.md)
- [morula](morula/README.md)
- [negative_fdr](negative_fdr/README.md)
- [zygote](zygote/README.md)

## Files

- `code_provenance.json`
- `input_manifest.tsv`
- `run_settings.json`
- `samples.txt`
- `submit.py`
- `submit_array.sh`
- `submit_correction.sh`
- `validation.json`
