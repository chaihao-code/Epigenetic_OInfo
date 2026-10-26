# cases/aging/neutrophil/data_info

Input manifests and retained metadata for the parent workflow. Repository paths in text manifests are rooted at .epigenetic_oinfo_root.

12 retained samples; original O-information and age-regression reference tables. This cleanup does not invent a new neutrophil triplet study.

## Inputs

The parent case's sample manifests, existing four-column `.bp` tracks (`chr start end signal`) and unchanged genome-specific annotation masks. Comparison/plot programs also read the existing result CSV/TSV files named in their source.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. For the shortest verified workflow, follow `demo/README.md`. For a case command preview, use `python tools/plan_oinfo.py --case <case> --sample <sample>`. Program CLI options are listed below where present.

## Outputs

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Files

- `IHECRE00000004.txt`
- `IHECRE00000094.txt`
- `IHECRE00000095.txt`
- `IHECRE00000109.txt`
- `IHECRE00000124.txt`
- `IHECRE00000159.txt`
- `IHECRE00000256.txt`
- `IHECRE00000265.txt`
- `IHECRE00000284.txt`
- `IHECRE00000308.txt`
- `IHECRE00000346.txt`
- `IHECRE00000353.txt`
