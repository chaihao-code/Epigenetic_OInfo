# cases/mouse_preimp/calculation_results

Directory index and supporting files for the parent workflow.

Seven early stages; nine-track full panels and their retained reduced panels. Historical triplet tables use the documented BH84 universe.

## Inputs

The parent case's sample manifests, existing four-column `.bp` tracks (`chr start end signal`) and unchanged genome-specific annotation masks. Comparison/plot programs also read the existing result CSV/TSV files named in their source.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. For the shortest verified workflow, follow `demo/README.md`. For a case command preview, use `python tools/plan_oinfo.py --case <case> --sample <sample>`. Program CLI options are listed below where present.

## Outputs

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Subdirectories

- [data_info](data_info/README.md)
- [results_entropy_n_bins_10](results_entropy_n_bins_10/README.md)
- [results_entropy_n_bins_3](results_entropy_n_bins_3/README.md)
- [results_entropy_n_bins_5](results_entropy_n_bins_5/README.md)
- [scripts](scripts/README.md)

