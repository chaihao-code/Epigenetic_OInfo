# cases/mouse_development/shared/hoi_code

Analysis programs or support modules for the parent workflow. O-information arithmetic is provided only by common/02_oinformation/compute_hoi_oinfo.py.

82 samples, seven signals and 35 triplets per sample/class. New triplet calls use 10,000 permutations and BH across all 35 lower-tail tests in each sample/class.

## Inputs

The parent case's sample manifests, existing four-column `.bp` tracks (`chr start end signal`) and unchanged genome-specific annotation masks. Comparison/plot programs also read the existing result CSV/TSV files named in their source.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. For the shortest verified workflow, follow `demo/README.md`. For a case command preview, use `python tools/plan_oinfo.py --case <case> --sample <sample>`. Program CLI options are listed below where present.

- `block_processor.py`: `python cases/mouse_development/shared/hoi_code/block_processor.py --help`; options include `--signal-matrix`, `--block-files`, `--block-groups`, `--output-dir`, `--n-max`, `--seed`, `--quiet`.
- `compute_hoi_oinfo.py` forwards imports and command-line execution to the single shared engine.
- `compute_redundancy_core.py`: `python cases/mouse_development/shared/hoi_code/compute_redundancy_core.py --help`; options include `--input`, `--output`, `--signals`, `--sample`, `--min-order`, `--z-threshold`, `--eta-threshold`, `--expected-blocks`.
- `compute_signal_correlation.py`: `python cases/mouse_development/shared/hoi_code/compute_signal_correlation.py --help`; options include `--input-file`, `--output`, `--n-bins`, `--method`, `--continuous`, `--block-files`, `--block-groups`, `--group-separator`, `--per-block`, `--n-max`.
- `convert_entropy_to_mi.py`: `python cases/mouse_development/shared/hoi_code/convert_entropy_to_mi.py --help`; options include `--input-dir`, `--output-dir`.
- `signal_loader.py`: `python cases/mouse_development/shared/hoi_code/signal_loader.py --help`; options include `--input`, `--output`, `--n-bins`, `--quiet`.

## Outputs

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Files

- `block_processor.py`
- `compute_hoi_oinfo.py`
- `compute_redundancy_core.py`
- `compute_signal_correlation.py`
- `convert_entropy_to_mi.py`
- `signal_loader.py`
