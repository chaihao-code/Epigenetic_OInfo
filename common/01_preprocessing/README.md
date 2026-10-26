# common/01_preprocessing

Analysis programs or support modules for the parent workflow. O-information arithmetic is provided only by common/02_oinformation/compute_hoi_oinfo.py.

## Inputs

See the files/subdirectories below and the parent README. Coordinate inputs use zero-based, half-open genomic intervals; signal lists use repository-root-relative paths. Read `docs/PATHS.md` for external-input restoration.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. For the shortest verified workflow, follow `demo/README.md`. For a case command preview, use `python tools/plan_oinfo.py --case <case> --sample <sample>`. Program CLI options are listed below where present.

- `bin_bigwig.py`: `python common/01_preprocessing/bin_bigwig.py --help`; options include `--input`, `--regions`, `--output`, `--bin-size`.
- `discretize_signal.py`: `python common/01_preprocessing/discretize_signal.py --help`; options include `--input`, `--bins`, `--output`, `--n-bins`, `--keep-autosomes`.

## Outputs

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Subdirectories

- [average_replicates_and_liftover](average_replicates_and_liftover/README.md)
- [generate_bin_avg](generate_bin_avg/README.md)
- [genome_annotation_pipeline](genome_annotation_pipeline/README.md)

## Files

- `bin_bigwig.py`
- `discretize_signal.py`
