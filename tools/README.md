# tools

Analysis programs or support modules for the parent workflow. O-information arithmetic is provided only by common/02_oinformation/compute_hoi_oinfo.py.

## Inputs

See the files/subdirectories below and the parent README. Coordinate inputs use zero-based, half-open genomic intervals; signal lists use repository-root-relative paths. Read `docs/PATHS.md` for external-input restoration.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. For the shortest verified workflow, follow `demo/README.md`. For a case command preview, use `python tools/plan_oinfo.py --case <case> --sample <sample>`. Program CLI options are listed below where present.

- `package_release.py`: `python tools/package_release.py --help`; options include `--write`.
- `plan_oinfo.py`: `python tools/plan_oinfo.py --help`; options include `--case`, `--sample`, `--jobs`, `--n-permutations`, `--n-bins`.

## Outputs

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Subdirectories

- [bin](bin/README.md)

## Files

- `package_release.py`
- `plan_oinfo.py`
- `restore_annotation_masks.py`
