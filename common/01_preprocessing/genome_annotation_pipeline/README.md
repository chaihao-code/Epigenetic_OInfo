# common/01_preprocessing/genome_annotation_pipeline

Analysis programs or support modules for the parent workflow. O-information arithmetic is provided only by common/02_oinformation/compute_hoi_oinfo.py.

## Inputs

See the files/subdirectories below and the parent README. Coordinate inputs use zero-based, half-open genomic intervals; signal lists use repository-root-relative paths. Read `docs/PATHS.md` for external-input restoration.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. For the shortest verified workflow, follow `demo/README.md`. For a case command preview, use `python tools/plan_oinfo.py --case <case> --sample <sample>`. Program CLI options are listed below where present.

- `06_purify_masks.py`: `python common/01_preprocessing/genome_annotation_pipeline/06_purify_masks.py --help`; options include `--genome`, `--dry-run`.

## Outputs

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Files

- `00_download_annotations.sh`
- `01_make_bins.sh`
- `02_extract_gene_features.sh`
- `03_ctcf_split.sh`
- `04_bed_to_bins.py`
- `05_run_bed_to_bins.sh`
- `06_purify_masks.py`
- `config.sh`
