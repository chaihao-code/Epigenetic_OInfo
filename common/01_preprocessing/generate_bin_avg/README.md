# common/01_preprocessing/generate_bin_avg

Analysis programs or support modules for the parent workflow. O-information arithmetic is provided only by common/02_oinformation/compute_hoi_oinfo.py.

## Inputs

See the files/subdirectories below and the parent README. Coordinate inputs use zero-based, half-open genomic intervals; signal lists use repository-root-relative paths. Read `docs/PATHS.md` for external-input restoration.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. For the shortest verified workflow, follow `demo/README.md`. For a case command preview, use `python tools/plan_oinfo.py --case <case> --sample <sample>`. Program CLI options are listed below where present.

## Outputs

Files listed below are already available results/inputs. New calculations are written only when their documented entrypoint is explicitly executed. Historical 1,000-permutation results remain historical reference data; all new triplet workflows default to 10,000.

## Files

- `batch_compute_bin_average.sh`
- `batch_convert_cpg_to_bedgraph.sh`
- `compute_bin_average.sh`
- `compute_weighted_methylation_average.sh`
- `convert_cpg_to_bedgraph.sh`
- `example_mouse_atac_binning.sh`
- `norm_using_fc_to_input.sh`
- `slurm_compute_bin_average.sh`
- `slurm_convert_cpg_to_bedgraph.sh`
