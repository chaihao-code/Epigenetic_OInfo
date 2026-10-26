# cases/bcells_and_CLL/dmr_repro_cemt/scripts

Analysis programs or support modules for the parent workflow. O-information arithmetic is provided only by common/02_oinformation/compute_hoi_oinfo.py.

17 epigenomic donors. RNA male contrast: 4 CLL-R versus 5 CLL-S, adjusting RNA method and IGHV. DMR: 13 CLL donors, M-value limma with age/sex. New per-sample triplet BH covers 1,400 tests; a separate global step covers all 23,800 tests.

## Inputs

Prepared four-coordinate/sample methylation matrix at `input/all_bins_1kb_cemt_g2_g3.csv`, the retained 13-donor age/sex metadata and optional GENCODE annotation. Beta values are on the original 0–1 scale.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. For the shortest verified workflow, follow `demo/README.md`. For a case command preview, use `python tools/plan_oinfo.py --case <case> --sample <sample>`. Program CLI options are listed below where present.

## Outputs

The canonical case wrapper writes limma results, BH q-values, delta-beta classifications and run diagnostics to its `output/` directory. Adjacent 1-kb windows are not merged.

## Files

- `check_input.py`
- `dmr_limma_1kb_cemt_g2_g3.R`
- `install_r_packages.R`
- `plot_dmr_volcano.py`
- `plot_igh_tracks_corrected.py`
- `run_dmr.R`
- `submit_dmr.slurm`
