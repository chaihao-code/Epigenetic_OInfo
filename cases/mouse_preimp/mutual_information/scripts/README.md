# cases/mouse_preimp/mutual_information/scripts

Analysis programs or support modules for the parent workflow. O-information arithmetic is provided only by common/02_oinformation/compute_hoi_oinfo.py.

Seven early stages; nine-track full panels and their retained reduced panels. Historical triplet tables use the documented BH84 universe.

## Inputs

The parent case's sample manifests, existing four-column `.bp` tracks (`chr start end signal`) and unchanged genome-specific annotation masks. Comparison/plot programs also read the existing result CSV/TSV files named in their source.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. For the shortest verified workflow, follow `demo/README.md`. For a case command preview, use `python tools/plan_oinfo.py --case <case> --sample <sample>`. Program CLI options are listed below where present.

- `plot_by_block_all_anchors.py`: `python cases/mouse_preimp/mutual_information/scripts/plot_by_block_all_anchors.py --help`; options include `--n-bins`, `--block`, `--anchor`.
- `plot_mi_ratio_hist.py`: `python cases/mouse_preimp/mutual_information/scripts/plot_mi_ratio_hist.py --help`; options include `--n-bins`.

## Outputs

Each script writes the result/figure paths declared in its source. New O-information runs write all-subset entropy, observed Omega/TC/DTC, a permutation summary and metadata; triplet-null runs also write a BH-corrected significance table. Local/sampled-bin exports are enabled with `--save-local-oinfo --save-sampled-bins`.

## Files

- `composite_plot.job.sh`
- `composite_t2.job.sh`
- `plot_by_block_all_anchors.py`
- `plot_composite_trend_hist_n3_n5_n10.py`
- `plot_composite_trend_hist_n3_n5_n10_transposed.py`
- `plot_joint_vs_pairwise_ratio_heatmap_v3.py`
- `plot_mi_n3_n10.job.sh`
- `plot_mi_ratio_hist.py`
- `rerun_heatmap.job.sh`
- `scan_block.job.sh`
- `scan_block_candidates.py`
