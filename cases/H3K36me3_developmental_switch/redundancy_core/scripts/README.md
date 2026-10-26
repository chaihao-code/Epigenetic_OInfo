# cases/H3K36me3_developmental_switch/redundancy_core/scripts

Analysis programs or support modules for the parent workflow. O-information arithmetic is provided only by common/02_oinformation/compute_hoi_oinfo.py.

The eight-signal BH56 comparison reuses preimplantation/late-development observations. Full-panel occurrence uses the retained sign-only rule. Replacement-background redundancy cores have a different null from permutation-tested triplets.

## Inputs

The parent case's sample manifests, existing four-column `.bp` tracks (`chr start end signal`) and unchanged genome-specific annotation masks. Comparison/plot programs also read the existing result CSV/TSV files named in their source.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. For the shortest verified workflow, follow `demo/README.md`. For a case command preview, use `python tools/plan_oinfo.py --case <case> --sample <sample>`. Program CLI options are listed below where present.

- `plot_full_signal_combination_occurrence.py`: `python cases/H3K36me3_developmental_switch/redundancy_core/scripts/plot_full_signal_combination_occurrence.py --help`; options include `--prepare`, `--core-task`, `--plot`.

## Outputs

Each script writes the result/figure paths declared in its source. New O-information runs write all-subset entropy, observed Omega/TC/DTC, a permutation summary and metadata; triplet-null runs also write a BH-corrected significance table. Local/sampled-bin exports are enabled with `--save-local-oinfo --save-sampled-bins`.

## Files

- `compute_full_signal_occurrence_cores.slurm`
- `compute_redundancy_cores.slurm`
- `finalize_analysis.py`
- `plot_combined_trajectories.py`
- `plot_core_occurrence_frequency_e11.5.py`
- `plot_core_occurrence_frequency_preimp.py`
- `plot_full_signal_combination_occurrence.py`
- `plot_redundancy_core_trajectories.py`
- `plot_synergistic_triplet_trajectories.py`
- `plot_trajectories.slurm`
- `plot_triplet_combination_occurrence_with_H3K9ac_H3K4me2.py`
- `run_core_task.py`
- `trajectory_utils.py`

## Redrawing retained figures

`plot_redundancy_core_trajectories.py --plots-only` and
`plot_synergistic_triplet_trajectories.py --plots-only` read the retained
`results/*_mark_counts_by_sample.tsv` tables, without recollecting per-sample cores.
`plot_full_signal_combination_occurrence.py --plots-only` reads the retained
`full_signal_occurrence/data/*_combination_occurrence_*_top25.tsv` and summary.
These modes do not recompute or revalidate the original statistical analyses.

`plot_triplet_combination_occurrence_with_H3K9ac_H3K4me2.py A --plots-only`
(and B/C) uses the retained full-panel post-E11.5 sign-only combination counts.
The classification is unchanged. The former default `D` had no definition and
always produced no figure; the default is now the valid class `A`.

The historical `plot_core_occurrence_frequency_preimp.py` and
`plot_core_occurrence_frequency_e11.5.py` require their separate
`all_supreme_cores_filtered.csv` inputs, which are absent in this release.
They are not substituted with results from another core analysis.
