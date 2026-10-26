# cases/H3K36me3_developmental_switch/synergistic_triplets/scripts

Analysis programs or support modules for the parent workflow. O-information arithmetic is provided only by common/02_oinformation/compute_hoi_oinfo.py.

The eight-signal BH56 comparison reuses preimplantation/late-development observations. Full-panel occurrence uses the retained sign-only rule. Replacement-background redundancy cores have a different null from permutation-tested triplets.

## Inputs

The parent case's sample manifests, existing four-column `.bp` tracks (`chr start end signal`) and unchanged genome-specific annotation masks. Comparison/plot programs also read the existing result CSV/TSV files named in their source.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. For the shortest verified workflow, follow `demo/README.md`. For a case command preview, use `python tools/plan_oinfo.py --case <case> --sample <sample>`. Program CLI options are listed below where present.

## Outputs

Each script writes the result/figure paths declared in its source. New O-information runs write all-subset entropy, observed Omega/TC/DTC, a permutation summary and metadata; triplet-null runs also write a BH-corrected significance table. Local/sampled-bin exports are enabled with `--save-local-oinfo --save-sampled-bins`.

## Files

- `plot_mark_occurrence_per_tissue_by_stage.py`
- `plot_original_development_heatmap.py`
- `plot_triplet_ABC_selected_regions_post_E11p5.py`
- `plot_triplet_ABC_stacked_post_E11p5.py`
- `plot_triplet_ABC_stacked_preimplantation.py`
- `plot_triplet_A_B_C_post_E11p5.py`
