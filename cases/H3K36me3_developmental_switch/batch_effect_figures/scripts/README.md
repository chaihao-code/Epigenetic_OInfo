# cases/H3K36me3_developmental_switch/batch_effect_figures/scripts

Analysis programs or support modules for the parent workflow. O-information arithmetic is provided only by common/02_oinformation/compute_hoi_oinfo.py.

The eight-signal BH56 comparison reuses preimplantation/late-development observations. Full-panel occurrence uses the retained sign-only rule. Replacement-background redundancy cores have a different null from permutation-tested triplets.

## Inputs

The parent case's sample manifests, existing four-column `.bp` tracks (`chr start end signal`) and unchanged genome-specific annotation masks. Comparison/plot programs also read the existing result CSV/TSV files named in their source.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. For the shortest verified workflow, follow `demo/README.md`. For a case command preview, use `python tools/plan_oinfo.py --case <case> --sample <sample>`. Program CLI options are listed below where present.

- `plot_batch_effect.py`: `python cases/H3K36me3_developmental_switch/batch_effect_figures/scripts/plot_batch_effect.py --help`; options include `--validate-inputs`, `--plots-only`.

## Outputs

Each script writes the result/figure paths declared in its source. New O-information runs write all-subset entropy, observed Omega/TC/DTC, a permutation summary and metadata; triplet-null runs also write a BH-corrected significance table. Local/sampled-bin exports are enabled with `--save-local-oinfo --save-sampled-bins`.

## Files

- `compare_stability_fdr_thresholds.py`
- `compute_bulk16.slurm`
- `intersection_fraction_of_each_set.py`
- `plot.slurm`
- `plot_batch_effect.py`
- `plot_support_without_h3k27me3_replacement.py`
- `plot_triplet_support_distribution.py`
