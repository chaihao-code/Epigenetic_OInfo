# cases/bcells_and_CLL/scripts

Analysis programs or support modules for the parent workflow. O-information arithmetic is provided only by common/02_oinformation/compute_hoi_oinfo.py.

17 epigenomic donors. RNA male contrast: 4 CLL-R versus 5 CLL-S, adjusting RNA method and IGHV. DMR: 13 CLL donors, M-value limma with age/sex. New per-sample triplet BH covers 1,400 tests; a separate global step covers all 23,800 tests.

## Inputs

Use the retained sample metadata, imported counts/TPM in the B-cell case and the frozen GO/GOA annotation snapshot under `data/annotations/go_human/`. O-information programs instead use that case's seven-track signal manifests and original hg38 masks.

## Usage

Run commands from the repository root and set `EPINFO_ROOT` to that root. For the shortest verified workflow, follow `demo/README.md`. For a case command preview, use `python tools/plan_oinfo.py --case <case> --sample <sample>`. Program CLI options are listed below where present.

- `summarize_c_triplet_recurrence.py`: `python cases/bcells_and_CLL/scripts/summarize_c_triplet_recurrence.py --help`; options include `--output-dir`.
- `summarize_triplet_significance.py`: `python cases/bcells_and_CLL/scripts/summarize_triplet_significance.py --help`; options include `--sample`.
- `validate_male_de_in_females.py`: `python cases/bcells_and_CLL/scripts/validate_male_de_in_females.py --help`; options include `--no-plots`.

## Outputs

Each script writes the result/figure paths declared in its source. New O-information runs write all-subset entropy, observed Omega/TC/DTC, a permutation summary and metadata; triplet-null runs also write a BH-corrected significance table. Local/sampled-bin exports are enabled with `--save-local-oinfo --save-sampled-bins`.

## Subdirectories

- [hoi_code](hoi_code/README.md)

## Files

- `aggregate_oinfo_triplets.slurm`
- `analyze_rna_de_male.py`
- `gsea_rna_male.py`
- `import_rna_rsem.R`
- `plot_abc_pooled_combined.py`
- `plot_abc_pooled_compact_preview.py`
- `plot_abc_triplet_classes.py`
- `plot_de_male_ighv_volcano_clear.py`
- `plot_gsea_eight_panel_v2_18x10cm.py`
- `plot_gsea_representative_bars.py`
- `plot_gsea_top20_bars.py`
- `resolve_signal_list.py`
- `run_final_rna.slurm`
- `submit_final_rna.sh`
- `submit_oinfo_hoi.slurm`
- `submit_oinfo_triplets.slurm`
- `summarize_c_triplet_recurrence.py`
- `summarize_triplet_significance.py`
- `validate_male_de_in_females.py`
