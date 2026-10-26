# Retained results and their source code

Paths are relative to this B-cell project. Scripts remain in their existing locations so that imports
and historical result paths continue to work. The source export preserves this layout.

| Result directory | Required generation code | Main inputs |
|---|---|---|
| `results_oinfo_n_bins_5/` | `scripts/submit_oinfo_hoi.slurm`, `scripts/resolve_signal_list.py`, four modules in `scripts/hoi_code/` | 17 signal lists, 119 signal files, hg38 block annotations |
| `results_oinfo_triplet_n_bins_5/` | `scripts/submit_oinfo_triplets.slurm`, same HOI/resolver modules, `scripts/summarize_triplet_significance.py`, `scripts/aggregate_oinfo_triplets.slurm` | Same signals/annotations; per-sample triplet null summaries |
| `ABC_triplet_classes/` | `scripts/plot_abc_triplet_classes.py`, `scripts/plot_abc_pooled_combined.py`, `scripts/plot_abc_pooled_compact_preview.py` | Combined triplet significance table, clinical and RNA sample metadata |
| `ABC_triplet_classes/C_pattern/` | `scripts/summarize_c_triplet_recurrence.py` | ABC classified tests and sample labels |
| `global_oinfo_and_triplets/` | `global_oinfo_and_triplets/plot_global_oinfo_and_ABC_triplets.py`, `reference_config.json`; imports the two main ABC plotting modules | Order-7 O-information, system permutation summaries, ABC pooled FDR counts and figures |
| `results_rna_de_slurm/` | `scripts/import_rna_rsem.R` | RSEM isoform files and clinical sample IDs; curated metadata/annotations are retained inputs |
| `results_rna_de_male_ighv/` | `scripts/analyze_rna_de_male.py`, `scripts/plot_de_male_ighv_volcano_clear.py` | Imported counts/TPM, sample metadata and gene annotation |
| `results_gsea_rna_male/` | `scripts/gsea_rna_male.py`, `scripts/plot_gsea_representative_bars.py`, `scripts/plot_gsea_top20_bars.py`, `scripts/plot_gsea_eight_panel_v2_18x10cm.py` | Male all-gene DE table, GO ontology and GOA annotations |
| `results_rna_de_female_validation/` | `scripts/validate_male_de_in_females.py` | Imported counts/TPM, sample metadata, male DE gene lists |
| `dmr_repro_cemt/output/` | `dmr_repro_cemt/scripts/dmr_limma_1kb_cemt_g2_g3.R`, `run_dmr.R`, `plot_dmr_volcano.py`, `plot_igh_tracks_corrected.py` | 13-donor 1 kb beta matrix, GENCODE v47 GTF |

`run_final_rna.slurm` and `submit_final_rna.sh` orchestrate the existing male RNA stages.
`dmr_repro_cemt/scripts/check_input.py`, `install_r_packages.R` and `submit_dmr.slurm`
provide input validation, environment installation and submission. `export_github.py` packages source files only.

## External data to restore

The code export does not include the following large files. Directory names and filenames must match.

- `$BCELLS_DATA_ROOT/aging/bcells/input_normed_binned_results/`: the 119 signal files named in `data_info/CEEHRC.CEMT*.txt`.
- `$BCELLS_DATA_ROOT/genome_features/hg38/`: 8 functional-class and 5 CpG/repeat-class `*_bins.tsv` files listed in the O-information submission scripts.
- `$BCELLS_DATA_ROOT/aging/bcells/rna_epiatlas_cll13/CEMT*.isoforms.results`: 13 files; public URLs in `data_info/ihec_rna_lookup/rna_isoforms_cll13.urls.txt`.
- `$BCELLS_DATA_ROOT/annotations/go_human/go-basic.obo` and `goa_human.gaf.gz`: restore the original snapshot, rather than silently replacing it with a different GO release. Provenance is in `data_info/external_input_provenance.json`.
- `results_rna_de_slurm/length_scaled_counts.tsv`, `gene_tpm_tximport.tsv`, `gene_effective_lengths.tsv` and `tximport.rds`: can be restored or rebuilt with `import_rna_rsem.R`. Small sample metadata, QC and gene annotation are included in the export.
- `dmr_repro_cemt/input/all_bins_1kb_cemt_g2_g3.csv` and `gencode.v47.annotation.gtf`: restore the retained beta matrix and GTF to these exact paths. The beta matrix is a prepared input; its upstream construction code was not present in this project.

To redraw without recomputing models, also restore the corresponding generated result directories.
To redraw the combined global/ABC figure, first regenerate ABC individual and combined figures.
`summary_zh.md` in C_pattern is a descriptive report: the new summary script regenerates the numerical
table and a concise computed report, rather than reconstructing the historical prose verbatim.

## Cleanup decisions

- Removed `scripts/compare_synergistic_triplets.py`: its only output, `results_triplet_R_vs_S/`, is absent.
- Removed `global_oinfo_and_triplets/plot_volcano_g2_vs_g3.py` and `plot_g2_vs_g3_order7_horizontal.py`: their standalone output directory is absent. The current combined figure read their constants through AST; these constants were extracted unchanged into `reference_config.json` before deletion.
- Kept the compact ABC preview and all three GSEA plotting variants because their outputs exist.
- Added a C-pattern recurrence generator because the retained 160-row table had no generation script.
- Copied only four necessary shared HOI modules into this project. Original shared source files elsewhere in the parent repository remain unchanged. `scripts/hoi_code/source_manifest.json` records the copied source paths and hashes.
- Removed local Python bytecode caches in the cleaned code directories. Installed environments and scientific results remain local and are excluded from the source export.

No statistical model was refitted, permutation job submitted, or existing scientific result overwritten during cleanup.
The old cleanup manifest describes an earlier directory state and remains local, excluded from the source export.
