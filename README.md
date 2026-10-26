# Epigenetic O-information

This release organizes the biological cases already present in the upstream repository. EMT, motif/TF-binding analyses and other local-only projects are excluded. The original genome windows and annotation memberships are fixed inputs and remain unchanged.

Supplementary exports, including the three English result workbooks and region-bin count tables, are kept locally under `supplementary_data/` and excluded from this GitHub snapshot. Browse [figures by case](cases/FIGURES.md). Only `figQ1_composition.png` is retained for thyroid deconvolution; thyroid GSEA and the other cases' figures are included.

## Start with the complete thyroid demo

Read [demo/README.md](demo/README.md). It uses one complete healthy human thyroid donor, CEMT0098 (30-year-old female), all seven signals, all original hg38 1-kb windows and all 40 genomic annotation classes. It can begin either with the supplied full prepared `.bp` tracks or the full original bigWigs. The six-donor thyroid case configuration is also retained. There is no cropped-window demo.

```bash
cd /path/to/Epigenetic_OInfo_release
export EPINFO_ROOT="$PWD"
python -m pip install -r demo/requirements.txt
python demo/run_demo.py --jobs 32
# Full bigWig route: restore the eight full files in demo/input/bigwig first.
python demo/run_demo.py --from-bigwig --jobs 32
```

The bigWig route also needs the UCSC `bigWigAverageOverBed` executable on PATH. The whole-genome demo is intended for a compute node. `demo/submit_demo.slurm` runs it with 32 CPU cores; set `EPINFO_PYTHON` and adjust cluster partition/time settings before submission.

## One supported O-information implementation

`common/02_oinformation/compute_hoi_oinfo.py` is the only numerical engine for global and local O-information. Former case engine paths are compatibility forwarders. The old entropy-only and entropy-to-O-information programs and their old launchers have been removed. New CLI runs require a positive permutation count.

The engine defaults to **10,000 permutations** for synergistic triplets, with independent track shuffling inside each existing class, a +1 empirical-P correction and lower-tail BH correction. BH includes every triplet in the chosen family before selecting negative Omega and q <= 0.05. The default family is one sample × class; B-cell per-sample and multi-sample global families are explicit alternatives. The thyroid demo defaults to its full-system permutation null; `--triplets` is an optional generic triplet/FDR demonstration.

## Layout

| Directory | Purpose |
|---|---|
| `common/01_preprocessing/` | Fixed-grid signal binning, Input normalization, discretization and annotation utilities |
| `common/02_oinformation/` | Single global/local O-information and permutation engine |
| `common/03_triplets_and_redundancy/` | Triplet BH correction and the separate replacement-background redundancy-core method |
| `common/04_local_oinformation/` | Thyroid local-regression adapters using the shared engine |
| `common/05_dmr/` | B-cell differential methylation reference implementation |
| `common/06_de_and_gsea/` | B-cell differential expression and case-specific GSEA reference implementations |
| `cases/` | The seven biological case directories already represented upstream |
| `demo/` | Complete thyroid donor, original grid, preprocessing and executable validation |
| `supplementary_data/` (local only) | Supplementary exports and final workbooks; excluded from GitHub |

Directory READMEs describe the workflows and inputs. Machine-specific external-input symlinks are omitted from the GitHub checkout; their locations are recorded in `config/local_input_links.tsv`. Complete prepared thyroid tracks and the unchanged grid are supplied as compressed files. Full bigWig source locations and their original download/conversion formats are documented in `demo/input/bigwig_sources.tsv`; methylation's public source is a CpG BED.gz, not a directly downloadable bigWig.

Numerical result files and supplied compressed inputs are included in Git. The complete 2,697,259-window B-cell DMR table is supplied as `all_bins_1kb_cemt_g2_g3_limma_results.csv.gz`; its plotting entrypoint reads the compressed file automatically when the uncompressed CSV is absent. Full-bigWig preprocessing outputs can be regenerated with the demo and are not duplicated in Git.

## Validation and historical results

The [full figure rerun audit](validation/plot_rerun_20261026/REPORT.md) lists each plotting entrypoint, newly generated figures, input dependencies and missing historical inputs. Its [TSV inventory](validation/plot_rerun_20261026/plot_inventory.tsv) and [file manifest](validation/plot_rerun_20261026/figure_files.tsv) support review. Reproduce the figure workflows with `python validation/plot_rerun_20261026/run_all.py`; mixed analysis/plotting programs use retained intermediate tables through the supplied plotting runners.

Run `python -m unittest discover -s tests -v`. Scientific checks cover known synergy/redundancy values, local/global equality, non-five-state encoding, BH family boundaries and the 10,000-permutation defaults. Actual thyroid execution evidence is written to `demo/results/<mode>/validation.json` when available. Older cohort result tables are not relabelled as newly recomputed 10,000-permutation analyses.

See [docs/PATHS.md](docs/PATHS.md), [docs/CODE_MAP.md](docs/CODE_MAP.md) and [the Chinese change notes](docs/整理说明_中文.md).
