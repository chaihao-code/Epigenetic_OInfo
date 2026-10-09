# Thyroid_Deconvolution

**Cell-composition confounding assessment** for the thyroid-aging O-information study.
Project path inside this repository: `results/aging/thyroid/deconvolution/`
(all relative paths below are rooted here).

## Purpose

Address the reviewer concern: can immune-cell infiltration / cell-composition
changes in aging thyroid explain the age-associated decline of bulk-tissue
O-information?

Scientific questions Q1–Q6 are defined in the task brief; conclusions are given
in `reports/` and are all based on actually executed analyses.

## Inputs

| Input | Description | Location / provenance |
|---|---|---|
| Methylation bigWig | 6 thyroid samples, hg38, 0–100 percent scale, **no per-CpG coverage** (CEMT0040/0042/0044/0098/0130/0132, 30–67 y) | `*CEEHRC.<sample>.WGB-Seq.methylation_profile.bigWig` (paths in `config.yaml`) |
| Histone-mark bigWig | H3K4me1/3, H3K9me3, H3K27me3/ac, H3K36me3 (Step 10 only) | same directory, `*.<mark>.signal_unstranded.bigWig` |
| Sample metadata | donor id, age, sex | `data/sample_metadata/analysis_samples.tsv` |
| O-information values | 40 genomic contexts × 6 samples (original HOI analysis) | `data/processed/oinfo_long.tsv`, prepared by `scripts/prep_oinfo_table.py` |
| Reference atlas | Loyfer et al., *Nature* 613:355–364 (2023), GEO GSE186458; 9 cell types × 3 donors, wgbstools `.hg38.beta` | downloaded by `scripts/02_download_reference.py` (gitignored, ~1.5G) |
| 1 kb bins | hg38 genomic bins of the original study | path in `config.yaml` |

## Outputs

| Output | Content |
|---|---|
| `results/tables/thyroid_cell_proportions.tsv` | Estimated cell fractions (Panel A/Model A primary; Panel B/Model B sensitivity) with chromosome-block bootstrap 95% CIs |
| `results/tables/oinfo_age_adjustment.tsv` | Per-context age slopes before/after immune-fraction adjustment, β_immune, VIF, LOO and bootstrap diagnostics (Q3/Q4) |
| `results/tables/oinfo_age_fdr_table.tsv` | Nominal p plus BH and two-stage BKY FDR for all 40 contexts (Q4) |
| `results/tables/simulation_*.tsv` | Deconvolution validation on simulated mixtures (known proportions) |
| `results/figures/` | figQ1 (composition + robustness), fig2 (age-effect adjustment panels), fig3 (volcano), fig4 (slope-ratio lollipop), plus QC figures |
| `results/qc/` | QC reports (reference identifiability, inventory, Step-9 divergence note) |
| `reports/` | `methods.md`, `methods_zh.md`, `results_summary.md`, `reviewer_response.md`, `limitations.md` |

Headline results (primary analysis): total immune fraction 1.2–10.8%,
thyroid epithelium 77.2–88.3%; after adjusting for immune fraction, the
age-associated Ω decline in CTCF-bound proximal enhancers persists
(attenuation ≤ 11%; 4 contexts nominal p < 0.05; CTCF+ proximal enhancer ×
CpG island FDR < 0.05 under both BH and BKY).

## Reproduce

**Minimal path for the Q4 conclusion** (deconvolution → regression →
significance + robustness):

```bash
bash scripts/run_q4_pipeline.sh        # all stages (3 4 5 6 7 8 9)
# or selected stages: bash scripts/run_q4_pipeline.sh 4 6 8
```

Prerequisites: reference data downloaded (`python3 scripts/02_download_reference.py`,
needs internet, ~2.4G) and `data/processed/oinfo_long.tsv` present.

Full pipeline (including Step 9/10 robustness, GO and histone analyses):

```bash
python3 scripts/01_inventory.py              # Step 1-2  data inventory
python3 scripts/02_download_reference.py     # Step 3    download GSE186458 reference
python3 scripts/03_build_reference.py        # Step 4a   reference bin-methylation matrix
python3 scripts/03b_reference_qc.py          # Step 4b   reference QC / external validation
python3 scripts/04_select_markers.py         # Step 4c   marker panels + immune-specific regions
python3 scripts/05_validate_deconvolution.py # Step 5    simulation validation
python3 scripts/06_deconvolve_samples.py     # Step 6    deconvolve 6 thyroid samples
python3 scripts/06_immune_age_association.py # Step 7    immune fraction ~ age (Q2)
python3 scripts/07_oinfo_age_adjustment.py   # Step 8    O-info age-effect adjustment (Q3/Q4)
python3 scripts/08_region_filtering.py       # Step 9    immune-region exclusion robustness
python3 scripts/09_marker_go_analysis.py     # Step 10   GO / immune-marker analysis (Q5/Q6)
python3 scripts/10_histone_association.py    # Step 10   histone-modification association
python3 scripts/17_fdr_table.py              #           nominal p + BH/BKY FDR table (Q4)
```

Figure scripts (re-runnable once tables exist): `12_fig_q1_composition.py`
(Q1 three-panel figure), `15_fig_volcano.py` (volcano), `16_fig_lollipop.py`
(slope-ratio lollipop).

Logs per step are written to `logs/` (gitignored).

## Data availability

Large files are **not** in this repository (see `.gitignore`): the GSE186458
reference betas (`data/reference/beta_hg38/`, ~1.5G), `hg38.2bit`,
`GSE186458_blocks` bed, and rebuildable `data/processed/*.npz` caches.
Reference files are re-obtained with `scripts/02_download_reference.py`;
intermediate files are rebuilt by the pipeline above.

## Key limitations (please read)

1. Thyroid methylation bigWigs are 0–100 percent scores **without coverage**;
   UXM_deconv (fragment/PAT-level) is not applicable. We used reference-based
   constrained NNLS with a beta linear-mixing approximation; this limitation is
   reported separately in the manuscripts.
2. The reference atlas (Loyfer 2023) is distributed on hg19; marker coordinates
   were converted to hg38.
3. Discretization (per-sample per-signal quintile) thresholds were not
   persisted; Step 9 recomputes them from raw signals following the original
   `signal_loader.assign_bins` logic to keep the five-bin framework consistent.

## Environment

See `environment.yml` (system Python 3.8 + external liftOver on the analysis
machine). Fixed random seeds are recorded in `config.yaml`.

## Directory layout

```
data/reference/       # downloaded GSE186458 reference (large files gitignored)
data/sample_metadata/ # sample metadata copy
data/processed/       # marker matrices, extracted methylation, O-info long table (small key intermediates)
scripts/              # 01-10 analysis + 12/15/16/17 figures & FDR + deconv_utils.py + run_q4_pipeline.sh
results/tables|figures|qc/
reports/              # methods.md, methods_zh.md, results_summary.md,
                      # reviewer_response.md, limitations.md
logs/                 # run logs (gitignored)
config.yaml           # all paths and parameters
environment.yml
```
