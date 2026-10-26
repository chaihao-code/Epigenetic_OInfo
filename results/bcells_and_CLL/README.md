# B-cell O-information, synergistic triplets, RNA-seq and methylation

Code for the analyses represented by the retained B-cell result directories.
The cohort contains 8 CLL-R, 5 CLL-S and 4 pB donors. RNA inputs cover the 13 CLL donors;
the primary differential expression model uses 4 male R and 5 male S donors.

See [CODE_MAP.md](CODE_MAP.md) for the complete code-to-result mapping and input requirements.
The GitHub source export includes code, documentation, signal manifests and small metadata inputs.
Large matrices, raw signals, generated figures/results, logs and Python/R environments remain local.
Historical result counts below describe the retained runs, not a fresh analysis performed during cleanup.

## Installation and paths

Use Python 3.11 (local recorded version: 3.11.6), R and, for array jobs, Slurm.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
export BCELLS_PROJECT_DIR="$PWD"
export BCELLS_DATA_ROOT=/absolute/path/to/data
export BCELLS_PYTHON="$PWD/.venv/bin/python"
```

Run these commands from the B-cell project root (the directory containing this README).
`BCELLS_DATA_ROOT` must contain `aging/bcells/`, `genome_features/hg38/` and `annotations/go_human/`.
In the original parent repository it defaults to that repository's `data/` directory.
`BCELLS_PROJECT_DIR` is used by batch jobs to locate this project after Slurm copies a script to its spool.
`BCELLS_PYTHON` selects the interpreter; jobs otherwise use the existing `.rna_de_env/bin/python`, if present,
then `python3`. Batch partition/time/CPU directives are examples from the original cluster; override them for your cluster.
Job logs default to Slurm's standard submission-directory filenames.

Arial is optional. Set `BCELLS_ARIAL_FONT` to an Arial TTF path to register it; most figures fall back to
DejaVu Sans when it is unavailable. Font substitution can change the appearance and layout.
The IGH track plot asks matplotlib for Arial and uses matplotlib's fallback if absent.

Recorded R environment: R 4.5.3, tximport 1.41.1, limma 3.66.0, data.table 1.18.6.1 and statmod 1.5.2.
The installer below obtains available package versions; it is not an exact historical R lockfile.

```bash
export R_LIBS_USER="$PWD/dmr_repro_cemt/.Rlib"
Rscript dmr_repro_cemt/scripts/install_r_packages.R
Rscript -e 'BiocManager::install("tximport", ask=FALSE, update=FALSE)'
```

The pinned Python versions were read from the retained environment. Dependency installation in a clean
machine has not been validated during this source cleanup.

## O-information and triplets

Seven tracks: H3K4me1, H3K4me3, H3K9me3, H3K27me3, H3K36me3, H3K27ac and DNA methylation.
Each sample uses 40 blocks (8 functional classes × 5 CpG/repeat classes), five discrete states,
at most 100,000 loci per block (seed 42), and 1,000 independent-track permutations (seed 1234).
The four required shared HOI modules are included under `scripts/hoi_code/`; their source hashes are recorded there.
The manifest resolver maps the retained signal lists onto `BCELLS_DATA_ROOT` without changing the inputs.

```bash
sbatch scripts/submit_oinfo_hoi.slurm
triplet_job_id=$(sbatch --parsable scripts/submit_oinfo_triplets.slurm)
sbatch --dependency=afterok:"$triplet_job_id" scripts/aggregate_oinfo_triplets.slurm
```

To summarize existing results and regenerate figures (these commands write into result directories):

```bash
python scripts/summarize_triplet_significance.py
python scripts/plot_abc_triplet_classes.py
python scripts/plot_abc_pooled_combined.py
python scripts/plot_abc_pooled_compact_preview.py
python scripts/summarize_c_triplet_recurrence.py
python global_oinfo_and_triplets/plot_global_oinfo_and_ABC_triplets.py
```

There are 35 triplets per block (A/B/C capacities: 12/19/4). Sample-level BH correction is applied
separately to each tail over 1,400 tests; global BH uses 23,800 tests across all 17 donors.
Original individual heatmaps use negative Ω and nominal lower-tail P<0.05; `_FDR005` figures and the
combined ABC figure use negative Ω and global lower-tail BH q<0.05. The global order-7 panel compares
R versus S with the retained directional Welch test P(R>S), with nominal P<0.05 coloring.
`reference_config.json` retains the group and axis constants extracted from the obsolete reference scripts.

## RNA differential expression and GO GSEA

Restore the large imported matrices to `results_rna_de_slurm/` (see CODE_MAP), or rebuild them
from the RSEM isoform files. Sample metadata, gene annotations and QC inputs are included in the source export.
The importer rebuilds expression/length matrices; it does not regenerate these curated metadata inputs.

```bash
Rscript scripts/import_rna_rsem.R "$PWD"
# Run locally, or submit with the helper below:
bash scripts/run_final_rna.slurm all
# bash scripts/submit_final_rna.sh all
python scripts/validate_male_de_in_females.py
```

Stages `de`, `gsea` and `figures` can replace `all`. The male DE model is
`RNA method + IGHV + R_vs_S`; positive log2FoldChange means higher in R.
The retained run has 38 genes at FDR<0.05 and 60 at FDR<0.10.
GSEA ranks all eligible genes by Wald statistic and uses the retained GO ontology/GOA snapshot.
There are four female R donors and no female S donors: female checks assess signature concordance,
not an independent female R/S differential expression test.
RNA URL lists are under `data_info/ihec_rna_lookup/`. The download helper downloads gene-level files
and QC only; tximport requires the separately listed **isoform** files.

## 1 kb differential methylation

Restore the two large files in `dmr_repro_cemt/input/` described in CODE_MAP, then run:

```bash
python dmr_repro_cemt/scripts/check_input.py
export R_LIBS_USER="$PWD/dmr_repro_cemt/.Rlib"
Rscript dmr_repro_cemt/scripts/run_dmr.R
python dmr_repro_cemt/scripts/plot_dmr_volcano.py
python dmr_repro_cemt/scripts/plot_igh_tracks_corrected.py
# Alternatively submit the model with:
# sbatch dmr_repro_cemt/scripts/submit_dmr.slurm
```

The supplied limma analysis is retained unchanged: M-values adjusted for age/sex, autosomes only,
minimum 8 nonmissing donors, BH q<0.05 and |delta_beta|>0.2 after its original rounding.
The retained run has 15 passing windows (R-high 11, S-high 4). Adjacent windows are not merged.
See [dmr_repro_cemt/README.md](dmr_repro_cemt/README.md) for model details and limitations.

## GitHub source export

```bash
python scripts/export_github.py /path/to/new/bcells_github
```

The destination must not exist. Export uses an explicit allowlist, writes `SOURCE_MANIFEST.json`
with per-file SHA256, and omits generated outputs and installed environments.
`.gitignore` also excludes these files if working directly in the analysis directory.
There is no GitHub upload or remote configuration in this workflow.
