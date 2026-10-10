"""Shared configuration for the immune-adjusted O-information pipeline.

All paths and analysis parameters live here so that scripts stay generic.
Edit this file to adapt the pipeline to another dataset.
"""

from pathlib import Path

# ---------------------------------------------------------------------------
# Project root (this file lives in <root>/scripts/)
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent

# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------
# Directory holding one subdirectory per genomic region, each produced by
# run_analysis_generic_excluded.py and containing:
#   oinfo_per_sample.bed   (per-bin local O-information for 6 samples + age regression)
#   extracted_signals.tsv  (per-bin discrete states 0-4 for 7 marks x 6 samples)
REGION_BASE = PROJECT_ROOT.parent / "local_oinfo_age_correlation"

# Cell-proportion table from the thyroid deconvolution project (long format:
# sample_id, panel, model, cell_type, fraction, ...). Panel A / Model A is the
# primary (independently screened 2500-marker, 9 fine-grained cell types) model.
PROPORTIONS_TSV = (
    PROJECT_ROOT.parent
    / "deconvolution/results/aging/thyroid/deconvolution/results/tables/thyroid_cell_proportions.tsv"
)
PANEL = "A"
MODEL = "A"

# Cell types that make up the immune compartment (must match names in
# PROPORTIONS_TSV for the chosen panel/model). Total immune fraction is their
# sum per sample.
IMMUNE_CELL_TYPES = ["T_CD4", "T_CD8", "B_cell", "NK", "Monocyte", "Granulocyte"]

# Gene annotation for TSS mapping (Gencode GTF, hg38) and window half-width.
GTF_PATH = Path(
    "/media/dell/data5/sequence_feature/sequence_entropy/data/genome_annotation"
    "/hg38/gencode.v47.annotation.gtf"
)
TSS_WINDOW = 2000  # bins within +/-2 kb of a TSS are assigned to that gene

# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------
OUT_REGRESSION = PROJECT_ROOT / "results/regression"
OUT_GSEA_PATHWAY = PROJECT_ROOT / "results/gsea_pathway"
OUT_GSEA_MARGINAL = PROJECT_ROOT / "results/gsea_marginal"
OUT_TABLES = PROJECT_ROOT / "results/tables"
LOG_DIR = PROJECT_ROOT / "logs"
for _d in (OUT_REGRESSION, OUT_GSEA_PATHWAY, OUT_GSEA_MARGINAL, OUT_TABLES, LOG_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Regression / GSEA parameters (mirroring the original unadjusted pipeline)
# ---------------------------------------------------------------------------
MARKS = ["H3K4me1", "H3K4me3", "H3K9me3", "H3K27me3", "H3K36me3", "H3K27ac", "WGB-Seq"]
SAMPLES = ["CEMT0098", "CEMT0130", "CEMT0132", "CEMT0042", "CEMT0044", "CEMT0040"]
AGES = [30, 37, 37, 46, 55, 67]

# Significance rule for the marginal (signal-change) GSEA, as in the original
# bin-level marginal analysis.
TREND_R_THRESHOLD = 0.3
TREND_P_THRESHOLD = 0.05
MIN_SET_SIZE = 5

PATHWAY_LIBRARIES = ["GO_Biological_Process_2021", "MSigDB_Hallmark_2020"]
PATHWAY_PERMUTATIONS = 1000
PATHWAY_MAX_SIZE = 1000
MARGINAL_PERMUTATIONS = 10000
MARGINAL_MAX_SIZE = 100000
GSEA_SEED = 42
# A GSEA term is called significant only if BOTH thresholds are met (the dual
# rule excludes degenerate flat-class artifacts, see methods text).
GSEA_P_CUTOFF = 0.05
GSEA_FDR_CUTOFF = 0.25
