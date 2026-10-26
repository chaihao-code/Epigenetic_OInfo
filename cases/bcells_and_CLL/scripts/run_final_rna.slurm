#!/usr/bin/env bash
#SBATCH --job-name=bcells_final_rna
#SBATCH --partition=C064M0256G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
set -euo pipefail
PROJECT_DIR="${BCELLS_PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)}}"
export BCELLS_DATA_ROOT="${BCELLS_DATA_ROOT:-$(cd "$PROJECT_DIR/../.." && pwd)/data}"
if [[ -z "${BCELLS_PYTHON:-}" && -x "$PROJECT_DIR/.rna_de_env/bin/python" ]]; then
    BCELLS_PYTHON="$PROJECT_DIR/.rna_de_env/bin/python"
fi
PY="${BCELLS_PYTHON:-python3}"
PROJECT="$PROJECT_DIR"
export MPLCONFIGDIR="$PROJECT/.cache/matplotlib"
export XDG_CACHE_HOME="$PROJECT/.cache"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
mkdir -p "$MPLCONFIGDIR"
stage=${1:-all}
case "$stage" in de|gsea|figures|all) ;; *) echo "Usage: $0 {de|gsea|figures|all}" >&2; exit 2;; esac
if [[ "$stage" == de || "$stage" == all ]]; then
 "$PY" "$PROJECT/scripts/analyze_rna_de_male.py"
fi
if [[ "$stage" == gsea || "$stage" == all ]]; then
 "$PY" "$PROJECT/scripts/gsea_rna_male.py"
fi
if [[ "$stage" == figures || "$stage" == all ]]; then
 for script in plot_de_male_ighv_volcano_clear.py plot_gsea_representative_bars.py plot_gsea_top20_bars.py plot_gsea_eight_panel_v2_18x10cm.py; do
  "$PY" "$PROJECT/scripts/$script"
 done
fi
