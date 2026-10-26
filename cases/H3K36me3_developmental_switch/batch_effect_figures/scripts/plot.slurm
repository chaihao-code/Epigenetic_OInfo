#!/bin/bash
#SBATCH --job-name=bulk16_FDR_stability
#SBATCH --partition=C064M0256G
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm-%A_%a.out
#SBATCH --error=slurm-%A_%a.err
export EPINFO_ROOT="${EPINFO_ROOT:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && while [ ! -f .epigenetic_oinfo_root ]; do cd .. || exit 1; done && pwd)}"
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export MPLCONFIGDIR=/tmp/bulk16-fdr-mpl XDG_CACHE_HOME=/tmp/bulk16-fdr-cache
python3 ${EPINFO_ROOT}/cases/H3K36me3_developmental_switch/batch_effect_figures/scripts/plot_batch_effect.py
