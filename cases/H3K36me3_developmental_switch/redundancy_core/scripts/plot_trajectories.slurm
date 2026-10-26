#!/bin/bash
#SBATCH --job-name=core_trajectories
#SBATCH --partition=C064M0256G
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm-%A_%a.out
#SBATCH --error=slurm-%A_%a.err
export EPINFO_ROOT="${EPINFO_ROOT:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && while [ ! -f .epigenetic_oinfo_root ]; do cd .. || exit 1; done && pwd)}"
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
export MPLCONFIGDIR=/tmp/bulk16-support-mpl XDG_CACHE_HOME=/tmp/bulk16-support-cache
cd ${EPINFO_ROOT}/cases/H3K36me3_developmental_switch/redundancy_core
python3 scripts/plot_synergistic_triplet_trajectories.py
python3 scripts/plot_redundancy_core_trajectories.py
python3 scripts/finalize_analysis.py
