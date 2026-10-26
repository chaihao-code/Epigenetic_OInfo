#!/bin/bash
#SBATCH --job-name=full_signal_cores
#SBATCH --partition=C064M0256G
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --array=0-72%24
#SBATCH --output=slurm-%A_%a.out
#SBATCH --error=slurm-%A_%a.err
export EPINFO_ROOT="${EPINFO_ROOT:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && while [ ! -f .epigenetic_oinfo_root ]; do cd .. || exit 1; done && pwd)}"
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
python3 ${EPINFO_ROOT}/cases/H3K36me3_developmental_switch/redundancy_core/scripts/plot_full_signal_combination_occurrence.py --core-task "$SLURM_ARRAY_TASK_ID"
