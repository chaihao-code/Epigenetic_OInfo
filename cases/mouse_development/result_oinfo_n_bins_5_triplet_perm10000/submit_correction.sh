#!/usr/bin/env bash
#SBATCH --job-name=dev_triplet10k_fdr
#SBATCH --partition=C064M0256G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=slurm-%A_%a.out
#SBATCH --error=slurm-%A_%a.err
export EPINFO_ROOT="${EPINFO_ROOT:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && while [ ! -f .epigenetic_oinfo_root ]; do cd .. || exit 1; done && pwd)}"
set -euo pipefail
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
python3 -u ${EPINFO_ROOT}/cases/mouse_development/result_oinfo_n_bins_5_triplet_perm10000/code/finish_analysis.py
