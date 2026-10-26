#!/bin/bash
#SBATCH --job-name=preimp_triplet_fdr
#SBATCH --partition=C064M0256G
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --output=slurm-%A_%a.out
#SBATCH --error=slurm-%A_%a.err
export EPINFO_ROOT="${EPINFO_ROOT:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && while [ ! -f .epigenetic_oinfo_root ]; do cd .. || exit 1; done && pwd)}"
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python3 ${EPINFO_ROOT}/cases/mouse_preimp/synergy_redundancy/synergistic_triplets/scripts/summarize_fdr.py
