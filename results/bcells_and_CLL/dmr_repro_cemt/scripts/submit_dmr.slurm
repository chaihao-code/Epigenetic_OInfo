#!/usr/bin/env bash
#SBATCH --job-name=bcells_dmr_limma
#SBATCH --partition=C064M0256G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
set -euo pipefail
PROJECT_DIR="${BCELLS_PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)}}"
PROJECT="$PROJECT_DIR/dmr_repro_cemt"
export R_LIBS_USER="$PROJECT/.Rlib"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
cd "$PROJECT"
mkdir -p output
"${BCELLS_RSCRIPT:-Rscript}" --vanilla scripts/run_dmr.R
