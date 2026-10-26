#!/usr/bin/env bash
#SBATCH --job-name=thyroid_bins_demo
#SBATCH --partition=C064M0256G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=01:00:00
#SBATCH --output=thyroid_demo_%j.out
#SBATCH --error=thyroid_demo_%j.err
set -euo pipefail
: "${EPINFO_ROOT:?Export the release root before submitting this job}"
EPINFO_PYTHON="${EPINFO_PYTHON:-python3}"
export PATH="$(dirname "$EPINFO_PYTHON"):$PATH"
export PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$EPINFO_ROOT"
"$EPINFO_PYTHON" -m unittest discover -s tests -v
"$EPINFO_PYTHON" demo/run_demo.py --n-permutations 10000 --jobs "${SLURM_CPUS_PER_TASK:-32}"
