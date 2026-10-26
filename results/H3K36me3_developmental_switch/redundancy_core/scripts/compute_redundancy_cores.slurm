#!/bin/bash
#SBATCH --job-name=development_cores
#SBATCH --partition=C064M0256G
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --array=0-72
#SBATCH --output=/lustre/home/2101110354/sequence_entropy/project/H3K36me3_developmental_switch/redundancy_core/logs/core_%A_%a.out
#SBATCH --error=/lustre/home/2101110354/sequence_entropy/project/H3K36me3_developmental_switch/redundancy_core/logs/core_%A_%a.err
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
/lustre/home/2101110354/anaconda3/envs/hoi_env/bin/python /lustre/home/2101110354/sequence_entropy/project/H3K36me3_developmental_switch/redundancy_core/scripts/run_core_task.py
