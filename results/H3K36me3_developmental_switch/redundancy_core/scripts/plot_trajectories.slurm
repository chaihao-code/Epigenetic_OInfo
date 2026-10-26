#!/bin/bash
#SBATCH --job-name=core_trajectories
#SBATCH --partition=C064M0256G
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --output=/lustre/home/2101110354/sequence_entropy/project/H3K36me3_developmental_switch/redundancy_core/logs/plots_%j.out
#SBATCH --error=/lustre/home/2101110354/sequence_entropy/project/H3K36me3_developmental_switch/redundancy_core/logs/plots_%j.err
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
export MPLCONFIGDIR=/tmp/bulk16-support-mpl XDG_CACHE_HOME=/tmp/bulk16-support-cache
cd /lustre/home/2101110354/sequence_entropy/project/H3K36me3_developmental_switch/redundancy_core
/lustre/home/2101110354/anaconda3/envs/hoi_env/bin/python scripts/plot_synergistic_triplet_trajectories.py
/lustre/home/2101110354/anaconda3/envs/hoi_env/bin/python scripts/plot_redundancy_core_trajectories.py
/lustre/home/2101110354/anaconda3/envs/hoi_env/bin/python scripts/finalize_analysis.py
