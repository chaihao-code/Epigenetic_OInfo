#!/usr/bin/env bash
#SBATCH --job-name=dev_triplet10k_fdr
#SBATCH --partition=C064M0256G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=/lustre/home/2101110354/sequence_entropy/project/mouse_development/result_oinfo_n_bins_5_triplet_perm10000/logs/correction_%j.out
#SBATCH --error=/lustre/home/2101110354/sequence_entropy/project/mouse_development/result_oinfo_n_bins_5_triplet_perm10000/logs/correction_%j.err
set -euo pipefail
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
/lustre/home/2101110354/anaconda3/envs/hoi_env/bin/python3 -u /lustre/home/2101110354/sequence_entropy/project/mouse_development/result_oinfo_n_bins_5_triplet_perm10000/code/finish_analysis.py
