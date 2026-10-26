#!/bin/bash
#SBATCH --job-name=preimp_triplet_fdr
#SBATCH --partition=C064M0256G
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --output=/lustre/home/2101110354/sequence_entropy/project/mouse_preimp/synergy_redundancy/synergistic_triplets/logs/fdr_%j.out
#SBATCH --error=/lustre/home/2101110354/sequence_entropy/project/mouse_preimp/synergy_redundancy/synergistic_triplets/logs/fdr_%j.err
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
/lustre/home/2101110354/anaconda3/envs/hoi_env/bin/python /lustre/home/2101110354/sequence_entropy/project/mouse_preimp/synergy_redundancy/synergistic_triplets/scripts/summarize_fdr.py
