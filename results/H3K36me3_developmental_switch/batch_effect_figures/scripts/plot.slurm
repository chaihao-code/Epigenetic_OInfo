#!/bin/bash
#SBATCH --job-name=bulk16_FDR_stability
#SBATCH --partition=C064M0256G
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=/lustre/home/2101110354/sequence_entropy/project/H3K36me3_developmental_switch/batch_effect_figures/logs/report_%j.out
#SBATCH --error=/lustre/home/2101110354/sequence_entropy/project/H3K36me3_developmental_switch/batch_effect_figures/logs/report_%j.err
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export MPLCONFIGDIR=/tmp/bulk16-fdr-mpl XDG_CACHE_HOME=/tmp/bulk16-fdr-cache
/lustre/home/2101110354/anaconda3/envs/hoi_env/bin/python /lustre/home/2101110354/sequence_entropy/project/H3K36me3_developmental_switch/batch_effect_figures/scripts/plot_batch_effect.py
