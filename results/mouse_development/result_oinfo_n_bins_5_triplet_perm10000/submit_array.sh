#!/usr/bin/env bash
#SBATCH --job-name=dev_triplet_p10k
#SBATCH --partition=C064M0256G
#SBATCH --array=1-82
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00
#SBATCH --output=/lustre/home/2101110354/sequence_entropy/project/mouse_development/result_oinfo_n_bins_5_triplet_perm10000/logs/slurm_%A_%a.out
#SBATCH --error=/lustre/home/2101110354/sequence_entropy/project/mouse_development/result_oinfo_n_bins_5_triplet_perm10000/logs/slurm_%A_%a.err
set -euo pipefail
ROOT=/lustre/home/2101110354/sequence_entropy/project/mouse_development/result_oinfo_n_bins_5_triplet_perm10000
mapfile -t samples < "$ROOT/samples.txt"
sample="${samples[$((SLURM_ARRAY_TASK_ID - 1))]}"
exec bash "$ROOT/$sample/job.sh"
