#!/usr/bin/env bash
#SBATCH --job-name=dev_triplet_p10k
#SBATCH --partition=C064M0256G
#SBATCH --array=1-82
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%A_%a.out
#SBATCH --error=slurm-%A_%a.err
export EPINFO_ROOT="${EPINFO_ROOT:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && while [ ! -f .epigenetic_oinfo_root ]; do cd .. || exit 1; done && pwd)}"
set -euo pipefail
ROOT=${EPINFO_ROOT}/cases/mouse_development/result_oinfo_n_bins_5_triplet_perm10000
mapfile -t samples < "$ROOT/samples.txt"
sample="${samples[$((SLURM_ARRAY_TASK_ID - 1))]}"
exec bash "$ROOT/$sample/job.sh"
