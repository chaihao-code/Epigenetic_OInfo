#!/bin/bash
#SBATCH --job-name=scan_block
#SBATCH --partition=C064M0256G
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --output=slurm-%A_%a.out
#SBATCH --error=slurm-%A_%a.err
export EPINFO_ROOT="${EPINFO_ROOT:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && while [ ! -f .epigenetic_oinfo_root ]; do cd .. || exit 1; done && pwd)}"
python3 ${EPINFO_ROOT}/cases/mouse_preimp/mutual_information/scripts/scan_block_candidates.py
