#!/bin/bash
#SBATCH --job-name=plot_mi_n3_n10
#SBATCH --partition=C064M0256G
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --output=slurm-%A_%a.out
#SBATCH --error=slurm-%A_%a.err

export EPINFO_ROOT="${EPINFO_ROOT:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && while [ ! -f .epigenetic_oinfo_root ]; do cd .. || exit 1; done && pwd)}"
set -e
PY=python3
MI=${EPINFO_ROOT}/cases/mouse_preimp/mutual_information/scripts

for n in 3 5 10; do
    echo "=== plot_by_block_all_anchors n=${n} ==="
    ${PY} ${MI}/plot_by_block_all_anchors.py --n-bins ${n}
    echo "=== plot_mi_ratio_hist n=${n} ==="
    ${PY} ${MI}/plot_mi_ratio_hist.py --n-bins ${n}
done

echo "=== joint vs pairwise ratio heatmap (n=3,5,10) ==="
${PY} ${MI}/plot_joint_vs_pairwise_ratio_heatmap_v3.py

echo "All plotting done at $(date)"
