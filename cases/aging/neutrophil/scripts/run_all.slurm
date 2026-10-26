#!/bin/bash
#SBATCH --job-name=neutrophil_pp
#SBATCH --partition=C064M0256G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=04:00:00
#SBATCH --output=slurm-%A_%a.out
#SBATCH --error=slurm-%A_%a.err

export EPINFO_ROOT="${EPINFO_ROOT:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && while [ ! -f .epigenetic_oinfo_root ]; do cd .. || exit 1; done && pwd)}"
set -uo pipefail

export PATH="${EPINFO_TOOLS_BIN:-$EPINFO_ROOT/tools/bin}:/lustre/home/2101110354/packages/bedtools:${PATH}"

SCRIPT_DIR="${EPINFO_ROOT}/cases/aging/neutrophil/scripts"

echo "Job started at $(date)"
echo "PATH check: $(which bedtools bigWigToBedGraph bigWigAverageOverBed bedGraphToBigWig | tr '\n' ' ')"

for step in bin_histone weighted_methylation divinput discretize make_datainfo; do
    echo ""
    echo "==================== Step: ${step} ($(date)) ===================="
    bash "${SCRIPT_DIR}/${step}.sh" || { echo "STEP FAILED: ${step}"; exit 1; }
done

echo ""
echo "All steps finished at $(date)"
