#!/bin/bash
#SBATCH --job-name=neutrophil_pp
#SBATCH --partition=C064M0256G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=04:00:00
#SBATCH --output=/lustre/home/2101110354/sequence_entropy/jobs/neutrophil_preprocess_%j.out
#SBATCH --error=/lustre/home/2101110354/sequence_entropy/jobs/neutrophil_preprocess_%j.err

set -uo pipefail

export PATH="/lustre/home/2101110354/anaconda3/envs/hoi_env/bin:/lustre/home/2101110354/packages/bedtools:${PATH}"

SCRIPT_DIR="/lustre/home/2101110354/sequence_entropy/project/aging/neutrophil/scripts"

echo "Job started at $(date)"
echo "PATH check: $(which bedtools bigWigToBedGraph bigWigAverageOverBed bedGraphToBigWig | tr '\n' ' ')"

for step in bin_histone weighted_methylation divinput discretize make_datainfo; do
    echo ""
    echo "==================== Step: ${step} ($(date)) ===================="
    bash "${SCRIPT_DIR}/${step}.sh" || { echo "STEP FAILED: ${step}"; exit 1; }
done

echo ""
echo "All steps finished at $(date)"
