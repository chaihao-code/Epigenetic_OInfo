#!/bin/bash
#SBATCH --job-name=oinfo_hoi_neutrophil
#SBATCH --partition=C064M0256G
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --array=0-11
#SBATCH --output=slurm-%A_%a.out
#SBATCH --error=slurm-%A_%a.err

export EPINFO_ROOT="${EPINFO_ROOT:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && while [ ! -f .epigenetic_oinfo_root ]; do cd .. || exit 1; done && pwd)}"
set -euo pipefail

PROJECT_DIR="${EPINFO_ROOT}/cases/aging/neutrophil"
OUTPUT_DIR="${PROJECT_DIR}/results_oinfo_n_bins_5"
PYTHON_SCRIPT="${EPINFO_ROOT}/common/02_oinformation/compute_hoi_oinfo.py"
PY="python3"
GENOME_FEATURES_DIR="${EPINFO_ROOT}/data/genome_features/hg38"

mkdir -p "${OUTPUT_DIR}"

mapfile -t INPUT_FILES < <(ls "${PROJECT_DIR}"/data_info/*.txt | sort)
INPUT_FILE="${INPUT_FILES[$SLURM_ARRAY_TASK_ID]}"
SAMPLE=$(basename "${INPUT_FILE}" .txt)

# functional annotation (updated, purified, mutually exclusive) x repeat class
BLOCK_GROUP_1="${GENOME_FEATURES_DIR}/promoter_purified_bins.tsv;${GENOME_FEATURES_DIR}/proximal_ctcf_bound_enhancer_purified_bins.tsv;${GENOME_FEATURES_DIR}/proximal_not_ctcf_bound_enhancer_purified_bins.tsv;${GENOME_FEATURES_DIR}/distal_ctcf_bound_enhancer_purified_bins.tsv;${GENOME_FEATURES_DIR}/distal_not_ctcf_bound_enhancer_purified_bins.tsv;${GENOME_FEATURES_DIR}/exon_purified_bins.tsv;${GENOME_FEATURES_DIR}/intron_purified_bins.tsv"
BLOCK_GROUP_2="${GENOME_FEATURES_DIR}/CpG_island_bins.tsv;${GENOME_FEATURES_DIR}/LINE_bins.tsv;${GENOME_FEATURES_DIR}/SINE_bins.tsv;${GENOME_FEATURES_DIR}/LTR_bins.tsv"
GROUP_SEPARATOR="::BLOCKGROUP::"

echo "Job started at: $(date), host: $(hostname), array id: $SLURM_ARRAY_TASK_ID"
echo "Sample: ${SAMPLE}, input: ${INPUT_FILE}"

"${PY}" "${PYTHON_SCRIPT}" \
    --input-file "${INPUT_FILE}" \
    --output-dir "${OUTPUT_DIR}/${SAMPLE}" \
    --sample-name "${SAMPLE}" \
    --n-bins 5 \
    --n-max 100000 \
    --random-seed 42 \
    --block-groups "${BLOCK_GROUP_1}${GROUP_SEPARATOR}${BLOCK_GROUP_2}" \
    --group-separator "${GROUP_SEPARATOR}" \
    --drop-signals "" \
    --n-permutations 10000 \
    --permutation-seed 1234 \
    --permutation-scope block \
    --permutation-level system \
    --n-jobs 32

echo "Job finished at: $(date)"
