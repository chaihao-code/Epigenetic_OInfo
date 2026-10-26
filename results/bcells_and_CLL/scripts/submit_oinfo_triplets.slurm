#!/bin/bash
#SBATCH --job-name=oinfo_triplet_bcells
#SBATCH --partition=C064M0256G
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --array=0-16

set -euo pipefail

PROJECT_DIR="${BCELLS_PROJECT_DIR:-${SLURM_SUBMIT_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)}}"
export BCELLS_DATA_ROOT="${BCELLS_DATA_ROOT:-$(cd "$PROJECT_DIR/../../.." && pwd)/data}"
if [[ -z "${BCELLS_PYTHON:-}" && -x "$PROJECT_DIR/.rna_de_env/bin/python" ]]; then
    BCELLS_PYTHON="$PROJECT_DIR/.rna_de_env/bin/python"
fi
PY="${BCELLS_PYTHON:-python3}"
OUTPUT_DIR="${PROJECT_DIR}/results_oinfo_triplet_n_bins_5"
PYTHON_SCRIPT="$PROJECT_DIR/scripts/hoi_code/compute_hoi_oinfo.py"
GENOME_FEATURES_DIR="${BCELLS_GENOME_FEATURES_DIR:-$BCELLS_DATA_ROOT/genome_features/hg38}"

mkdir -p "${OUTPUT_DIR}"

mapfile -t INPUT_FILES < <(printf '%s\n' "${PROJECT_DIR}"/data_info/CEEHRC.CEMT*.txt | sort)
[[ ${#INPUT_FILES[@]} -eq 17 ]] || { echo "Expected 17 sample manifests" >&2; exit 1; }
INPUT_FILE="${INPUT_FILES[$SLURM_ARRAY_TASK_ID]}"
SAMPLE=$(basename "${INPUT_FILE}" .txt)
RESOLVED_INPUT=$(mktemp)
trap 'rm -f "$RESOLVED_INPUT"' EXIT
"$PY" "$PROJECT_DIR/scripts/resolve_signal_list.py" "$INPUT_FILE" "$BCELLS_DATA_ROOT" > "$RESOLVED_INPUT"
INPUT_FILE="$RESOLVED_INPUT"

# functional annotation (purified, mutually exclusive) x repeat class, 8 x 5 = 40 blocks
BLOCK_GROUP_1="${GENOME_FEATURES_DIR}/promoter_purified_bins.tsv;${GENOME_FEATURES_DIR}/proximal_ctcf_bound_enhancer_purified_bins.tsv;${GENOME_FEATURES_DIR}/proximal_not_ctcf_bound_enhancer_purified_bins.tsv;${GENOME_FEATURES_DIR}/distal_ctcf_bound_enhancer_purified_bins.tsv;${GENOME_FEATURES_DIR}/distal_not_ctcf_bound_enhancer_purified_bins.tsv;${GENOME_FEATURES_DIR}/exon_purified_bins.tsv;${GENOME_FEATURES_DIR}/intron_purified_bins.tsv;${GENOME_FEATURES_DIR}/None_func_bins.tsv"
BLOCK_GROUP_2="${GENOME_FEATURES_DIR}/CpG_island_bins.tsv;${GENOME_FEATURES_DIR}/LINE_bins.tsv;${GENOME_FEATURES_DIR}/SINE_bins.tsv;${GENOME_FEATURES_DIR}/LTR_bins.tsv;${GENOME_FEATURES_DIR}/None_rep_bins.tsv"
GROUP_SEPARATOR="::BLOCKGROUP::"

echo "Job started at: $(date), host: $(hostname), array id: $SLURM_ARRAY_TASK_ID"
echo "Sample: ${SAMPLE}, input: ${INPUT_FILE}"

echo "=== genome blocks (8 func x 5 repeat = 40) ==="
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
    --n-permutations 1000 \
    --permutation-seed 1234 \
    --permutation-scope block \
    --permutation-level all \
    --permutation-max-order 3 \
    --n-jobs "${SLURM_CPUS_PER_TASK:-32}"


"${PY}" "${PROJECT_DIR}/scripts/summarize_triplet_significance.py" --sample "${SAMPLE}"
echo "Job finished at: $(date)"
