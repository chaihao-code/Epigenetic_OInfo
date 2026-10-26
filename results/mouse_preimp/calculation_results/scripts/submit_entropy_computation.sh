#!/bin/bash
#
# Submit All-Subset Entropy analysis jobs to SLURM cluster
# Computes joint entropy H(S) for all non-empty subsets S of features
#
# Usage:
#   bash submit_entropy_jobs.sh [OPTIONS]
#

set -e

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
INPUT_DIR="${PROJECT_DIR}/data_info"
RESULTS_DIR="${PROJECT_DIR}/results_entropy_n_bins_5"
PYTHON_SCRIPT="/lustre/home/2101110354/sequence_entropy/script/hoi_code/compute_hoi_entropy.py"

# Genome feature directory (mm10 for mouse embryo development)
GENOME_FEATURES_DIR="/lustre/home/2101110354/sequence_entropy/data/genome_features/mm10"

# Default block groups
BLOCK_GROUP_1="${GENOME_FEATURES_DIR}/promoter_bins.tsv;${GENOME_FEATURES_DIR}/proximal_ctcf_bound_enhancer_purified_bins.tsv;${GENOME_FEATURES_DIR}/proximal_not_ctcf_bound_enhancer_purified_bins.tsv;${GENOME_FEATURES_DIR}/distal_ctcf_bound_enhancer_purified_bins.tsv;${GENOME_FEATURES_DIR}/distal_not_ctcf_bound_enhancer_purified_bins.tsv;${GENOME_FEATURES_DIR}/exon_purified_bins.tsv;${GENOME_FEATURES_DIR}/intron_purified_bins.tsv"
BLOCK_GROUP_2="${GENOME_FEATURES_DIR}/CpG_island_bins.tsv;${GENOME_FEATURES_DIR}/LINE_bins.tsv;${GENOME_FEATURES_DIR}/SINE_bins.tsv;${GENOME_FEATURES_DIR}/LTR_bins.tsv"

# Group separator for passing to Python
GROUP_SEPARATOR="::BLOCKGROUP::"

# SLURM configuration
PARTITION="C064M1024G"
TIME="02:00:00"
NODES=1
NTASKS=1
CPUS_PER_TASK=32
MEM="32G"

# HOI parameters
N_MAX="--n-max 100000"
METHOD="binning"

# JAX memory control
DISABLE_JIT=1

# Block groups
custom_block_group_1=""
custom_block_group_2=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --partition)
            PARTITION="$2"
            shift 2
            ;;
        --time)
            TIME="$2"
            shift 2
            ;;
        --mem)
            MEM="$2"
            shift 2
            ;;
        --cpus)
            CPUS_PER_TASK="$2"
            shift 2
            ;;
        --n-max)
            N_MAX="--n-max $2"
            shift 2
            ;;
        --use-all)
            N_MAX=""
            shift
            ;;
        --method)
            METHOD="$2"
            shift 2
            ;;
        --block-group-1)
            shift
            custom_block_group_1=""
            while [[ $# -gt 0 ]] && [[ ! "$1" =~ ^-- ]]; do
                if [ -z "$custom_block_group_1" ]; then
                    custom_block_group_1="$1"
                else
                    custom_block_group_1="${custom_block_group_1};$1"
                fi
                shift
            done
            ;;
        --block-group-2)
            shift
            custom_block_group_2=""
            while [[ $# -gt 0 ]] && [[ ! "$1" =~ ^-- ]]; do
                if [ -z "$custom_block_group_2" ]; then
                    custom_block_group_2="$1"
                else
                    custom_block_group_2="${custom_block_group_2};$1"
                fi
                shift
            done
            ;;
        --disable-jit)
            DISABLE_JIT=1
            shift
            ;;
        --help|-h)
            echo "Usage: bash submit_entropy_jobs.sh [OPTIONS]"
            echo ""
            echo "Computes joint entropy H(S) for all non-empty subsets of features"
            echo ""
            echo "SLURM Options:"
            echo "  --partition NAME    SLURM partition (default: cpu1,cpu2)"
            echo "  --time TIME         Time limit (default: 36:00:00)"
            echo "  --mem MEM           Memory per job (default: 32G)"
            echo "  --cpus N            CPUs per task (default: 20, also sets --n-jobs)"
            echo ""
            echo "Analysis Options:"
            echo "  --n-max N           Maximum samples per file (default: 100000)"
            echo "  --use-all           Use all samples (no limit)"
            echo "  --method METHOD     Entropy method (default: binning)"
            echo ""
            echo "Block Group Options:"
            echo "  --block-group-1 F...  Group 1 files (mutually exclusive)"
            echo "  --block-group-2 F...  Group 2 files (Cartesian product)"
            echo ""
            echo "Memory Options:"
            echo "  --disable-jit       Disable JAX JIT (slower but less memory)"
            echo ""
            echo "Example:"
            echo "  bash submit_entropy_jobs.sh --use-all"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Use custom block groups if provided
if [ -n "$custom_block_group_1" ]; then
    BLOCK_GROUP_1="$custom_block_group_1"
fi
if [ -n "$custom_block_group_2" ]; then
    BLOCK_GROUP_2="$custom_block_group_2"
fi

# Build block groups string
if [ -n "$BLOCK_GROUP_2" ]; then
    BLOCK_GROUPS="${BLOCK_GROUP_1}${GROUP_SEPARATOR}${BLOCK_GROUP_2}"
else
    BLOCK_GROUPS="${BLOCK_GROUP_1}"
fi

# Check if sbatch is available
if ! command -v sbatch &> /dev/null; then
    echo "Error: sbatch command not found. Are you on a SLURM cluster?"
    exit 1
fi

# Get list of input files
INPUT_FILES=($(ls -1 "${INPUT_DIR}"/*.txt 2>/dev/null))
N_FILES=${#INPUT_FILES[@]}

if [ $N_FILES -eq 0 ]; then
    echo "Error: No input files found in ${INPUT_DIR}"
    exit 1
fi

echo "======================================================================"
echo "All-Subset Entropy Job Submission - Mouse Embryo Development"
echo "======================================================================"
echo "Number of jobs to submit: ${N_FILES}"
echo "SLURM partition: ${PARTITION}"
echo "Time limit: ${TIME}"
echo "Memory: ${MEM}"
echo "CPUs per task: ${CPUS_PER_TASK}"
echo ""
echo "Entropy parameters:"
echo "  method=${METHOD}"
[ -n "$N_MAX" ] && echo "  n_max=${N_MAX#--n-max }" || echo "  n_max=All (unlimited)"
echo ""
echo "Output format:"
echo "  - *_entropy_all_subsets.csv : H(S) for all non-empty subsets S"
echo ""
echo "Block Groups:"
echo "  Group 1 (mutually exclusive):"
echo "    $(echo "$BLOCK_GROUP_1" | tr ';' '\n' | sed 's/^/      - /')"
if [ -n "$BLOCK_GROUP_2" ]; then
    echo "  Group 2:"
    echo "    $(echo "$BLOCK_GROUP_2" | tr ';' '\n' | sed 's/^/      - /')"
fi
echo "======================================================================"
echo ""
read -p "Press Enter to continue or Ctrl+C to cancel..."

# Submit jobs
for input_file in "${INPUT_FILES[@]}"; do
    sample_name=$(basename "$input_file" .txt)
    output_dir="${RESULTS_DIR}/${sample_name}"
    job_name="entropy_${sample_name}"
    
    mkdir -p "$output_dir"
    
    cat > "${output_dir}/job.sh" << EEOF
#!/bin/bash
#SBATCH --job-name=${job_name}
#SBATCH --partition=${PARTITION}
#SBATCH --time=${TIME}
#SBATCH --nodes=${NODES}
#SBATCH --ntasks=${NTASKS}
#SBATCH --cpus-per-task=${CPUS_PER_TASK}
#SBATCH --output=${output_dir}/slurm_%j.out
#SBATCH --error=${output_dir}/slurm_%j.err

echo "Job started at: \$(date)"
echo "Running on host: \$(hostname)"
echo "SLURM job ID: \$SLURM_JOB_ID"

# Set JAX memory limits
export JAX_PLATFORM_NAME=cpu
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_PYTHON_CLIENT_ALLOCATOR=platform
export XLA_FLAGS="--xla_cpu_use_thunk_runtime=false"

if [ "${DISABLE_JIT}" == "1" ]; then
    echo "WARNING: JAX JIT disabled for memory saving"
    export JAX_DISABLE_JIT=true
fi


/lustre/home/2101110354/anaconda3/envs/hoi_env/bin/python3 ${PYTHON_SCRIPT} \
    --input-file "${input_file}" \
    --output-dir "${output_dir}" \
    --sample-name "${sample_name}" \
    --n-bins 5 \
    --method ${METHOD} \
    --resume \
    --n-jobs ${CPUS_PER_TASK} \
    ${N_MAX} \
    --block-groups "${BLOCK_GROUPS}" \
    --group-separator "${GROUP_SEPARATOR}" \
    2>&1 | tee "${output_dir}/analysis.log"

echo "Job finished at: \$(date)"
EEOF

    echo "Submitting job for: ${sample_name}"
    sbatch "${output_dir}/job.sh"
    sleep 0.5
done

echo ""
echo "======================================================================"
echo "All entropy jobs submitted!"
echo "Check status with: squeue -u \$USER"
echo "Results will be saved to: ${RESULTS_DIR}"
echo "======================================================================"
