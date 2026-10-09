#!/bin/bash
# =============================================================================
# 批量提交 pairwise correlation 任务到 SLURM
# 处理项目：mouse_after_E11.5 和 mouse_preimp
# 复用两个项目 entropy 脚本中的 input files、block groups 和 n-bins=5 配置
# 默认：离散模式 + 按 genome region 笛卡尔积分块计算（输出 4 列）
# 可选：--overall 只输出总体相关（3 列）
# =============================================================================

set -euo pipefail

# -----------------------------------------------------------------------------
# 配置
# -----------------------------------------------------------------------------
CORR_SCRIPT="/lustre/home/2101110354/sequence_entropy/script/hoi_code/compute_signal_correlation.py"
GENOME_FEATURES_DIR="/lustre/home/2101110354/sequence_entropy/data/genome_features/mm10"
GROUP_SEPARATOR="::BLOCKGROUP::"

# 与 entropy 脚本相同的 block groups
BLOCK_GROUP_1="${GENOME_FEATURES_DIR}/promoter_bins.tsv;${GENOME_FEATURES_DIR}/proximal_ctcf_bound_enhancer_purified_bins.tsv;${GENOME_FEATURES_DIR}/proximal_not_ctcf_bound_enhancer_purified_bins.tsv;${GENOME_FEATURES_DIR}/distal_ctcf_bound_enhancer_purified_bins.tsv;${GENOME_FEATURES_DIR}/distal_not_ctcf_bound_enhancer_purified_bins.tsv;${GENOME_FEATURES_DIR}/exon_purified_bins.tsv;${GENOME_FEATURES_DIR}/intron_purified_bins.tsv"
BLOCK_GROUP_2="${GENOME_FEATURES_DIR}/CpG_island_bins.tsv;${GENOME_FEATURES_DIR}/LINE_bins.tsv;${GENOME_FEATURES_DIR}/SINE_bins.tsv;${GENOME_FEATURES_DIR}/LTR_bins.tsv"
BLOCK_GROUPS="${BLOCK_GROUP_1}${GROUP_SEPARATOR}${BLOCK_GROUP_2}"

# SLURM 资源（correlation 单线程，资源比 entropy 小）
PARTITION="C064M0256G"
TIME="00:10:00"
CPUS_PER_TASK=1

# 输出目录名
OUTPUT_DIR_NAME="results_pairwise_correlation_n_bins_5"

# 要处理的项目及其 input 目录
# 格式：project_dir|relative_input_dir
PROJECT_CONFIGS=(
    "/lustre/home/2101110354/sequence_entropy/project/mouse_after_E11.5|data/input_files"
    "/lustre/home/2101110354/sequence_entropy/project/mouse_preimp|data_info"
)

# 模式：默认 per-block，--overall 为总体相关
OVERALL=0
if [ "${1:-}" == "--overall" ]; then
    OVERALL=1
    echo "Mode: overall correlation (3-column output)"
else
    echo "Mode: per-block correlation (4-column output: block feature1 feature2 correlation)"
fi

# -----------------------------------------------------------------------------
# 检查 sbatch
# -----------------------------------------------------------------------------
if ! command -v sbatch &> /dev/null; then
    echo "Error: sbatch command not found. Are you on a SLURM cluster?"
    exit 1
fi

# -----------------------------------------------------------------------------
# 提交任务
# -----------------------------------------------------------------------------
total_submitted=0
total_skipped=0

for config in "${PROJECT_CONFIGS[@]}"; do
    IFS='|' read -r project_dir input_subdir <<< "$config"
    input_dir="${project_dir}/${input_subdir}"
    output_dir="${project_dir}/${OUTPUT_DIR_NAME}"
    logs_dir="${output_dir}/logs"

    mkdir -p "${logs_dir}"

    echo ""
    echo "Project: ${project_dir}"
    echo "  Input:  ${input_dir}"
    echo "  Output: ${output_dir}"

    shopt -s nullglob
    input_files=("${input_dir}"/*.txt)
    shopt -u nullglob

    if [ ${#input_files[@]} -eq 0 ]; then
        echo "  Warning: no .txt files found in ${input_dir}, skipping"
        continue
    fi

    for input_file in "${input_files[@]}"; do
        sample_name=$(basename "${input_file}" .txt)
        output_file="${output_dir}/${sample_name}.tsv"
        job_script="${logs_dir}/job_${sample_name}.sh"

        # resume：已存在则跳过
        if [ -f "${output_file}" ]; then
            echo "  [skip] ${sample_name} (output exists)"
            ((total_skipped++)) || true
            continue
        fi

        if [ "${OVERALL}" -eq 1 ]; then
            extra_args=""
            job_name="corr_overall_${sample_name}"
        else
            extra_args="--per-block --block-groups \"${BLOCK_GROUPS}\" --group-separator \"${GROUP_SEPARATOR}\""
            job_name="corr_block_${sample_name}"
        fi

        cat > "${job_script}" << EOF
#!/bin/bash
#SBATCH --job-name=${job_name}
#SBATCH --partition=${PARTITION}
#SBATCH --time=${TIME}
#SBATCH --cpus-per-task=${CPUS_PER_TASK}
#SBATCH --output=${logs_dir}/${sample_name}_%j.out
#SBATCH --error=${logs_dir}/${sample_name}_%j.err

/lustre/home/2101110354/anaconda3/envs/hoi_env/bin/python3.11 "${CORR_SCRIPT}" \
    -i "${input_file}" \
    -o "${output_file}" \
    -n 5 \
    ${extra_args}
EOF

        echo "  [submit] ${sample_name}"
        sbatch "${job_script}"
        ((total_submitted++)) || true
        sleep 0.2
    done
done

# -----------------------------------------------------------------------------
# 总结
# -----------------------------------------------------------------------------
echo ""
echo "======================================================================"
echo "All correlation jobs submitted!"
echo "  Submitted: ${total_submitted}"
echo "  Skipped:   ${total_skipped}"
echo "Monitor: squeue -u \$USER"
echo "======================================================================"
