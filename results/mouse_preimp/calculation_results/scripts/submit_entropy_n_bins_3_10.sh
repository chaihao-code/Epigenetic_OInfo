#!/bin/bash
# Submit entropy computation for n_bins=3 and n_bins=10 (mouse preimplantation)
# Each job: compute_hoi_entropy.py -> merge block csvs -> convert_entropy_to_mi.py
set -e

PROJECT_DIR="/lustre/home/2101110354/sequence_entropy/project/mouse_preimp/calculation_results"
INPUT_DIR="${PROJECT_DIR}/data_info"
PYTHON_SCRIPT="/lustre/home/2101110354/sequence_entropy/script/hoi_code/compute_hoi_entropy.py"
CONVERT_SCRIPT="/lustre/home/2101110354/sequence_entropy/script/hoi_code/convert_entropy_to_mi.py"
PY="/lustre/home/2101110354/anaconda3/envs/hoi_env/bin/python3"

GENOME_FEATURES_DIR="/lustre/home/2101110354/sequence_entropy/data/genome_features/mm10"
BLOCK_GROUP_1="${GENOME_FEATURES_DIR}/promoter_bins.tsv;${GENOME_FEATURES_DIR}/proximal_ctcf_bound_enhancer_purified_bins.tsv;${GENOME_FEATURES_DIR}/proximal_not_ctcf_bound_enhancer_purified_bins.tsv;${GENOME_FEATURES_DIR}/distal_ctcf_bound_enhancer_purified_bins.tsv;${GENOME_FEATURES_DIR}/distal_not_ctcf_bound_enhancer_purified_bins.tsv;${GENOME_FEATURES_DIR}/exon_purified_bins.tsv;${GENOME_FEATURES_DIR}/intron_purified_bins.tsv"
BLOCK_GROUP_2="${GENOME_FEATURES_DIR}/CpG_island_bins.tsv;${GENOME_FEATURES_DIR}/LINE_bins.tsv;${GENOME_FEATURES_DIR}/SINE_bins.tsv;${GENOME_FEATURES_DIR}/LTR_bins.tsv"
GROUP_SEPARATOR="::BLOCKGROUP::"

PARTITION="C064M0256G"
TIME="04:00:00"
CPUS=32

STAGES="zygote 2cell 4cell 8cell morula ICM TE"
N_BINS_LIST="3 10"

for N_BINS in ${N_BINS_LIST}; do
    RESULTS_DIR="${PROJECT_DIR}/results_entropy_n_bins_${N_BINS}"
    MI_DIR="${PROJECT_DIR}/../mutual_information/pairwise_mi_n${N_BINS}"
    for stage in ${STAGES}; do
        input_file="${INPUT_DIR}/${stage}.txt"
        output_dir="${RESULTS_DIR}/${stage}"
        job_name="ent_n${N_BINS}_${stage}"
        mkdir -p "${output_dir}"
        cat > "${output_dir}/job.sh" << EEOF
#!/bin/bash
#SBATCH --job-name=${job_name}
#SBATCH --partition=${PARTITION}
#SBATCH --time=${TIME}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --output=${output_dir}/slurm_%j.out
#SBATCH --error=${output_dir}/slurm_%j.err

echo "Job started at: \$(date)"
export JAX_PLATFORM_NAME=cpu
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_PYTHON_CLIENT_ALLOCATOR=platform
export XLA_FLAGS="--xla_cpu_use_thunk_runtime=false"
export JAX_DISABLE_JIT=true

${PY} ${PYTHON_SCRIPT} \
    --input-file "${input_file}" \
    --output-dir "${output_dir}" \
    --sample-name "${stage}" \
    --n-bins ${N_BINS} \
    --method binning \
    --resume \
    --n-jobs ${CPUS} \
    --n-max 100000 \
    --block-groups "${BLOCK_GROUP_1}${GROUP_SEPARATOR}${BLOCK_GROUP_2}" \
    --group-separator "${GROUP_SEPARATOR}" \
    2>&1 | tee "${output_dir}/analysis.log"

# Merge per-block entropy_all_subsets.csv into one stage-level csv
${PY} - <<'PYEOF'
import glob, os
out = "${output_dir}/${stage}_entropy_all_subsets.csv"
files = sorted(glob.glob("${output_dir}/*/*_entropy_all_subsets.csv"))
with open(out, "w") as fo:
    first = True
    for f in files:
        with open(f) as fi:
            header = fi.readline()
            if first:
                fo.write(header)
                first = False
            fo.write(fi.read())
print(f"Merged {len(files)} block files -> {out}")
PYEOF

# Convert to MI format
${PY} ${CONVERT_SCRIPT} -i "${RESULTS_DIR}" -o "${MI_DIR}"

echo "Job finished at: \$(date)"
EEOF
        echo "Submitting: n_bins=${N_BINS} stage=${stage}"
        sbatch "${output_dir}/job.sh"
        sleep 0.5
    done
done

echo "All jobs submitted."
