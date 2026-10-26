#!/bin/bash
#SBATCH --job-name=ent_n10_zygote
#SBATCH --partition=C064M0256G
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=/lustre/home/2101110354/sequence_entropy/project/mouse_preimp/calculation_results/results_entropy_n_bins_10/zygote/slurm_%j.out
#SBATCH --error=/lustre/home/2101110354/sequence_entropy/project/mouse_preimp/calculation_results/results_entropy_n_bins_10/zygote/slurm_%j.err

echo "Job started at: $(date)"
export JAX_PLATFORM_NAME=cpu
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_PYTHON_CLIENT_ALLOCATOR=platform
export XLA_FLAGS="--xla_cpu_use_thunk_runtime=false"
export JAX_DISABLE_JIT=true

/lustre/home/2101110354/anaconda3/envs/hoi_env/bin/python3 /lustre/home/2101110354/sequence_entropy/script/hoi_code/compute_hoi_entropy.py     --input-file "/lustre/home/2101110354/sequence_entropy/project/mouse_preimp/calculation_results/data_info/zygote.txt"     --output-dir "/lustre/home/2101110354/sequence_entropy/project/mouse_preimp/calculation_results/results_entropy_n_bins_10/zygote"     --sample-name "zygote"     --n-bins 10     --method binning     --resume     --n-jobs 32     --n-max 100000     --block-groups "/lustre/home/2101110354/sequence_entropy/data/genome_features/mm10/promoter_bins.tsv;/lustre/home/2101110354/sequence_entropy/data/genome_features/mm10/proximal_ctcf_bound_enhancer_purified_bins.tsv;/lustre/home/2101110354/sequence_entropy/data/genome_features/mm10/proximal_not_ctcf_bound_enhancer_purified_bins.tsv;/lustre/home/2101110354/sequence_entropy/data/genome_features/mm10/distal_ctcf_bound_enhancer_purified_bins.tsv;/lustre/home/2101110354/sequence_entropy/data/genome_features/mm10/distal_not_ctcf_bound_enhancer_purified_bins.tsv;/lustre/home/2101110354/sequence_entropy/data/genome_features/mm10/exon_purified_bins.tsv;/lustre/home/2101110354/sequence_entropy/data/genome_features/mm10/intron_purified_bins.tsv::BLOCKGROUP::/lustre/home/2101110354/sequence_entropy/data/genome_features/mm10/CpG_island_bins.tsv;/lustre/home/2101110354/sequence_entropy/data/genome_features/mm10/LINE_bins.tsv;/lustre/home/2101110354/sequence_entropy/data/genome_features/mm10/SINE_bins.tsv;/lustre/home/2101110354/sequence_entropy/data/genome_features/mm10/LTR_bins.tsv"     --group-separator "::BLOCKGROUP::"     2>&1 | tee "/lustre/home/2101110354/sequence_entropy/project/mouse_preimp/calculation_results/results_entropy_n_bins_10/zygote/analysis.log"

# Merge per-block entropy_all_subsets.csv into one stage-level csv
/lustre/home/2101110354/anaconda3/envs/hoi_env/bin/python3 - <<'PYEOF'
import glob, os
out = "/lustre/home/2101110354/sequence_entropy/project/mouse_preimp/calculation_results/results_entropy_n_bins_10/zygote/zygote_entropy_all_subsets.csv"
files = sorted(glob.glob("/lustre/home/2101110354/sequence_entropy/project/mouse_preimp/calculation_results/results_entropy_n_bins_10/zygote/*/*_entropy_all_subsets.csv"))
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
/lustre/home/2101110354/anaconda3/envs/hoi_env/bin/python3 /lustre/home/2101110354/sequence_entropy/script/hoi_code/convert_entropy_to_mi.py -i "/lustre/home/2101110354/sequence_entropy/project/mouse_preimp/calculation_results/results_entropy_n_bins_10" -o "/lustre/home/2101110354/sequence_entropy/project/mouse_preimp/mutual_information/pairwise_mi_n10"

echo "Job finished at: $(date)"
