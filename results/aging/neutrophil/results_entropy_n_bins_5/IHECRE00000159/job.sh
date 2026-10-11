#!/bin/bash
#SBATCH --job-name=entropy_IHECRE00000159
#SBATCH --partition=C064M0256G
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=/lustre/home/2101110354/sequence_entropy/project/aging/neutrophil/results_entropy_n_bins_5/IHECRE00000159/slurm_%j.out
#SBATCH --error=/lustre/home/2101110354/sequence_entropy/project/aging/neutrophil/results_entropy_n_bins_5/IHECRE00000159/slurm_%j.err

echo "Job started at: $(date)"
echo "Running on host: $(hostname)"
echo "SLURM job ID: $SLURM_JOB_ID"

# Set JAX memory limits
export JAX_PLATFORM_NAME=cpu
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export XLA_PYTHON_CLIENT_ALLOCATOR=platform
export XLA_FLAGS="--xla_cpu_use_thunk_runtime=false"

if [ "1" == "1" ]; then
    echo "WARNING: JAX JIT disabled for memory saving"
    export JAX_DISABLE_JIT=true
fi


/lustre/home/2101110354/anaconda3/envs/hoi_env/bin/python3 /lustre/home/2101110354/sequence_entropy/script/hoi_code/compute_hoi_entropy.py     --input-file "/lustre/home/2101110354/sequence_entropy/project/aging/neutrophil/data_info/IHECRE00000159.txt"     --output-dir "/lustre/home/2101110354/sequence_entropy/project/aging/neutrophil/results_entropy_n_bins_5/IHECRE00000159"     --sample-name "IHECRE00000159"     --n-bins 5     --method binning     --resume     --n-jobs 32     --n-max 100000     --block-groups "/lustre/home/2101110354/sequence_entropy/data/genome_features/hg38/promoter_bins.tsv;/lustre/home/2101110354/sequence_entropy/data/genome_features/hg38/proximal_ctcf_bound_enhancer_purified_bins.tsv;/lustre/home/2101110354/sequence_entropy/data/genome_features/hg38/proximal_not_ctcf_bound_enhancer_purified_bins.tsv;/lustre/home/2101110354/sequence_entropy/data/genome_features/hg38/distal_ctcf_bound_enhancer_purified_bins.tsv;/lustre/home/2101110354/sequence_entropy/data/genome_features/hg38/distal_not_ctcf_bound_enhancer_purified_bins.tsv;/lustre/home/2101110354/sequence_entropy/data/genome_features/hg38/exon_purified_bins.tsv;/lustre/home/2101110354/sequence_entropy/data/genome_features/hg38/intron_purified_bins.tsv::BLOCKGROUP::/lustre/home/2101110354/sequence_entropy/data/genome_features/hg38/CpG_island_bins.tsv;/lustre/home/2101110354/sequence_entropy/data/genome_features/hg38/LINE_bins.tsv;/lustre/home/2101110354/sequence_entropy/data/genome_features/hg38/SINE_bins.tsv;/lustre/home/2101110354/sequence_entropy/data/genome_features/hg38/LTR_bins.tsv"     --group-separator "::BLOCKGROUP::"     2>&1 | tee "/lustre/home/2101110354/sequence_entropy/project/aging/neutrophil/results_entropy_n_bins_5/IHECRE00000159/analysis.log"

echo "Job finished at: $(date)"
