#!/bin/bash
#SBATCH --job-name=mouse_preimp_oinfo
#SBATCH --partition=C064M0256G
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --array=0-6
#SBATCH --output=/lustre/home/2101110354/sequence_entropy/project/mouse_preimp/synergy_redundancy/results_n_bins_5/logs/%A_%a.out
#SBATCH --error=/lustre/home/2101110354/sequence_entropy/project/mouse_preimp/synergy_redundancy/results_n_bins_5/logs/%A_%a.err
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
cd "/lustre/home/2101110354/sequence_entropy/project/mouse_preimp"
line=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" synergy_redundancy/scripts/oinfo_samples.tsv)
IFS=$'\t' read -r sample input_file <<< "$line"
test -n "$sample"
out="/lustre/home/2101110354/sequence_entropy/project/mouse_preimp/synergy_redundancy/results_n_bins_5/$sample"
mkdir -p "$out"
echo "Started $(date -Is): $sample; job=$SLURM_JOB_ID"
/lustre/home/2101110354/anaconda3/envs/hoi_env/bin/python /lustre/home/2101110354/sequence_entropy/script/hoi_code/compute_hoi_oinfo.py \
  --input-file "$input_file" --output-dir "$out" --sample-name "$sample" \
  --n-bins 5 --n-max 100000 --random-seed 42 --drop-signals '' \
  --n-permutations 10000 --permutation-seed 1234 --permutation-scope block --permutation-level system --n-jobs "$SLURM_CPUS_PER_TASK" \
  --block-groups '/lustre/home/2101110354/sequence_entropy/data/genome_features/mm10/promoter_bins.tsv;/lustre/home/2101110354/sequence_entropy/data/genome_features/mm10/proximal_ctcf_bound_enhancer_purified_bins.tsv;/lustre/home/2101110354/sequence_entropy/data/genome_features/mm10/proximal_not_ctcf_bound_enhancer_purified_bins.tsv;/lustre/home/2101110354/sequence_entropy/data/genome_features/mm10/distal_ctcf_bound_enhancer_purified_bins.tsv;/lustre/home/2101110354/sequence_entropy/data/genome_features/mm10/distal_not_ctcf_bound_enhancer_purified_bins.tsv;/lustre/home/2101110354/sequence_entropy/data/genome_features/mm10/exon_purified_bins.tsv;/lustre/home/2101110354/sequence_entropy/data/genome_features/mm10/intron_purified_bins.tsv::BLOCKGROUP::/lustre/home/2101110354/sequence_entropy/data/genome_features/mm10/CpG_island_bins.tsv;/lustre/home/2101110354/sequence_entropy/data/genome_features/mm10/LINE_bins.tsv;/lustre/home/2101110354/sequence_entropy/data/genome_features/mm10/SINE_bins.tsv;/lustre/home/2101110354/sequence_entropy/data/genome_features/mm10/LTR_bins.tsv' --group-separator '::BLOCKGROUP::'
echo "Finished $(date -Is): $sample"
