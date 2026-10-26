#!/bin/bash
#SBATCH --job-name=preimp_triplet_null
#SBATCH --partition=C064M0256G
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --array=0-6
#SBATCH --output=slurm-%A_%a.out
#SBATCH --error=slurm-%A_%a.err
export EPINFO_ROOT="${EPINFO_ROOT:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && while [ ! -f .epigenetic_oinfo_root ]; do cd .. || exit 1; done && pwd)}"
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
cd "${EPINFO_ROOT}/cases/mouse_preimp"
line=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" synergy_redundancy/scripts/oinfo_samples.tsv)
IFS=$'\t' read -r sample input_file <<< "$line"
test -n "$sample"
out="${EPINFO_ROOT}/cases/mouse_preimp/synergy_redundancy/synergistic_triplets/permutation_results/$sample"
mkdir -p "$out"
echo "Started $(date -Is): $sample; job=$SLURM_JOB_ID"
python3 ${EPINFO_ROOT}/common/02_oinformation/compute_hoi_oinfo.py \
  --input-file "$input_file" --output-dir "$out" --sample-name "$sample" \
  --n-bins 5 --n-max 100000 --random-seed 42 --drop-signals '' \
  --n-permutations 10000 --permutation-seed 1234 --permutation-scope block --permutation-level all --permutation-max-order 3 --n-jobs "$SLURM_CPUS_PER_TASK" \
  --block-groups "${EPINFO_ROOT}/data/genome_features/mm10/promoter_bins.tsv;${EPINFO_ROOT}/data/genome_features/mm10/proximal_ctcf_bound_enhancer_purified_bins.tsv;${EPINFO_ROOT}/data/genome_features/mm10/proximal_not_ctcf_bound_enhancer_purified_bins.tsv;${EPINFO_ROOT}/data/genome_features/mm10/distal_ctcf_bound_enhancer_purified_bins.tsv;${EPINFO_ROOT}/data/genome_features/mm10/distal_not_ctcf_bound_enhancer_purified_bins.tsv;${EPINFO_ROOT}/data/genome_features/mm10/exon_purified_bins.tsv;${EPINFO_ROOT}/data/genome_features/mm10/intron_purified_bins.tsv::BLOCKGROUP::${EPINFO_ROOT}/data/genome_features/mm10/CpG_island_bins.tsv;${EPINFO_ROOT}/data/genome_features/mm10/LINE_bins.tsv;${EPINFO_ROOT}/data/genome_features/mm10/SINE_bins.tsv;${EPINFO_ROOT}/data/genome_features/mm10/LTR_bins.tsv" --group-separator '::BLOCKGROUP::'
echo "Finished $(date -Is): $sample"
