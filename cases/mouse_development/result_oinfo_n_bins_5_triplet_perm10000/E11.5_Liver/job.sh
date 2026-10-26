#!/usr/bin/env bash
export EPINFO_ROOT="${EPINFO_ROOT:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && while [ ! -f .epigenetic_oinfo_root ]; do cd .. || exit 1; done && pwd)}"
set -euo pipefail
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
echo "Started $(date --iso-8601=seconds); sample=E11.5_Liver; job=${SLURM_JOB_ID:-local}"
python3 ${EPINFO_ROOT}/cases/mouse_development/result_oinfo_n_bins_5_triplet_perm10000/code/compute_hoi_oinfo.py --input-file ${EPINFO_ROOT}/cases/mouse_development/data/input_files/E11.5_Liver.txt --output-dir ${EPINFO_ROOT}/cases/mouse_development/result_oinfo_n_bins_5_triplet_perm10000/E11.5_Liver --sample-name E11.5_Liver --n-bins 5 --n-max 100000 --random-seed 42 --n-jobs 32 --block-groups "${EPINFO_ROOT}/data/genome_features/mm10/promoter_bins.tsv;${EPINFO_ROOT}/data/genome_features/mm10/proximal_ctcf_bound_enhancer_purified_bins.tsv;${EPINFO_ROOT}/data/genome_features/mm10/proximal_not_ctcf_bound_enhancer_purified_bins.tsv;${EPINFO_ROOT}/data/genome_features/mm10/distal_ctcf_bound_enhancer_purified_bins.tsv;${EPINFO_ROOT}/data/genome_features/mm10/distal_not_ctcf_bound_enhancer_purified_bins.tsv;${EPINFO_ROOT}/data/genome_features/mm10/exon_purified_bins.tsv;${EPINFO_ROOT}/data/genome_features/mm10/intron_purified_bins.tsv::BLOCKGROUP::${EPINFO_ROOT}/data/genome_features/mm10/CpG_island_bins.tsv;${EPINFO_ROOT}/data/genome_features/mm10/LINE_bins.tsv;${EPINFO_ROOT}/data/genome_features/mm10/SINE_bins.tsv;${EPINFO_ROOT}/data/genome_features/mm10/LTR_bins.tsv" --drop-signals '' --n-permutations 10000 --permutation-seed 1234 --permutation-scope block --permutation-level all --permutation-min-order 3 --permutation-max-order 3
echo "Completed $(date --iso-8601=seconds)"
