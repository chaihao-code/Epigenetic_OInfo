#!/usr/bin/env bash
export EPINFO_ROOT="${EPINFO_ROOT:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && while [ ! -f .epigenetic_oinfo_root ]; do cd .. || exit 1; done && pwd)}"
set -euo pipefail

BW_DIR="${EPINFO_ROOT}/data/mouse_preimplementation/original/methylation"
OUT_DIR="${EPINFO_ROOT}/data/mouse_preimplementation/original/methylation_rep_averaged"
CHROM_SIZES="${EPINFO_ROOT}/common/01_preprocessing/average_replicates_and_liftover/mm9.chrom.sizes"

mkdir -p logs ${OUT_DIR}

groups=$(
    ls ${BW_DIR}/*.bw \
    | grep '_single_cell_' \
    | sed -E 's#.*/[^_]+_(.+)_single_cell_[0-9]+\.([^.]+(\.[^.]+)*)\.bw#\1.\2#' \
    | sort | uniq
)

for g in ${groups}; do
    GROUP=${g%%.*}
    KMER=${g#*.}

    echo "Submitting ${GROUP}.${KMER}"

    sbatch \
        -J avg_${GROUP}_${KMER} \
        average_one_replicate_group.sh \
        ${BW_DIR} \
        ${OUT_DIR} \
        ${GROUP} \
        ${KMER}
done

