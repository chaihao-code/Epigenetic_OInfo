#!/usr/bin/env bash
set -euo pipefail

BW_DIR="/lustre/home/2101110354/sequence_entropy/data/mouse_preimplementation/original/methylation"
OUT_DIR="/lustre/home/2101110354/sequence_entropy/data/mouse_preimplementation/original/methylation_rep_averaged"
CHROM_SIZES="/lustre/home/2101110354/sequence_entropy/script/preprocessing/average_replicates_and_liftover/mm9.chrom.sizes"

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

