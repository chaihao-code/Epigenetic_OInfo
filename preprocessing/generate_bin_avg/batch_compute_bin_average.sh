#!/usr/bin/env bash
set -euo pipefail

BW_DIR="$1"
OUT_DIR="$2"
CHROM_SIZES="$3"
BIN_SIZE="${4:-1000}"
STEP_SIZE="${5:-$BIN_SIZE}"

mkdir -p logs

shopt -s nullglob
for BW in "$BW_DIR"/*.bw "$BW_DIR"/*.bigWig "$BW_DIR"/*.bigwig "$BW_DIR"/*.bedGraph; do
    echo "Submitting job for $(basename "$BW")"
    sbatch /lustre/home/2101110354/sequence_entropy/script/preprocessing/generate_bin_avg/slurm_compute_bin_average.sh \
        "$BW" \
        "$CHROM_SIZES" \
        "$OUT_DIR" \
        "$BIN_SIZE" \
        "$STEP_SIZE"
done
