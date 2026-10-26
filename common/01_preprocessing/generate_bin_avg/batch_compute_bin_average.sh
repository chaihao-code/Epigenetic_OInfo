#!/usr/bin/env bash
export EPINFO_ROOT="${EPINFO_ROOT:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && while [ ! -f .epigenetic_oinfo_root ]; do cd .. || exit 1; done && pwd)}"
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
    sbatch ${EPINFO_ROOT}/common/01_preprocessing/generate_bin_avg/slurm_compute_bin_average.sh \
        "$BW" \
        "$CHROM_SIZES" \
        "$OUT_DIR" \
        "$BIN_SIZE" \
        "$STEP_SIZE"
done
