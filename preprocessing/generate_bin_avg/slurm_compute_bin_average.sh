#!/bin/bash
#SBATCH -p C064M1024G
#SBATCH -J bw2bin
#SBATCH -c 1
#SBATCH -t 24:00:00
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err

module load deeptools
module load bedtools

BW_FILE="$1"
CHROM_SIZES="$2"
OUT_DIR="$3"
BIN_SIZE="${4:-1000}"
STEP_SIZE="${5:-$BIN_SIZE}"

bash /lustre/home/2101110354/sequence_entropy/script/preprocessing/generate_bin_avg/compute_bin_average.sh "$BW_FILE" "$CHROM_SIZES" "$OUT_DIR" "$BIN_SIZE" "$STEP_SIZE"
