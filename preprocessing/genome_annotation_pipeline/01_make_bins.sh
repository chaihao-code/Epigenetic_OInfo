#!/usr/bin/env bash
# 01 生成 1-kb bins
set -euo pipefail
source "$(dirname "$0")/config.sh"

bedtools makewindows -g "$CHROM_SIZES" -w "$BIN_SIZE" > "$BINS"
echo "[done] $BINS : $(wc -l < "$BINS") bins"
