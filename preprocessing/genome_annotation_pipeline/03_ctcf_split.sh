#!/usr/bin/env bash
# 03 用 CTCF 结合与否划分 proximal / distal enhancer（仅 hg38 有 ENCODE enhancer）
set -euo pipefail
source "$(dirname "$0")/config.sh"
[ "$GENOME" = "hg38" ] || { echo "CTCF 划分仅适用于 hg38"; exit 0; }

for pair in "proximal:$ENH_PROX" "distal:$ENH_DIST"; do
    name=${pair%%:*}; bed=${pair#*:}
    bedtools intersect -a "$bed" -b "$CTCF" -u \
        > "$ANN_DIR/${name}_enhancer.ctcf_bound.bed"
    bedtools intersect -a "$bed" -b "$CTCF" -v \
        > "$ANN_DIR/${name}_enhancer.not_ctcf_bound.bed"
done
echo "[done] ctcf_bound / not_ctcf_bound splits in $ANN_DIR"
