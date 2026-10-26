#!/usr/bin/env bash

# ==============================
# Convert per-CpG BED to bedGraph
# ==============================
#
# 输入：
#   IN_BED    : 每一行对应一个 cytosine 位点的 BED 文件
#               （通常来自 Bismark / NOMe-seq）
#
# 输出：
#   OUT_BG    : bedGraph 文件
#               格式：
#               chr  start  end  methylation_fraction
#
# 参数：
#   THRESHOLD : 最小 coverage（推荐 8 或 10）
#
# 功能：
#   1. 过滤低 coverage CpG
#   2. 计算 methylation fraction = meth / coverage
#   3. 输出为连续信号用于：
#        - bedtools map
#        - genome bin averaging
#        - entropy / MI analysis
#
# ==============================

set -euo pipefail

IN_BED="$1"
OUT_BG="$2"
THRESHOLD="${3:-8}"   # 默认 coverage >= 8

echo "[INFO] Input BED       : $IN_BED"
echo "[INFO] Output bedGraph : $OUT_BG"
echo "[INFO] Coverage cutoff : $THRESHOLD"

awk -v threshold="$THRESHOLD" '
BEGIN { OFS="\t" }
{
    cov  = $10      # 真正的 coverage
    meth = $11      # 11列是百分比, 需要转成0-1的小数
    if (cov >= threshold) {
        print $1, $2, $3, meth
    }
}' "$IN_BED" > "$OUT_BG"

echo "[DONE] bedGraph written to $OUT_BG"
