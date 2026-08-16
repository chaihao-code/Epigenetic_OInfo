#!/usr/bin/env bash
#SBATCH --job-name=cpg2bg
#SBATCH -p cpu1,cpu2,cpu_short
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

IN_BED_GZ="$1"
OUT_DIR="$2"
THRESHOLD="$3"

# SCRIPT="/lustre/grp/gyqlab/chaih/sequence_entropy/code/generate_bin_avg/convert_cpg_to_bedgraph.sh"

BASENAME=$(basename "$IN_BED_GZ" .bed.gz)
OUT_BG="${OUT_DIR}/${BASENAME}.bedGraph"

echo "[INFO] Processing $IN_BED_GZ"

mkdir -p "$OUT_DIR"

# 流式解压 + 转换（不产生中间文件）
zcat "$IN_BED_GZ" | \
awk -v threshold="$THRESHOLD" '
BEGIN{OFS="\t"}
{
    cov=$5
    meth=$11
    if (cov >= threshold) {
        print $1,$2,$3,meth
    }
}' > "$OUT_BG"

echo "[DONE] $OUT_BG"
