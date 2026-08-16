#!/usr/bin/env bash
# ============================================================================
# 程序名称：WGBS 甲基化 coverage 加权分箱
# 功能描述：遍历输入目录中所有成对的 WGB-Seq coverage / methylation bigWig，
#          按 coverage 阈值过滤后，生成 (methylation * coverage) 的乘积 bigWig，
#          并对 coverage 与乘积分别做 1 kb（可自定义）分箱平均。
#          最终用户可以自行将乘积分箱结果除以 coverage 分箱结果，得到加权甲基化。
#
# 输入参数：
#   $1: 输入目录（包含 *.WGB-Seq.coverage.bigWig 与 *.WGB-Seq.methylation.bigWig）
#   $2: 染色体大小文件 chrom.sizes
#   $3: 输出目录
#   $4: bin size（默认 1000 bp）
#   $5: step size（默认等于 bin size）
#   $6: coverage 阈值（默认 5，低于该值的 CpG 会被过滤）
#
# 输出文件（每个样本）：
#   ${prefix}.WGB-Seq.methXcov.thresh${T}.bigWig          # 乘积 bigWig
#   ${prefix}.WGB-Seq.coverage.thresh${T}.bigWig          # 过滤后的 coverage bigWig
#   ${prefix}.WGB-Seq.methXcov.thresh${T}.binSIZE.stepSIZE.bp   # 乘积分箱平均
#   ${prefix}.WGB-Seq.coverage.thresh${T}.binSIZE.stepSIZE.bp   # coverage 分箱平均
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

IN_DIR="$1"
CHROM_SIZES="$2"
OUT_DIR="$3"
BIN_SIZE="${4:-1000}"
STEP_SIZE="${5:-$BIN_SIZE}"
COV_THRESHOLD="${6:-5}"

mkdir -p "$OUT_DIR"

# 优先使用项目里已有的 UCSC 二进制工具（避免 anaconda 中 libssl 不兼容问题）
UCSC_DIR="${UCSC_TOOLS:-/media/dell/data1/chaihao/ucsc_binary_tools}"
if [[ -d "$UCSC_DIR" ]]; then
    export PATH="$UCSC_DIR:$PATH"
fi

# 依赖检查
for cmd in bigWigToBedGraph bedGraphToBigWig bedtools bigWigAverageOverBed; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "[ERROR] Required command not found: $cmd" >&2
        exit 1
    fi
done

# 按样本核心 ID（Blueprint.IHECRE...WGB-Seq）配对 coverage 和 methylation
# 注意：不同文件类型的数字前缀可能不同，例如
#   126764.Blueprint.IHECRE00000004.3.WGB-Seq.coverage.bigWig
#   126765.Blueprint.IHECRE00000004.3.WGB-Seq.methylation.bigWig
# 因此不能简单替换后缀，需要用去掉数字前缀后的核心字符串配对。
get_core() {
    local f
    f=$(basename "$1")
    # 去掉开头的 "数字."
    f="${f#[0-9]*.}"
    # 去掉结尾的 .methylation.bigWig 或 .coverage.bigWig
    f="${f%.methylation.bigWig}"
    f="${f%.coverage.bigWig}"
    printf '%s' "$f"
}

shopt -s nullglob

PAIRED=0
for METH_BW in "$IN_DIR"/*.WGB-Seq.methylation.bigWig; do
    [[ -f "$METH_BW" ]] || continue

    METH_CORE=$(get_core "$METH_BW")
    COV_BW=""
    for c in "$IN_DIR"/*.WGB-Seq.coverage.bigWig; do
        [[ -f "$c" ]] || continue
        if [[ "$(get_core "$c")" == "$METH_CORE" ]]; then
            COV_BW="$c"
            break
        fi
    done

    echo "$COV_BW"
    if [[ -z "$COV_BW" ]]; then
        echo "[WARN] No matching coverage file for $METH_BW (core: $METH_CORE), skipping" >&2
        continue
    fi

    PAIRED=$((PAIRED + 1))

    # prefix 使用 methylation 文件的名称（去掉后缀），保持与单文件版一致
    BASENAME=$(basename "$METH_BW" .methylation.bigWig)
    echo "[INFO] Processing $BASENAME ..."
    echo "[INFO]   Coverage   : $COV_BW"
    echo "[INFO]   Methylation: $METH_BW"

    TMPD=$(mktemp -d)

    COV_BG="$TMPD/cov.bg"
    METH_BG="$TMPD/meth.bg"
    UNION_BG="$TMPD/union.bg"
    PROD_BG="$TMPD/prod.bg"
    COV_FILT_BG="$TMPD/cov_filt.bg"

    PROD_BW="$OUT_DIR/${BASENAME}.methXcov.thresh${COV_THRESHOLD}.bigWig"
    COV_FILT_BW="$OUT_DIR/${BASENAME}.coverage.thresh${COV_THRESHOLD}.bigWig"

    # 1. bigWig -> bedGraph
    echo "[INFO]   Converting bigWig to bedGraph ..."
    bigWigToBedGraph "$COV_BW" "$COV_BG"
    bigWigToBedGraph "$METH_BW" "$METH_BG"

    # 2. unionbedg 把 coverage 和 methylation 按坐标对齐
    #    输出：chrom start end cov meth
    echo "[INFO]   Aligning coverage and methylation by coordinate ..."
    bedtools unionbedg -g "$CHROM_SIZES" -i "$COV_BG" "$METH_BG" > "$UNION_BG"

    # 3. 按 coverage 阈值过滤，并生成乘积 bedGraph 和过滤后的 coverage bedGraph
    echo "[INFO]   Filtering CpGs with coverage >= $COV_THRESHOLD and computing meth*cov ..."
    awk -v t="$COV_THRESHOLD" 'BEGIN{OFS="\t"} ($4 >= t) {print $1, $2, $3, $4 * $5}' "$UNION_BG" > "$PROD_BG"
    awk -v t="$COV_THRESHOLD" 'BEGIN{OFS="\t"} ($4 >= t) {print $1, $2, $3, $4}' "$UNION_BG" > "$COV_FILT_BG"

    # 4. bedGraph -> bigWig
    echo "[INFO]   Writing product and filtered coverage bigWigs ..."
    bedGraphToBigWig "$PROD_BG" "$CHROM_SIZES" "$PROD_BW"
    bedGraphToBigWig "$COV_FILT_BG" "$CHROM_SIZES" "$COV_FILT_BW"

    # 5. 用 compute_bin_average.sh 计算 1 kb（或自定义）分箱平均
    echo "[INFO]   Computing bin averages (bin=$BIN_SIZE, step=$STEP_SIZE) ..."
    bash "$SCRIPT_DIR/compute_bin_average.sh" \
        "$PROD_BW" "$CHROM_SIZES" "$OUT_DIR" "$BIN_SIZE" "$STEP_SIZE"

    bash "$SCRIPT_DIR/compute_bin_average.sh" \
        "$COV_FILT_BW" "$CHROM_SIZES" "$OUT_DIR" "$BIN_SIZE" "$STEP_SIZE"

    # 清理本样本临时文件
    rm -rf "$TMPD"

done

if [[ $PAIRED -eq 0 ]]; then
    echo "[WARN] No paired WGB-Seq coverage/methylation files found in $IN_DIR"
    exit 0
fi

echo "[DONE] All samples processed. Results in $OUT_DIR"
