#!/usr/bin/env bash
# 程序名称：基因组信号数据分箱处理工具
# 功能描述：将基因组连续信号数据（bigWig 或 bedGraph 格式）按照固定窗口大小（bin）进行分箱，
#           计算每个窗口内的平均信号值，生成标准化的 bed 格式文件。
# 主要用途：用于基因组学数据分析，如 ChIP-seq、ATAC-seq、RNA-seq 等信号数据的标准化处理，
#           便于后续的可视化、比较和统计分析。
# 输入参数：
#   $1: 输入文件路径（支持 .bw/.bigwig/.bigWig/.bedGraph 格式）
#   $2: 染色体大小文件路径（格式：染色体名称<tab>染色体长度）
#   $3: 输出目录路径
# 输出文件：${输出目录}/${输入文件名}.bin1000.bp
#           格式：染色体<tab>起始位置<tab>终止位置<tab>平均信号值

set -euo pipefail

IN_FILE="$1"
CHROM_SIZES="$2"
OUT_DIR="$3"
BIN_SIZE="${4:-1000}"
STEP_SIZE="${5:-$BIN_SIZE}"

mkdir -p "$OUT_DIR"

#######################################
# 1. 判断文件类型 & 处理后缀
#######################################
FILENAME=$(basename "$IN_FILE")

case "$FILENAME" in
    *.bw)
        FILE_TYPE="bigwig"
        BASENAME="${FILENAME%.bw}"
        ;;
    *.bigwig)
        FILE_TYPE="bigwig"
        BASENAME="${FILENAME%.bigwig}"
        ;;
    *.bigWig)
        FILE_TYPE="bigwig"
        BASENAME="${FILENAME%.bigWig}"
        ;;
    *.bedGraph)
        FILE_TYPE="bedgraph"
        BASENAME="${FILENAME%.bedGraph}"
        ;;
    *)
        echo "[ERROR] Unsupported file type: $FILENAME" >&2
        exit 1
        ;;
esac

OUT_BED="${OUT_DIR}/${BASENAME}.bin${BIN_SIZE}.step${STEP_SIZE}.bp"

TMP_BINS=$(mktemp)

#######################################
# 2. 生成 genome bins
#######################################
echo "[INFO] Generating genome bins (${BIN_SIZE} bp)..."
bedtools makewindows \
    -g "$CHROM_SIZES" \
    -w "$BIN_SIZE" \
    -s "$STEP_SIZE" \
| awk 'BEGIN{OFS="\t"}{
    print $1, $2, $3, $1":"$2"-"$3
}' > "$TMP_BINS"

#######################################
# 3. 根据文件类型计算 bin 平均值
#######################################
if [[ "$FILE_TYPE" == "bigwig" ]]; then
    echo "[INFO] Detected BigWig format"
    echo "$IN_FILE"
    bigWigAverageOverBed \
        "$IN_FILE" \
        "$TMP_BINS" \
        "${OUT_BED}.tmp"

    echo "CREATED"
    # name size covered sum mean0 mean
    awk 'BEGIN{OFS="\t"}{
        split($1,a,":|-");
        print a[1], a[2], a[3], $6
    }' "${OUT_BED}.tmp" > "$OUT_BED"

    rm "${OUT_BED}.tmp"
    rm "$TMP_BINS"

elif [[ "$FILE_TYPE" == "bedgraph" ]]; then
    echo "[INFO] Detected bedGraph format"

    #######################################
    # 1. bedGraph → 1bp expansion
    #######################################
    awk 'BEGIN{OFS="\t"}
    {
        len=$3-$2
        if(len==1){
            print
        }else{
            for(i=0;i<len;i++){
                print $1,$2+i,$2+i+1,$4
            }
        }
    }' "$IN_FILE" > "$IN_FILE".expanded

    #######################################
    # 2. 排序
    #######################################
    sort -k1,1 -k2,2n "$IN_FILE".expanded > "$IN_FILE".sorted
    sort -k1,1 -k2,2n "$TMP_BINS" > "$TMP_BINS".sorted

    #######################################
    # 3. bin 内平均
    #######################################
    bedtools map \
        -a "$TMP_BINS".sorted \
        -b "$IN_FILE".sorted \
        -c 4 \
        -o mean \
        -null 0 \
    | awk 'BEGIN{OFS="\t"}{
        split($4,a,":|-");
        print a[1], a[2], a[3], $5
    }' > "$OUT_BED"

    #######################################
    # 4. 清理
    #######################################
    rm "$TMP_BINS"
    rm "$IN_FILE".expanded
    rm "$IN_FILE".sorted
    rm "$TMP_BINS".sorted

fi

echo "[DONE] Output written to $OUT_BED"

