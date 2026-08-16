#!/usr/bin/env bash
#SBATCH -J avg_bw
#SBATCH -c 1
#SBATCH -t 04:00:00
#SBATCH -o logs/%x_%j.out
#SBATCH -e logs/%x_%j.err

set -euo pipefail

BW_DIR=$1
OUT_DIR=$2
GROUP=$3
KMER=$4

SCRIPT_DIR="/lustre/home/2101110354/sequence_entropy/script/preprocessing/average_replicates_and_liftover"
CHAIN="/lustre/home/2101110354/sequence_entropy/script/preprocessing/average_replicates_and_liftover/mm9ToMm10.over.chain.gz"

TMP_DIR=${OUT_DIR}/tmp_${GROUP}_${KMER}
mkdir -p ${TMP_DIR} ${OUT_DIR} logs

echo "Processing ${GROUP}.${KMER}"

files=$(ls ${BW_DIR}/*_${GROUP}_single_cell_*.${KMER}.bw 2>/dev/null || true)

n=$(echo ${files} | wc -w)
if [ "${n}" -eq 0 ]; then
    echo "No files found, exit."
    exit 0
fi

echo "Found ${n} files"

# -----------------------------------------------
# 新流程：bigWig -> bedGraph -> liftOver -> 对 mapped bedGraph 求平均
# -----------------------------------------------
lifted_bg=${OUT_DIR}/${GROUP}.${KMER}.mm10.bedGraph
lifted_bgs=()

echo "Converting bigWig to bedGraph and liftOver mm9 -> mm10..."
for f in ${files}; do
    base=$(basename ${f} .bw)
    bg="${TMP_DIR}/${base}.bedGraph"
    mapped_bg="${TMP_DIR}/${base}.mm10.bedGraph"

    bigWigToBedGraph ${f} ${bg}

    # 拆分为 1 bp 区间
    split_bg="${TMP_DIR}/${base}.split1bp.bedGraph"
    awk -v OFS="\t" '{
        for (i = $2; i < $3; i++) {
            print $1, i, i+1, $4
        }
    }' ${bg} > ${split_bg}

    liftOver \
        ${split_bg} \
        ${CHAIN} \
        stdout \
        /dev/null \
    | grep -v '^lambda' \
    | awk -v OFS="\t" '($3 - $2) == 1' \
    | sort -k1,1 -k2,2n \
    | awk -v OFS="\t" '!seen[$1,$2,$3]++' \
    > ${mapped_bg}

    lifted_bgs+=("${mapped_bg}")
done

echo "Merging with unionbedg and averaging by actual coverage..."
union_bg=${TMP_DIR}/${GROUP}.${KMER}.union.bedGraph

echo "Merging with unionbedg..."
bedtools unionbedg -i "${lifted_bgs[@]}" > ${union_bg}

echo "Averaging by actual coverage..."
awk -v OFS="\t" '
{
    sum = 0;
    n = 0;
    for (i = 4; i <= NF; i++) {
        if ($i != "." && $i != "NA") {
            sum += $i;
            n++;
        }
    }
    if (n > 0) {
        print $1, $2, $3, sum / n;
    }
}' ${union_bg} > ${lifted_bg}

echo "Output written to ${lifted_bg}"
