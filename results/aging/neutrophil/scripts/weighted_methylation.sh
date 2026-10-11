#!/usr/bin/env bash
# Step 2: WGBS 甲基化 coverage 加权分箱
#   过滤 3 <= cov <= 200，空 bin（无合格 CpG）记 0
#   输出：加权甲基化 bp -> binned_results/ + input_normed_binned_results/
#         QC: CpG 计数 / coverage 总和 bp -> methylation_qc/
set -uo pipefail

H="/lustre/home/2101110354/anaconda3/envs/hoi_env/bin"
DATA_DIR="/lustre/home/2101110354/sequence_entropy/data/aging/neutrophil"
OUT_BINS="${DATA_DIR}/binned_results"
OUT_NORMED="${DATA_DIR}/input_normed_binned_results"
QC_DIR="${DATA_DIR}/methylation_qc"
CHROM_SIZES="/lustre/home/2101110354/sequence_entropy/data/genome_features/hg38/hg38.chrom.sizes"
MIN_COV=3
MAX_COV=200

mkdir -p "$OUT_BINS" "$OUT_NORMED" "$QC_DIR"

meth_one() {
    local meth_bw="$1"
    local base core cov_bw
    base=$(basename "$meth_bw" .bigWig)
    core="${base#[0-9]*.}"; core="${core%.WGB-Seq.methylation}"

    cov_bw=$(ls /lustre/home/2101110354/sequence_entropy/data/aging/neutrophil/*.WGB-Seq.coverage.bigWig 2>/dev/null | while read -r c; do
        cb=$(basename "$c" .bigWig)
        cc="${cb#[0-9]*.}"; cc="${cc%.WGB-Seq.coverage}"
        [[ "$cc" == "$core" ]] && { echo "$c"; break; }
    done)
    if [[ -z "$cov_bw" ]]; then echo "[FAIL] $base (no coverage pair)"; return 1; fi

    local meth_out="${OUT_BINS}/${base}.bin1000.step1000.bp"
    local cpg_out="${QC_DIR}/${base}.cpg_count.bin1000.step1000.bp"
    local cov_out="${QC_DIR}/${base}.cov_sum.bin1000.step1000.bp"
    if [[ -s "$meth_out" && -s "$cpg_out" && -s "$cov_out" ]]; then
        echo "[SKIP] $base"; return 0
    fi

    local T; T=$(mktemp -d)
    local DATA_DIR="/lustre/home/2101110354/sequence_entropy/data/aging/neutrophil"
    local CHROM_SIZES="/lustre/home/2101110354/sequence_entropy/data/genome_features/hg38/hg38.chrom.sizes"
    local H="/lustre/home/2101110354/anaconda3/envs/hoi_env/bin"
    local MIN_COV=3 MAX_COV=200

    local rc=0
    "${H}/bigWigToBedGraph" "$meth_bw" "${T}/meth.bg" 2>/dev/null
    "${H}/bigWigToBedGraph" "$cov_bw" "${T}/cov.bg" 2>/dev/null

    # 坐标对齐：两文件的区间化方式不同（bigWig 存储合并了部分相邻同值区间），
    # 不能 paste。把 coverage 展开成 1bp，再用 bedtools map 取重叠的甲基化值。
    awk 'BEGIN{OFS="\t"} {len=$3-$2; if(len==1){print; next}
        for(i=0;i<len;i++) print $1,$2+i,$2+i+1,$4}' "${T}/cov.bg" \
        | sort -k1,1 -k2,2n -S 4G > "${T}/cov.1bp.sorted"
    sort -k1,1 -k2,2n -S 4G "${T}/meth.bg" > "${T}/meth.sorted"

    bedtools map -a "${T}/cov.1bp.sorted" -b "${T}/meth.sorted" \
        -c 4 -o mean -null -1 | \
    awk -v t1="$MIN_COV" -v t2="$MAX_COV" 'BEGIN{OFS="\t"}
        ($5 < 0) { next }                              # 无甲基化记录的位点
        ($4 >= t1 && $4 <= t2) {
            print $1,$2,$3,$4*$5 > "'"${T}"'/prod.bg"
            print $1,$2,$3,$4    > "'"${T}"'/covf.bg"
            print $1,$2,$3,1     > "'"${T}"'/cnt.bg"
        }'

    # 1kb 窗口（与 compute_bin_average.sh 同一 grid）
    bedtools makewindows -g "${CHROM_SIZES}" -w 1000 -s 1000 \
        | awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' | sort -k1,1 -k2,2n -S 4G > "${T}/windows.bed"

    for pair in "prod.bg sum" "covf.bg sum" "cnt.bg count"; do
        set -- $pair
        sort -k1,1 -k2,2n -S 4G "${T}/$1" > "${T}/$1.sorted"
        bedtools map -a "${T}/windows.bed" -b "${T}/$1.sorted" -c 4 -o "$2" -null 0 \
            > "${T}/$1.mapped"
    done

    # 加权甲基化 = sum(meth*cov)/sum(cov)，无合格 CpG 的 bin 记 0
    paste "${T}/prod.bg.mapped" "${T}/covf.bg.mapped" "${T}/cnt.bg.mapped" | \
    awk 'BEGIN{OFS="\t"}
        { w = ($8 > 0) ? $4/$8 : 0;
          print $1,$2,$3,w  > "'"${T}"'/meth.weighted"
          print $1,$2,$3,$12 > "'"${T}"'/cpg.count"
          print $1,$2,$3,$8  > "'"${T}"'/cov.sum" }' || rc=1

    if [[ $rc -eq 0 ]]; then
        cp "${T}/meth.weighted" "$meth_out"
        cp "${T}/cpg.count" "$cpg_out"
        cp "${T}/cov.sum" "$cov_out"
        cp "${T}/meth.weighted" "${OUT_NORMED}/$(basename "$meth_out")"
        echo "[DONE] $base"
    else
        echo "[FAIL] $base"
    fi
    rm -rf "$T"
    return $rc
}
export -f meth_one
export OUT_BINS OUT_NORMED QC_DIR

ls "${DATA_DIR}"/*.WGB-Seq.methylation.bigWig \
    | xargs -P 12 -I{} bash -c 'meth_one "$@"' _ {}

n_done=$(ls "${OUT_BINS}" | grep -c "WGB-Seq.methylation.bin1000.step1000.bp$")
echo "[SUMMARY] weighted methylation: ${n_done} files (expect 12)"
[[ "$n_done" -ge 12 ]]
