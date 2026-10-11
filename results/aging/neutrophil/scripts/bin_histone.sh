#!/usr/bin/env bash
# Step 1: 组蛋白 pval_signal + Input 分箱（bin1000.step1000，hg38 grid）
set -uo pipefail

DATA_DIR="/lustre/home/2101110354/sequence_entropy/data/aging/neutrophil"
OUT_DIR="${DATA_DIR}/binned_results"
CHROM_SIZES="/lustre/home/2101110354/sequence_entropy/data/genome_features/hg38/hg38.chrom.sizes"
BIN_SCRIPT="/lustre/home/2101110354/sequence_entropy/script/preprocessing/generate_bin_avg/compute_bin_average.sh"

mkdir -p "$OUT_DIR"

bin_one() {
    local f="$1" outdir="$2"
    local base out
    base=$(basename "$f" .bigWig)
    out="${outdir}/${base}.bin1000.step1000.bp"
    if [[ -s "$out" ]]; then echo "[SKIP] $base"; return 0; fi
    if bash /lustre/home/2101110354/sequence_entropy/script/preprocessing/generate_bin_avg/compute_bin_average.sh \
        "$f" /lustre/home/2101110354/sequence_entropy/data/genome_features/hg38/hg38.chrom.sizes "$outdir" 1000 1000 >/dev/null 2>&1 \
        && [[ -s "$out" ]]; then
        echo "[DONE] $base"
    else
        echo "[FAIL] $base"
        return 1
    fi
    return 0
}
export -f bin_one

ls "${DATA_DIR}"/*.bigWig | grep -v "WGB-Seq" \
    | xargs -P 16 -I{} bash -c 'bin_one "$@"' _ {} "${OUT_DIR}"

n_done=$(ls "${OUT_DIR}" | grep -c "bin1000.step1000.bp$")
echo "[SUMMARY] binned_results: ${n_done} files (expect 84)"
[[ "$n_done" -ge 84 ]]
