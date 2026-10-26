#!/usr/bin/env bash
# Step 5: 生成 project/aging/neutrophil/data_info/<IHECRE ID>.txt
# 每个 donor 7 行：6 个 divInput 组蛋白 + 1 个加权甲基化 bp（与 thyroid data_info 格式对齐）
export EPINFO_ROOT="${EPINFO_ROOT:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && while [ ! -f .epigenetic_oinfo_root ]; do cd .. || exit 1; done && pwd)}"
set -uo pipefail

DATA_DIR="${EPINFO_ROOT}/data/aging/neutrophil"
INFO_DIR="${EPINFO_ROOT}/cases/aging/neutrophil/data_info"
NORMED="${DATA_DIR}/input_normed_binned_results"
mkdir -p "$INFO_DIR"

MARKS=(H3K4me1 H3K4me3 H3K9me3 H3K27me3 H3K36me3 H3K27ac)

for input_bp in "${DATA_DIR}"/binned_results/*.Input.ctl_raw.bin1000.step1000.bp; do
    ibase=$(basename "$input_bp" .Input.ctl_raw.bin1000.step1000.bp)
    core="${ibase#[0-9]*.}"                       # Blueprint.IHECRE00000004.3
    donor=$(echo "$core" | grep -oE "IHECRE[0-9]+")
    [[ -z "$donor" ]] && { echo "[FAIL] cannot parse donor from $ibase"; exit 1; }

    txt="${INFO_DIR}/${donor}.txt"
    : > "$txt"
    n=0
    for mark in "${MARKS[@]}"; do
        f=$(ls "${NORMED}"/*.${mark}.pval_signal.bin1000.step1000.divInput.bp 2>/dev/null | grep "$core")
        [[ -z "$f" ]] && { echo "[FAIL] ${donor} missing ${mark}"; exit 1; }
        echo "$f" >> "$txt"; n=$((n+1))
    done
    mf=$(ls "${NORMED}"/*.WGB-Seq.methylation.bin1000.step1000.bp 2>/dev/null | grep "$core")
    [[ -z "$mf" ]] && { echo "[FAIL] ${donor} missing methylation"; exit 1; }
    echo "$mf" >> "$txt"; n=$((n+1))
    echo "[DONE] ${donor}.txt (${n} lines)"
done

n_txt=$(ls "${INFO_DIR}"/*.txt | wc -l)
echo "[SUMMARY] data_info: ${n_txt} files (expect 12)"
[[ "$n_txt" -ge 12 ]]
