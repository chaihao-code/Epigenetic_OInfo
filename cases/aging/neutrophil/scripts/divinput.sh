#!/usr/bin/env bash
# Step 3: divInput 归一化（signal/(input+0.01)，与 thyroid/scripts/norm.sh 一致）
export EPINFO_ROOT="${EPINFO_ROOT:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && while [ ! -f .epigenetic_oinfo_root ]; do cd .. || exit 1; done && pwd)}"
set -uo pipefail

DATA_DIR="${EPINFO_ROOT}/data/aging/neutrophil"
BIN_DIR="${DATA_DIR}/binned_results"
OUT_DIR="${DATA_DIR}/input_normed_binned_results"
mkdir -p "$OUT_DIR"

fails=0
for input_bp in "${BIN_DIR}"/*.Input.ctl_raw.bin1000.step1000.bp; do
    ibase=$(basename "$input_bp" .Input.ctl_raw.bin1000.step1000.bp)
    core="${ibase#[0-9]*.}"   # e.g. Blueprint.IHECRE00000004.3
    for sig_bp in "${BIN_DIR}"/*.pval_signal.bin1000.step1000.bp; do
        sbase=$(basename "$sig_bp" .pval_signal.bin1000.step1000.bp)
        score="${sbase#[0-9]*.}"; score="${score%.*}"   # 去掉 mark 后缀，得到 donor core
        [[ "$score" != "$core" ]] && continue
        out="${OUT_DIR}/$(basename "$sig_bp" .bp).divInput.bp"
        [[ -s "$out" ]] && continue
        paste "$sig_bp" "$input_bp" | \
        awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4/($8+0.01)}' > "$out" \
            || { echo "[FAIL] $sbase"; fails=$((fails+1)); }
        echo "[DONE] $(basename "$out")"
    done
done

n_done=$(ls "${OUT_DIR}" | grep -c "divInput.bp$")
echo "[SUMMARY] divInput: ${n_done} files (expect 72)"
[[ "$n_done" -ge 72 && "$fails" -eq 0 ]]
