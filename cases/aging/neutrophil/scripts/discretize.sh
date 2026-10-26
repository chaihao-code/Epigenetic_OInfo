#!/usr/bin/env bash
# Step 4: divInput 信号离散化（n_bins=5，仅常染色体，与 thyroid 约定一致）
export EPINFO_ROOT="${EPINFO_ROOT:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && while [ ! -f .epigenetic_oinfo_root ]; do cd .. || exit 1; done && pwd)}"
set -uo pipefail

PY="python3"
DISCRETIZE="${EPINFO_ROOT}/common/01_preprocessing/discretize_signal.py"
DATA_DIR="${EPINFO_ROOT}/data/aging/neutrophil"
IN_DIR="${DATA_DIR}/input_normed_binned_results"
OUT_DIR="${DATA_DIR}/input_normed_binned_discreted"
mkdir -p "$OUT_DIR"

disc_one() {
    local f="$1" outdir="$2"
    local base out
    base=$(basename "$f" .divInput.bp)
    out="${outdir}/${base}.divInput.discrete.bed"
    [[ -s "$out" ]] && { echo "[SKIP] $base"; return 0; }
    if python3 \
        ${EPINFO_ROOT}/common/01_preprocessing/discretize_signal.py \
        -i "$f" -o "$out" -n 5 --keep-autosomes > /dev/null 2>&1 && [[ -s "$out" ]]; then
        echo "[DONE] $base"
    else
        echo "[FAIL] $base"; return 1
    fi
    return 0
}
export -f disc_one

ls "${IN_DIR}"/*.divInput.bp | xargs -P 16 -I{} bash -c 'disc_one "$@"' _ {} "${OUT_DIR}"

n_done=$(ls "${OUT_DIR}" | grep -c "discrete.bed$")
echo "[SUMMARY] discrete: ${n_done} files (expect 72)"
[[ "$n_done" -ge 72 ]]
