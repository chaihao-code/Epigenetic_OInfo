#!/usr/bin/env bash
# Step 4: divInput 信号离散化（n_bins=5，仅常染色体，与 thyroid 约定一致）
set -uo pipefail

PY="/lustre/home/2101110354/anaconda3/envs/hoi_env/bin/python"
DISCRETIZE="/lustre/home/2101110354/sequence_entropy/script/save_discrete_result/discretize_signal.py"
DATA_DIR="/lustre/home/2101110354/sequence_entropy/data/aging/neutrophil"
IN_DIR="${DATA_DIR}/input_normed_binned_results"
OUT_DIR="${DATA_DIR}/input_normed_binned_discreted"
mkdir -p "$OUT_DIR"

disc_one() {
    local f="$1" outdir="$2"
    local base out
    base=$(basename "$f" .divInput.bp)
    out="${outdir}/${base}.divInput.discrete.bed"
    [[ -s "$out" ]] && { echo "[SKIP] $base"; return 0; }
    if /lustre/home/2101110354/anaconda3/envs/hoi_env/bin/python \
        /lustre/home/2101110354/sequence_entropy/script/save_discrete_result/discretize_signal.py \
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
