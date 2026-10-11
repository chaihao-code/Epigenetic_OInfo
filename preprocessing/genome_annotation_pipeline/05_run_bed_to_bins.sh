#!/usr/bin/env bash
# 05 对每个特征生成 *_bins.tsv
# 用法: GENOME=hg38 bash 05_run_bed_to_bins.sh / GENOME=mm10 bash ...
set -euo pipefail
source "$(dirname "$0")/config.sh"
P="$(dirname "$0")/04_bed_to_bins.py"

run() { python3 "$P" "$ANN_DIR/$1" "$BINS" "$2" "$OUT_DIR"; }

run promoters.bed promoter
run exons.bed exon
run introns.bed intron
run cpg_island.bed CpG_island
run LINE.bed LINE
run SINE.bed SINE
run LTR.bed LTR

if [ "$GENOME" = "hg38" ]; then
    run proximal_enhancer.ctcf_bound.bed proximal_ctcf_bound_enhancer
    run proximal_enhancer.not_ctcf_bound.bed proximal_not_ctcf_bound_enhancer
    run distal_enhancer.ctcf_bound.bed distal_ctcf_bound_enhancer
    run distal_enhancer.not_ctcf_bound.bed distal_not_ctcf_bound_enhancer
fi
echo "[done] all *_bins.tsv in $OUT_DIR"
