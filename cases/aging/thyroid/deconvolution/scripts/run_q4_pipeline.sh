#!/usr/bin/env bash
# run_q4_pipeline.sh — Q4 结论的最小可复现流程
#   去卷积(Panel A/Model A 主分析) → 免疫比例 → 40 区间回归
#   → 显著性(nominal p + BH/BKY FDR) + 稳健性(模拟验证 / bootstrap / LOO / panel 敏感性)
#
# 输入前提(不在本脚本内生成):
#   1) config.yaml 与 data/sample_metadata/analysis_samples.tsv (6 例甲状腺样本)
#   2) 6 个样本的甲基化 bigWig (路径见 config.yaml)
#   3) data/processed/oinfo_long.tsv — 原始 O-information 结果
#      (40 区间 x 6 样本; 由原研究 HOI 分析经 scripts/prep_oinfo_table.py 整理)
#   4) 参考数据: 首次运行需先执行 scripts/02_download_reference.py
#      (GSE186458 的 27 个 .hg38.beta, 约 27 x 60MB)
#
# 用法: bash scripts/run_q4_pipeline.sh [stage ...]
#   不带参数 = 跑全部; 例: bash scripts/run_q4_pipeline.sh 4 6 7 只跑 04/06/07
# 日志: logs/q4_pipeline.log
set -euo pipefail
cd "$(dirname "$0")/.."
mkdir -p logs
LOG=logs/q4_pipeline.log

run() { echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOG"; "$@" 2>&1 | tee -a "$LOG"; }

STAGES="${@:-3 4 5 6 7 8 9}"
for s in $STAGES; do
  case $s in
    3) run python3 scripts/03_build_reference.py ;;   # 参考 bin 甲基化矩阵 (约 10-20 min)
    4) run python3 scripts/04_select_markers.py ;;    # marker panel + R 矩阵 (约 5 min)
    5) run python3 scripts/05_validate_deconvolution.py ;;  # 模拟验证 (约 3 min)
    6) run python3 scripts/06_deconvolve_samples.py ;;      # 6 样本去卷积 + bootstrap CI (约 5 min)
    7) run python3 scripts/06_immune_age_association.py ;;  # Q2: 免疫比例~年龄 (秒级)
    8) run python3 scripts/07_oinfo_age_adjustment.py ;;    # Q3/Q4: 40 区间回归 + LOO/bootstrap/panel 敏感性 (约 3 min)
    9) run python3 scripts/17_fdr_table.py ;;               # FDR 表 (秒级)
    *) echo "unknown stage: $s" >&2; exit 1 ;;
  esac
done

echo "[$(date +%H:%M:%S)] done. 关键输出:" | tee -a "$LOG"
echo "  results/tables/thyroid_cell_proportions.tsv   (免疫比例 + 95% CI)" | tee -a "$LOG"
echo "  results/tables/oinfo_age_adjustment.tsv       (β_age 调整前后/β_immune/LOO/bootstrap)" | tee -a "$LOG"
echo "  results/tables/oinfo_age_fdr_table.tsv        (nominal p + BH/BKY FDR)" | tee -a "$LOG"
