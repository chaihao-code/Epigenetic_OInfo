# 互信息分析

- `pairwise_mi_n3/`、`pairwise_mi_n5/`、`pairwise_mi_n10/`：各阶段的自身熵、两两 MI 和联合 MI 表。
- `results_pairwise_correlation_n_bins_5/`：原始相关性结果和历史任务文件。
- `scripts/`：分区 MI 变化曲线、联合/两两 MI 比值热图、比值分布及合成图绘制代码。
- `plots/`：已有图、绘图数据和分箱敏感性结果。
- `logs/`：原绘图任务日志；保留历史路径用于追溯。

绘图示例：从本目录运行 `python scripts/plot_joint_vs_pairwise_ratio_heatmap_v3.py`。
此操作会覆盖相应已有图片。MI 来源为 calculation_results 中的熵估计。
