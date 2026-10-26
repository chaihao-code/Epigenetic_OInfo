# 总 O-information 与 ABC triplet 合并图

成图：`plots/global_oinfo_and_ABC_triplets_18x8p5cm.pdf`，并提供同名 PNG。
画布为 18 × 8.5 cm；格内数字沿用 3.5 pt，其余文字为 8 pt。
上方棒棒糖图、下方连续排列的 Triplet A/B/C 共用最底部的一套完整区域标签。
图例与行标签使用 CLL-R、CLL-S；pB 标签保留。
40 列依次为 8 类功能区域，每类按 CpG island、LINE、SINE、LTR、None 排列。

上方读取 `../results_oinfo_n_bins_5/` 中 13 位 CLL 供者的原始七信号总 O-information，
即 `_oinfo.csv` 的 `order=7` 行。每项均核对对应系统置换汇总中的 `observed_hoi`。
R=8、S=5，分别对应参考脚本的 Group2、Group3。
分组与区域顺序保存在本目录 `reference_config.json`，提取自原两个参考脚本：Welch t 检验后，按照 R>S 的方向取 P 值；
未改为零分布中心化指标，也未改变参考脚本的名义 P<0.05 着色条件。
纵轴显示真实 `log10(p value)`，刻度为负数，轴方向反转，因此越显著的点越靠上。
图上省略检验方向说明，完整定义及双侧 P 值保存在数据表与验证记录中。

下方直接读取 `../ABC_triplet_classes/plots/ABC_pooled_average_counts.tsv`，
保持 Ω<0 且全局 BH FDR<0.05 的原有计数均值，BH 范围为全部 23,800 项检验。
R、S、pB 各自按供者等权平均；格内数字、YlOrRd 配色及 A/B/C 的 0–12、0–19、0–4 色轴保持原样。
色条标签显示为 `mean occurence`，对应上述计数均值。
原 ABC 计数表与单独成图均保留，脚本用 SHA256 核对其内容未变。

在 B 细胞项目根目录复现（先按主 README 恢复输入）：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python global_oinfo_and_triplets/plot_global_oinfo_and_ABC_triplets.py
```

`data/global_oinfo_order7_by_sample.tsv` 保存逐供者数值；
`data/global_oinfo_order7_R_vs_S_by_region.tsv` 保存 40 个区域的组均值、效应方向和 P 值；
`data/global_oinfo_and_ABC_triplets_18x8p5cm_validation.json` 记录来源哈希、区域对齐、尺寸、字号及文字边界检查。
