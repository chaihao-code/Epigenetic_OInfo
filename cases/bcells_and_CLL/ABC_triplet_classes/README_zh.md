# B 细胞 ABC 协同 triplet 热图

主图在 plots/A_sample_heatmap、B_sample_heatmap、C_sample_heatmap（PNG/PDF）。
纵轴依次为 8 个 R、5 个 S、4 个 pB 样本，每行一个样本；R/S 仅为标签，作为平行 CLL 样本处理，未做 R/S 差异检验。
横轴为 8 种功能区域 × 5 种重复序列/CpG 注释，沿用参考图的顺序、配色及单元格计数。
黑色横线仅分隔 13 个 CLL 样本与 4 个 pB 样本。

主图与参考图一致，统计 Ω<0 且下尾置换 p_down<0.05 的 triplet 数量；每项有 1000 次置换。
后缀 _FDR005 的图改为 Ω<0 且 q_down_global_BH<0.05；其 FDR 沿用全部 23800 项检验，主图没有把名义 P 值写成 FDR。
_pooled_average_heatmap 分别展示 R（n=8）、S（n=5）、pB（n=4）三组均值，各组内部按样本等权平均。
R 和 S 分开汇总，不合并为一个 CLL 均值；这里仅做描述性绘图，不进行 R/S 组间显著性比较。
合并图 plots/ABC_pooled_average_heatmap 已更新为 Ω<0 且全局 BH FDR<0.05，
覆盖全部 23800 项检验，在按符号、类别、组别筛选前完成校正。
其绘图数值由原始检验表重新汇总，ABC_pooled_average_counts.tsv 的 mean_count 为 FDR 后计数均值，
mean_count_nominal 保留名义 P 值均值供对照。该图单独运行 ../scripts/plot_abc_pooled_combined.py 重画；
该脚本独立复核 BH 值及样本级计数，保留原有 18 cm 宽度、正方形格子和配色。
合并图行标签仅显示 R、S、pB，不显示 n；格内数字保持 3.5 pt，其余文字统一为 8 pt。
三个色条共用一个文字标签，画布为 18 × 8.035 cm。

ABC 分类沿用参考定义，令 T={H3K9me3,H3K27me3,H3K36me3}：
- A：含 DNA 甲基化，且至少包含 T 中一个标记。
- B：不含 DNA 甲基化，且至少包含 T 中一个标记。
- C：不含 T 中任何标记。

本数据只有 7 个信号，ABC 的组合上限分别为 12、19、4，总计 35；
参考图使用 8 个信号，上限为 15、31、10，不能直接套用。
各类原始计数受组合上限影响，counts_by_sample_region_class.tsv 同时提供各类比例。

pB 对应 CEMT0129、CEMT0134、CEMT0135、CEMT0138；本地原始 metadata 注释为 B cells、disease=None，来源为 Tonsil。
sample_labels.tsv 保留完整样本映射。所有单项检验、分类和计数均已导出，可从现有数据重画，无需重新提交计算任务。
复现脚本：../scripts/plot_abc_triplet_classes.py。
