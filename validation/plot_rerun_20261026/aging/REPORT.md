# 衰老案例作图重跑

14个历史入口，当前保留成功4个，按用户要求移除10个；当前8幅、8个图文件。

甲状腺去卷积仅保留figQ1_composition.png；其他去卷积图片已按要求删除，默认重跑不会再生成。甲状腺GSEA和其他项目图片保留。

本次所有数值输入均为release内已有中间表和元数据，未从原工作目录恢复数据；没有重算O-information、模拟混合或重新去卷积。混合入口通过run_aging_plots.py调用原绘图函数或执行原main中的绘图AST段，详细命令与实际读取文件哈希见report.json。

| 程序 | 内容 | 幅/文件 | 输入核实与本次方式 | 输出 |
|---|---|---:|---|---|
| [04_select_markers.py](../../../cases/aging/thyroid/deconvolution/scripts/04_select_markers.py) | 参考细胞甲基化marker热图：Panel A/B | 0/0 | 原入口从参考beta/npz重选marker；本次从reference_matrix_A/B中间表调用heatmap() | [removed_by_request](logs) |
| [05_validate_deconvolution.py](../../../cases/aging/thyroid/deconvolution/scripts/05_validate_deconvolution.py) | 模拟混合比例与估计免疫比例验证图（4面板） | 0/0 | 原入口重新模拟/NNLS；本次读取simulation_results.tsv执行原绘图段 | [removed_by_request](logs) |
| [06_deconvolve_samples.py](../../../cases/aging/thyroid/deconvolution/scripts/06_deconvolve_samples.py) | 6供体细胞组成与模型敏感性图（4面板） | 0/0 | 原入口读取甲基化bigWig并重新去卷积/bootstrap；本次读取细胞比例表调用figures() | [removed_by_request](logs) |
| [06_immune_age_association.py](../../../cases/aging/thyroid/deconvolution/scripts/06_immune_age_association.py) | 免疫细胞比例随年龄变化及回归（4面板） | 0/0 | 中间细胞比例表+样本年龄；本次执行轻量回归及绘图 | [removed_by_request](logs) |
| [07_oinfo_age_adjustment.py](../../../cases/aging/thyroid/deconvolution/scripts/07_oinfo_age_adjustment.py) | 免疫调整前后年龄斜率、置信区间及稳健性（6面板） | 0/0 | 原入口从O-info及比例表重新回归/bootstrap；本次读取oinfo_age_adjustment.tsv及散点源表执行绘图段 | [removed_by_request](logs) |
| [08_region_filtering.py](../../../cases/aging/thyroid/deconvolution/scripts/08_region_filtering.py) | 排除免疫特异区域与matched downsampling敏感性（3面板） | 0/0 | 原入口读离散轨迹重新计算O-information；本次读取filtering_values/region_filtering中间表执行绘图段 | [removed_by_request](logs) |
| [09_marker_go_analysis.py](../../../cases/aging/thyroid/deconvolution/scripts/09_marker_go_analysis.py) | 免疫/甲状腺marker表达、年龄关联和免疫区间重叠（3面板） | 0/0 | 原入口读RNA bigWig、GTF、离散轨迹/注释并做GO分析；本次只读immune_marker_analysis/immune_region_enrichment表执行绘图段 | [removed_by_request](logs) |
| [10_histone_association.py](../../../cases/aging/thyroid/deconvolution/scripts/10_histone_association.py) | 6种组蛋白信号的年龄斜率；信号与免疫比例关系（各6面板） | 0/0 | 原入口优先bin_means_cache缓存，缺缓存读取组蛋白bigWig；本次只读histone_immune_association及细胞比例表执行绘图段 | [removed_by_request](logs) |
| [12_fig_q1_composition.py](../../../cases/aging/thyroid/deconvolution/scripts/12_fig_q1_composition.py) | 细胞组成及Panel/Model稳健性（3面板） | 1/1 | 中间细胞比例表+样本元数据；纯绘图 | [success](../../../cases/aging/thyroid/deconvolution/results/figures) |
| [15_fig_volcano.py](../../../cases/aging/thyroid/deconvolution/scripts/15_fig_volcano.py) | 免疫调整后的年龄斜率/P值火山图 | 0/0 | 中间oinfo_age_fdr_table.tsv；纯绘图 | [removed_by_request](logs) |
| [16_fig_lollipop.py](../../../cases/aging/thyroid/deconvolution/scripts/16_fig_lollipop.py) | 年龄斜率调整前后相对差的棒棒糖图 | 0/0 | 中间oinfo_age_fdr_table.tsv；纯绘图 | [removed_by_request](logs) |
| [plot_oinfo_vs_age.py](../../../cases/aging/neutrophil/scripts/plot_oinfo_vs_age.py) | 40区域及CpG区域O-info与年龄散点：全部/仅成人，共4图 | 4/4 | 中间O-information CSV+供体元数据；拟合图中趋势线 | [success](../../../cases/aging/neutrophil/plots) |
| [plot_paper_figures.py](../../../cases/aging/neutrophil/scripts/plot_paper_figures.py) | O-info供体×区域热图；年龄斜率最大的6区域散点 | 2/2 | 中间O-information CSV+供体元数据；重新选择斜率最大的区域 | [success](../../../cases/aging/neutrophil/plots) |
| [plot_volcano_oinfo_age.py](../../../cases/aging/neutrophil/scripts/plot_volcano_oinfo_age.py) | 成人中性粒细胞O-info年龄斜率/P值火山图 | 1/1 | 中间oinfo_age_regression_adults.tsv；纯绘图 | [success](../../../cases/aging/neutrophil/plots) |

## 运行与核验

```bash
python validation/plot_rerun_20261026/aging/run_aging_plots.py
python validation/plot_rerun_20261026/aging/make_report.py
```

PNG逐张解码，PDF检查文件头和EOF，全部图记录SHA256。原程序05的模拟验证图原来固定标注sigma=0.03，但实际选sigma=0.05数据；已改为按实际数据标注。Panel A/B heatmap保留中间表marker原始顺序。Arial字体已作为真实文件保存到data/fonts/arial.ttf。

原入口中06_deconvolve_samples、09_marker_go_analysis会读取原始/参考bigWig等文件，10_histone_association缺缓存时也会读取bigWig；08_region_filtering从离散轨迹重算O-info。不能把这些原入口当成单纯作图命令。

保留的图片位于各案例结果目录，路径见上表；日志在logs/。[机器清单](report.json)，[缩略图](contact_sheet.jpg)。
