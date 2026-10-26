# B-cell/CLL 作图重跑记录

11个入口全部成功，新生成19幅图（19 PNG + 19 PDF）。[全部新图缩略预览](contact_sheet.jpg)。PNG逐个解码验证；PDF检查文件头/EOF。所有输出位于当前release目录。

原release缺少多数计算中间表。本次从原 sibling workspace 复制65个CSV/TSV/GMT/RNK文件（448,053,695 bytes），没有复制历史图片。[输入恢复清单](restored_inputs.json)包含原路径、目标路径和SHA256；[外部输入链接清单](external_input_links.json)标明仍经现有软链接读取的beta矩阵和GTF。

| 程序 | 新生成内容 | 幅数/文件数 | 输入性质 | 输出目录 |
|---|---|---:|---|---|
| `plot_dmr_volcano.py` | 13名CLL供体1 kb差异甲基化火山图 | 1/2 | 中间结果；纯绘图，全部2,697,259行 | [dmr_repro_cemt/output](../../../cases/bcells_and_CLL/dmr_repro_cemt/output) |
| `plot_igh_tracks_corrected.py` | IGH区域甲基化轨迹、组差异和coverage条带 | 1/2 | 已处理beta矩阵+GENCODE原始参考注释；不是只读分析结果表，现场计算组均值 | [dmr_repro_cemt/output](../../../cases/bcells_and_CLL/dmr_repro_cemt/output) |
| `plot_de_male_ighv_volcano_clear.py` | 男性RNA method/IGHV校正差异表达火山图 | 1/2 | 中间结果；纯绘图 | [results_rna_de_male_ighv](../../../cases/bcells_and_CLL/results_rna_de_male_ighv) |
| `plot_gsea_eight_panel_v2_18x10cm.py` | 8条代表GO通路的GSEA富集曲线组合图 | 1/2 | GSEA中间结果+冻结基因集；重建running ES曲线，无置换 | [results_gsea_rna_male/male_RNA_method_IGHV_adjusted/eight_panel_GSEA_v2_18x10cm](../../../cases/bcells_and_CLL/results_gsea_rna_male/male_RNA_method_IGHV_adjusted/eight_panel_GSEA_v2_18x10cm) |
| `plot_gsea_representative_bars.py` | 代表GO通路NES/FDR条形图 | 1/2 | 中间结果；纯绘图 | [results_gsea_rna_male/male_RNA_method_IGHV_adjusted/representative_pathways](../../../cases/bcells_and_CLL/results_gsea_rna_male/male_RNA_method_IGHV_adjusted/representative_pathways) |
| `plot_gsea_top20_bars.py` | Top20 GO通路NES/FDR条形图 | 1/2 | 中间结果；纯绘图 | [results_gsea_rna_male/male_RNA_method_IGHV_adjusted/top20_pathways](../../../cases/bcells_and_CLL/results_gsea_rna_male/male_RNA_method_IGHV_adjusted/top20_pathways) |
| `validate_male_de_in_females.py` | 38个男性DE基因在4名女性CLL-R中的表达热图与投影图 | 1/2 | 已导入表达矩阵+中间DE表；混合分析/绘图，重新归一化并拟合男性投影模型 | [results_rna_de_female_validation](../../../cases/bcells_and_CLL/results_rna_de_female_validation) |
| `plot_abc_triplet_classes.py` | ABC三元组：A/B/C每样本nominal与FDR热图、各组均值热图，共9幅 | 9/18 | 中间结果+样本元数据；重新聚合计数 | [ABC_triplet_classes/plots](../../../cases/bcells_and_CLL/ABC_triplet_classes/plots) |
| `plot_abc_pooled_combined.py` | ABC三类的R/S/pB组均值组合热图 | 1/2 | 中间结果；独立验证23,800-test BH后聚合 | [ABC_triplet_classes/plots](../../../cases/bcells_and_CLL/ABC_triplet_classes/plots) |
| `plot_abc_pooled_compact_preview.py` | ABC组均值紧凑版组合热图 | 1/2 | 中间结果；纯绘图 | [ABC_triplet_classes/plots](../../../cases/bcells_and_CLL/ABC_triplet_classes/plots) |
| `plot_global_oinfo_and_ABC_triplets.py` | 全系统O-information组间差异与ABC三元组组合图 | 1/2 | 中间结果；重新做40区域Welch检验 | [global_oinfo_and_triplets/plots](../../../cases/bcells_and_CLL/global_oinfo_and_triplets/plots) |

## 输入核实

纯绘图程序没有读取FASTQ、原始bigWig或重新计算O-information。大部分输入是已有O-information、triplet显著性、DE、GSEA、DMR结果。以下入口需要明确区分：

- `plot_igh_tracks_corrected.py`读取全基因组prepared beta矩阵及GENCODE GTF参考注释，计算局部组均值和差异；GTF不属于分析中间文件。两者是现有外部软链接。
- `validate_male_de_in_females.py`从已导入counts/TPM重做归一化和男性模型拟合再画图，属于分析与作图混合入口。
- `plot_global_oinfo_and_ABC_triplets.py`从已有order-7数值重做40区域Welch检验；`plot_abc_pooled_combined.py`重新验证全局BH，均不运行基因组计算。
- GSEA 8面板从已有排名和GMT重建running ES曲线，并断言与已有ES一致，不运行GSEA置换。

## 运行证据与修复

最终逐项命令、退出码、时长、结果路径与SHA256在[report.json](report.json)；[run.py](run.py)可在release中复跑，先并行8项，再按依赖顺序画ABC组合、紧凑版和global组合图。日志在[logs](logs)。

实际执行发现并修复：ABC两个入口缺少`os`导入；ABC pooled布局元数据引用未定义`df`；IGH入口未自动创建输出目录。科学定义未修改。

初始探测与输入恢复有时间重叠：8项初始失败（日志保留），其余3项探测时已可见恢复后的文件。正式重跑是在恢复后完成，并对失败修复后再次执行。

历史triplet和GSEA输入中的1,000次置换仍保留原义；本次只重画和必要表层统计，没有将其标成新的10,000次计算。
