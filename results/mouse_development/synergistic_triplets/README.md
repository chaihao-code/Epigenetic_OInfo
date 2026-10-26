# 10,000 次置换的 BH-FDR 结果

82 个样本全部完整，共 3,280 个 sample × block 校正家族、114,800 个三阶项。每个家族只对全部 35 个三阶项的 `p_down` 做 BH，不预筛显著项，不合并其他尾部检验。显著负向项要求 `observed_hoi < 0` 且 `q_down_sample_block_bh <= 0.05`。

校正前显著负向项 24,769 个，BH35 后保留 **20,505** 个。

## 文件

- `all_triplets_fdr_10000.tsv.gz`：全部三阶项、原始 p 值、BH q 值及显著性标记。
- `significant_negative_triplets_10000.tsv`：BH35 后保留的显著负向项。
- `summary_by_sample.tsv`、`summary_by_stage.tsv`、`summary_by_block.tsv`、`summary_by_triplet.tsv`、`summary_by_triplet_class.tsv`、`summary_by_sample_block.tsv`：各层面数量汇总。
- `analysis_metadata.json`：分析口径及一万次 FDR 汇总。
- `input_snapshot.tsv`：输入文件及 SHA256；`validation.json`：全量 BH 独立公式核对结果。

早期样本中的 K27ac、K4M3、K4me1 分别统一为 H3K27ac、H3K4me3、H3K4me1；WGBS 统一为 CpG_methyl。`subset_names` 保留原始名称。匹配键为 sample × block × 规范化信号组合，不依赖 subset_idx。

## 重算

在项目根目录运行：

```bash
/usr/bin/python3 project/mouse_development/script/recalculate_triplet_fdr.py
```

该脚本只读取已有置换汇总表，不重新运行置换计算。

## 一万次结果的绘图

下面两个脚本默认读取本目录的 `all_triplets_fdr_10000.tsv.gz`，结果直接保存到本目录：

```bash
OPENBLAS_NUM_THREADS=1 python project/mouse_development/script/plot_special_triplet_f1.py
OPENBLAS_NUM_THREADS=1 python project/mouse_development/script/plot_triplet_counts_corrected_by_stage.py
```

- `H3K9me3_H3K27me3_X_F1_heatmap.png/.pdf`：特殊组合在 Triplet B 内的 F1 热图；X 为 ATAC、H3K27ac、H3K4me1、H3K4me3。各格是时期内各组织 F1 的算术均值。
- `f1_by_sample_block.tsv`、`f1_by_stage_block.tsv`：F1 的样本级混淆计数与绘图数值。
- `triplet_ABC_stacked.png/.pdf`：A/B/C 三类显著负向项数量的合并热图，各格为时期内组织均值。
- `triplet_A/B/C_corrected_synergy_count_by_stage.png/.pdf`：三类单独热图（文件名中 A/B/C 分别替换）。
- `total_corrected_synergy_count_by_stage.png/.pdf`、`total_synergistic_triplet_count_by_stage.png/.pdf`：总数量热图和发育时期折线图。
- `stage_coverage.tsv`、`triplet_counts_plot_metadata.json`、`special_triplet_f1_metadata.json`：时期样本覆盖和绘图口径。

合并热图 `triplet_ABC_stacked.png/.pdf` 的画布为 18 × 19.5 cm，格内数字保持 6 pt，
其余文字统一为 Arial 8 pt。`triplet_ABC_stacked_layout.json` 保存尺寸、绘图矩阵、
色轴、来源哈希及文字边界检查。仅使用已有时期汇总表重绘该图：

```bash
OPENBLAS_NUM_THREADS=1 python project/mouse_development/script/plot_triplet_counts_corrected_by_stage.py --stacked-only
```

按用户要求，一千次、十万次结果目录及一万次与十万次的比较文件均已删除。保留一万次原始结果、FDR 结果及绘图文件。`deleted_runs.json` 记录删除范围。重算脚本只处理一万次结果。

## 精选区域的七信号 ABC 叠加轨迹

`triplet_ABC_selected_region_trajectories_7signals.png/.pdf` 使用完整的
zygote–P0 时间轴，每个区域叠加 A/B/C 三条曲线，排成两行三列。
上排依次为 exon/SINE、distal enhancer CTCF−/SINE、proximal enhancer
CTCF+/CpG island；下排为 intron/SINE、distal enhancer CTCF−/LINE、None/None。
Exon/intron 的 SINE 区域提供 gene body 对照；中列的 distal enhancer
CTCF− 的 SINE/LINE 提供相同功能区域内的重复序列注释对照。
None/None 对应原始 `None_AND_None` 区块，在七信号植入前阶段仍有
A/B 检出，原肠胚及晚期没有显著 ABC；八信号晚期配套图中亦为零。
注释区块为交集；SINE 区块不表示已排除 CpG island。
画布宽 18 cm、高约 9.67 cm，全部文字为 Arial 8 pt；整体和每个小图的
高度均压缩至上一版的 3/4，小图高/宽约 0.5606。所有区域标题使用单行，
两排共用一个纵轴标题。上下排共享横轴，完整的
17 个时间标签仅在下排显示，全部为单行竖排。早期阶段使用分类坐标，
晚期维持实际日距（P0 按 19.5 天），每一天与早期一个分类间隔同宽。

只纳入 ATAC、DNA methylation、H3K9me3、H3K27me3、H3K4me1、H3K4me3、
H3K27ac 七信号，不包含 H3K36me3。A/B/C 分母分别为 9/16/10。
保留原七信号图的全部 82 个样本：各时期跨全部可用组织平均，
植入前每期 n=1，原肠胚阶段 n=2–4，E11.5–P0 每期 n=8–12。
误差线沿用总体 SD / sqrt(n)，反映样本间差异。筛选为
O-information < 0 且 BH35 q <= 0.05（10,000 次置换）。

```bash
python project/mouse_development/script/plot_triplet_ABC_selected_regions.py
```

同前缀的 `_by_sample.tsv`、`_summary.tsv`、`_metadata.json` 保存样本计数、
绘图数值和来源/校验记录。脚本验证七信号全集、重新核对 BH35，
并逐项比对已有全区域计数表，同时检查文字边界和时间标签重叠。

按用户要求，旧版 3 × 4 ABC 轨迹图及配套汇总、元数据、校验文件和日志，
以及选区探索生成的临时比较表，已清理；只保留上述最终轨迹图及复现文件。
热图、全区域计数、原始检验结果和源代码保留。两项目的删除清单见
`../../H3K36me3_developmental_switch/synergistic_triplets/data/trajectory_cleanup_manifest.json`。

## ABC 与特殊组合 F1 的 a/b 合并图

`triplet_ABC_and_F1_panels_ab.png/.pdf` 宽 18 cm、高 22.5 cm。
面板 a 包含 Triplet A/B/C 三张热图；面板 b 展示
H3K9me3 + H3K27me3 + X 的 F1，X 为 ATAC、H3K27ac、H3K4me1、H3K4me3。
四张热图按相同的 17 个时期、40 个区域对齐，共用最底部的区域标签。
A/B 标题明确列出 `H3K9me3 or H3K27me3`。
所有文字统一为 Arial 8 pt；四张热图均不显示格内数字，时期标签不显示 `n=...`。
区域标签使用 `Promoter`、`Proximal\nenh. CTCF+` 等统一格式。
面板 a 与 b 之间的矩阵间距为 1 cm，a 内部矩阵间距为 0.32 cm。
顶部方法说明移至本文：10,000 次置换，sample × region 内 BH35 q≤0.05，
共 82 个样本，按时期内样本/组织求均值。
计数、F1 数值和各自色轴沿用现有汇总表；F1 保持每样本计算后在时期内求均值的口径。
原始数据表及两张独立图均保留。

```bash
OPENBLAS_NUM_THREADS=1 python project/mouse_development/script/plot_triplet_ABC_and_F1_panels.py
```

同前缀的 `_layout.json` 保存来源哈希、绘图矩阵、尺寸、字号及区域对齐和文字边界检查。

## 七信号在显著协同 triplet 中的出现次数：完整发育过程

`synergistic_triplet_mark_trajectories_7signals.pdf/.png`：参考八信号项目的
单信号参与次数轨迹图，尺寸 8.5 × 5 cm，全部文字为 Arial 8 pt。
17 个时期标签竖排，图例四列两行，纵轴标题分行显示。
使用相同的七种信号颜色，不包含 H3K36me3。

每个样本内，对全部 40 类区域的显著协同 triplet 计数，
某个信号每参与一个显著 triplet 就累计一次；相同组合在不同区域分别计数。
沿用 `observed_hoi < 0` 且 sample × block 内 BH35 q≤0.05、10,000 次置换的结果。
不使用 A/B/C 各类组合数作分母，也不重新运行置换。

保留全部 82 个样本，覆盖 zygote–P0 的 17 个时期，包括 E6.5/E7.0/E7.5；
ICM、TE 分开，晚期按各时期全部可用组织取均值。
误差线参照所给单信号图，使用样本 SD（ddof=1）/sqrt(n) 的 SEM，单样本时期 SEM=0。
这与八信号参考图的 ICM/TE 合并、晚期固定五组织的取样不同，
也与本目录 ABC 轨迹图所用的总体 SD 计算口径不同。

同前缀 `_by_sample.tsv` 保存 82×7=574 条计数，`_summary.tsv` 保存 17×7=119 条均值和 SEM；
`_metadata.json` 保存来源哈希、阶段样本数、计数口径、BH35 重算核对及排版检查。
脚本还验证每个样本七信号计数之和等于其显著 triplet 总数的三倍，
并与现有 A+B+C 的全区域时期总数逐项比对。

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /lustre/home/2101110354/anaconda3/envs/hoi_env/bin/python3 project/mouse_development/script/plot_synergistic_triplet_mark_trajectories.py
```

仅调整绘图规格时可加 `--plot-only`，读取已保存的单信号计数汇总。
