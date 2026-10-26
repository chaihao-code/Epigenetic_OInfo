# CEMT Group2 / Group3：1 kb 差异甲基化复跑

## 原始内容与方法

本目录原有一个约 455 MB 的输入 CSV、一个 R 程序和空的 `output/`。

- `input/all_bins_1kb_cemt_g2_g3.csv`：13 位供者的基因组窗口甲基化 beta 值（预期 0–1）；染色体末端窗口可能短于 1 kb。
- `scripts/dmr_limma_1kb_cemt_g2_g3.R`：用户提供的原始分析程序，本次不修改。
- Group2 为 CLL-R，8 人；Group3 为 CLL-S，5 人。delta_beta 定义为 Group2−Group3。
- 保留 chr1–chr22、至少 8 个样本不缺失的窗口。
- M-value = log2((beta + 1e-4)/(1−beta + 1e-4))。
- limma 模型为 `~ group + age + sex`，使用 `eBayes(trend=TRUE)` 和全部保留窗口的 BH 校正。
- 命中条件为 |delta_beta| > 0.2 且 q < 0.05。delta_beta 是两组原始 beta 均值之差；p/q 来自调整年龄和性别的 M-value 模型，两者不是同一个效应估计。

本次复跑沿用原程序细节：delta_beta 先四舍五入至四位小数再应用阈值，CEMT0097 的年龄按原脚本取 90（原临床表为 90+），不额外调整 IGHV 或 RNA 方法。全部 S 为男性，性别与分组部分混杂。结果是差异甲基化窗口，程序不合并相邻窗口形成连续 DMR。

## 环境与运行

在 B 细胞项目根目录运行；Python、R 与输入配置见主 README。历史运行使用 R 4.5.3，limma 3.66.0、data.table 1.18.6.1。原 limma 程序保持不变。

```bash
export BCELLS_PROJECT_DIR="$PWD"
export R_LIBS_USER="$PWD/dmr_repro_cemt/.Rlib"
Rscript dmr_repro_cemt/scripts/install_r_packages.R
python dmr_repro_cemt/scripts/check_input.py
Rscript dmr_repro_cemt/scripts/run_dmr.R
# 如需 Slurm：sbatch dmr_repro_cemt/scripts/submit_dmr.slurm
```

可用 `BCELLS_RSCRIPT` 为 Slurm 指定 Rscript。输入 CSV 和 GTF 保留在本地，GitHub 代码导出不包含它们。复跑需要先恢复到 `input/`。

`run_dmr.R` 只包装原程序执行，并保存会话版本、样本/设计矩阵、数量汇总和 BH 校正一致性核对。首次运行时输出目录为空，没有覆盖既有分析结果。

## 结果

- `output/all_bins_1kb_cemt_g2_g3_limma_results.csv`：全部保留窗口的效应、p/q 值与命中标记。
- `output/cemt_g2_vs_g3_dmrs_1kb_limma.csv`：满足原条件的窗口，按 |delta_beta| 降序。
- `output/run_summary.tsv`：检验数量、显著窗口数、方向与核对结果。
- `output/sample_metadata.tsv`、`output/design_matrix.tsv`、`output/session_info.txt`：设计和环境记录。
- `input_validation.json`：输入范围、覆盖情况、输入与原脚本 SHA256。
- `submission.json`、`logs/`：此次提交与运行记录。

## 本次复跑结果

Job `15896213` 正常完成（exit 0:0），耗时 2 分 22 秒。原分析脚本 SHA256 与运行前一致。

| 项目 | 窗口数 |
|---|---:|
| 输入 | 2,886,545 |
| 常染色体 | 2,875,012 |
| 至少 8 个样本不缺失，保留进入模型 | 2,697,259 |
| 有有效组效应 P 值 | 2,697,256 |
| FDR < 0.05 | 29 |
| FDR < 0.05 且绝对 delta_beta > 0.2 | **15** |
| Group2 / CLL-R 甲基化较高 | **11** |
| Group3 / CLL-S 甲基化较高 | **4** |

limma 提示 70 个窗口存在部分系数 NA；其中最终组效应 P 值缺失的窗口为 3 个，没有被判为显著。BH 校正值通过重新计算核对。完整运行核对见 `output/result_verification.json`；原参数、年龄/性别模型及四舍五入后应用效应阈值的行为均未改变。

## 火山图

`output/cemt_R_vs_S_dmr_volcano.pdf/.png` 使用带方向的 delta_beta（R−S）作为横轴，−log10(原始 P 值)为纵轴。全部有效窗口都参与绘图，没有随机抽样；PDF 中灰色背景点栅格化，以控制文件大小。

红色为 q<0.05 且 delta_beta>0.2，蓝色为 q<0.05 且 delta_beta<−0.2，橙色为 q<0.05 但效应未达阈值。竖虚线为 ±0.2。横虚线为本次 BH 显著项中最大的原始 P 值，是本数据集的边界，不是把 q=0.05 当作原始 P=0.05。

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python dmr_repro_cemt/scripts/plot_dmr_volcano.py
```

29 个 FDR 显著窗口的绘图坐标见 `output/volcano_FDR005_windows.tsv`，绘图口径与计数见 `output/volcano_metadata.json`。

## IGH 区域轨迹图

```bash
python dmr_repro_cemt/scripts/plot_igh_tracks_corrected.py
```

读取相同的 1 kb beta 矩阵及 `input/gencode.v47.annotation.gtf`，生成 `output/chr14_105800000_106500000_methylation_tracks_corrected.pdf/.png`。
