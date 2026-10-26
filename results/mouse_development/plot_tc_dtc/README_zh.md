# O-information / TC / DTC 发育轨迹图

总 O-information 绘图入口：`plot_total_oinfo_trajectories_18cm.py`。

- `combined_oinfo_trajectories_3x2_n5.pdf/.png`：七信号总系统 O-information，18 × 8 cm；
  沿用下列 TC/DTC 图的六类区域、82 个样本、时间坐标、组织颜色及实虚线。
- 直接读取 `tc_dtc_overall_values_for_plot.tsv` 的 `hoi` 列，逐条核对 `hoi = tc - dtc`，
  不将 triplet 平均值或局部 O-information 作为总系统结果。
- 主标题 10 pt，其余文字 8 pt。保持右侧图例，合并重复的植入前图例条目并减小行距以适应压缩高度。
- `oinfo_overall_values_for_plot.tsv` 保存实际绘图数值，
  `combined_oinfo_trajectories_3x2_n5_validation.json` 保存来源、指标核对及排版检查。

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /lustre/home/2101110354/anaconda3/envs/hoi_env/bin/python3 project/mouse_development/plot_tc_dtc/plot_total_oinfo_trajectories_18cm.py
```

绘图入口：`plot_tc_dtc_stacked_18cm.py`。

- `combined_tc_trajectories_3x2_n5.pdf/.png`：独立 TC 图，18 × 9.75 cm。
- `combined_dtc_trajectories_3x2_n5.pdf/.png`：独立 DTC 图，18 × 9.75 cm。
- `combined_tc_dtc_trajectories_6x2_n5.pdf/.png`：TC 在上、DTC 在下，18 × 20 cm，共用一套右侧图例。

TC、DTC 两个总标题为 Arial 10 pt，其余文字为 8 pt。两部分中间间隔 0.5 cm。
区域标题统一采用 `Distal enhancer CTCF- / CpG island` 格式。
每个绘图区宽 6.3 cm，为图例保留更宽空间，右侧至少留 0.4 cm 边距。
时间标签采用单行竖排。六个区域、发育树、时间点、线条颜色、实虚线和 TC/DTC 数值沿用原图。
为使收窄后的图仍能容纳 8 pt 标签，植入前/原肠胚的分类间隔为 1.25 单位；
器官发生阶段保留每天两个单位，ICM 和 TE 共用 E3.5 位置。
图例中较长的组织名分行显示。

数据读取 `../../permutation_test/2026-09-28/results_mouse_development/` 中已有的
七信号总系统结果（order=7、n_missing=0），不重算熵或 TC/DTC。
保留原图覆盖的 82 个样本和原始缺失样本记录，来源与
`../../mouse_development_fine_te/results_oinfo_n_bins_5/` 的对应区域结果进行数值核对。

`tc_dtc_overall_values_for_plot.tsv` 保存实际使用的原始数值；
`combined_tc_dtc_trajectories_6x2_n5_validation.json` 保存输入哈希、区域核对、布局、字号、
时间标签边界及曲线坐标/样式检查。`tc_dtc_missing_inputs.csv` 保留原图的覆盖记录。

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /lustre/home/2101110354/anaconda3/envs/hoi_env/bin/python3 project/mouse_development/plot_tc_dtc/plot_tc_dtc_stacked_18cm.py
```
