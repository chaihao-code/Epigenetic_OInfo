# O-information 重算验证报告 (Step 9 前置 QC)

日期: 2026-10-09

## 验证设计

用 `scripts/hoi_recompute.py` 从磁盘上的离散化矩阵
(`data/aging_thyroid_input_normed_binned_discreted/*.discrete.bed`, 2026-06-01 生成)
严格按原管线逻辑(`compute_hoi_entropy.py` + `compute_oinfo_from_entropy.py` +
`block_processor.cartesian_product_blocks`)重算全部 6 样本 × 40 区域的全变量
O-information(hoi = tc − dtc),与 `results_oinfo_n_bins_5` 中存储的原始结果逐值比对。

## 结果

1. **区域定义与 bin 宇宙完全一致**: 40 个区域的 bin 数与原始 metadata 中的
   n_samples 全部一致(13 个大区域原始端经 n_max=100000, seed=42 子采样,
   本实现复现同一 pandas `df.sample(random_state=42)` 协议)。
2. **DNA methylation 轨道熵逐位精确一致**(例: proximal_ctcf_bound_enhancer_purified_bins_AND_None,
   H(WGBS) = 1.711810 vs 存储 1.711810)。
3. **6 个组蛋白轨道熵不一致**: 表明 2026-06-24 原始熵计算所用的 divInput
   bin-average 输入文件(`.divInput.bp`,位于已失效的 /media/dell/data4)与
   2026-06-01 生成离散化矩阵所用的输入不同(文件已全盘搜索确认不可恢复;
   两种 assign_bins 实现逐行一致,排除离散化代码差异)。
4. 存储 hoi  vs 重算 hoi(240 对): Pearson r = 0.848 (p = 1.4e-67),
   平均 |Δ| = 0.38 bit。

## 决策(记录)

- Step 9 的"原始 / 过滤 / 匹配删减"三个臂**统一使用本重算实现与 6 月 1 日离散化
  矩阵**,保证三臂估计器完全一致,内部可比;回答的科学问题是年龄效应是否在
  免疫区域排除后稳健,估计器的绝对值差异不影响三臂对比。
- 论文中的原始 O-info 值(Step 8 调整分析)仍使用存储的原始结果,不重算。
- 存储的组蛋白离散输入不可恢复,
  重算值与存储值高度相关(r=0.85)但不完全相同。
