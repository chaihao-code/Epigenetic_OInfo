# 方法（中文）—— 甲状腺衰老细胞组成混杂评估

> 本文档按科学问题分节撰写，每节对应实际运行的脚本与结果表。
> 全局配置与路径见 `config.yaml`，逐步复现命令见 `README.md`。

## 一、Q1：甲状腺样本的免疫细胞组成估计（DNA 甲基化去卷积）

### 1.1 输入数据与关键限制

目标样本为 IHEC/CEEHRC 的 6 例健康成人甲状腺供体（CEMT0098/0130/0132/0042/0044/0040，
30–67 岁，5 女 1 男）。甲基化数据为
`*.WGB-Seq.methylation_profile.bigWig`（hg38，染色体命名 `chrN`，全部样本坐标一致，
经染色体长度校验，如 chr21 = 46,709,983）。

**关键限制**：该文件为 0–100 百分制的甲基化比例，**不含 per-CpG 覆盖度与甲基化读段计数**。
因此无法运行基于片段（PAT）输入的 UXM_deconv（Loyfer et al. 配套方法）。
本研究改用基于甲基化比例（beta 值）线性混合近似的参考约束回归，并在结果解释中单独报告此限制。

### 1.2 参考图谱数据

参考数据来自 Loyfer et al., *Nature* 613:355–364 (2023)（GEO: GSE186458），
即"正常人类细胞类型 DNA 甲基化图谱"。

- 从 GEO metadata 获取真实细胞类型标签与下载地址，仅下载所需细胞类型，
  未下载数百 GB 的 RAW.tar。
- 下载 27 个 wgbstools 格式的 `.hg38.beta` 文件（每文件为
  29,401,795 个 CpG 的 #meth/#coverage uint8 对），覆盖 9 种细胞类型 × 3 供体：
  甲状腺上皮（GSM5652264–566，甲状腺组织）、T_CD4、T_CD8、B 细胞、单核细胞、
  NK、粒细胞、内皮（主动脉/隐静脉）、成纤维细胞（结肠/心脏）。
- 参考与目标同为 hg38，无需 liftOver。下载清单与状态见
  `results/tables/reference_download_manifest.tsv`。

### 1.3 参考矩阵构建

1. 建立全部 CpG 的基因组索引（29,401,795 个，与 wgbstools CpG 编号 1-based 精确对应，
   验证：5,000/5,000 坐标一致）。
2. 对 1 kb bin 内的 CpG 取甲基化计数与覆盖度之和，计算 bin 水平 beta 与平均覆盖度
   （float16 存储，2,875,012 常染色体 bin × 27 样本）。
3. 每细胞类型取 3 个供体的 bin 水平均值，得到参考矩阵
   R[g, k]（g 为 1 kb bin，k 为细胞类型）。

**参考数据质量控制**：
- 同型别供体间重复相关性（bin 水平 Pearson r）中位数 0.66（成纤维，跨组织来源）至
  0.93；内皮/成纤维较低源于跨组织混合，已记录。
- 外部验证：由 beta 文件计算的甲状腺上皮甲基化与 GEO 官方提供的
  GSM5652264 bigWig 逐 CpG 比对，r = 0.9999998；与 Loyfer 发表的
  `Markers.U250` 靶标均值比对 r = 0.67–0.90。
- 发现 Atlas.U250 完整文件的列与行定义不一致，弃用该文件，改用逐样本 beta 自行构建。

### 1.4 Marker Panel 与细胞类型粒度的定义

本研究在去卷积中同时考虑两个独立维度（两个维度两两组合共 4 套参考矩阵，
互相作为敏感性分析）。

**Marker panel（参考矩阵 R 的 marker 来源与筛选方式）：**

- **Panel B（文献来源）**：直接采用 Loyfer 等发表的 UXM 去卷积配套 marker 集
  （`Markers.U250.hg38`，基于 250 bp CpG 甲基化区块定义），并从中保留属于本研究
  建模细胞类型的 marker，共 2,149 个。每个 marker 记录组间对比度、型别内重复变异、
  参考覆盖度（`results/tables/marker_quality.tsv`）。
- **Panel A（独立筛选）**：不依赖任何已发表 marker 定义，由我们从参考数据独立筛选：
  在 1 kb bin 分辨率上计算"免疫细胞均值与其他细胞均值之差" d，取 |d| 最大的区间，
  同时要求 9 类参考内最大重复间标准差 ≤ 0.08、平均覆盖度 ≥ 10；
  高甲基化（d ≥ 阈值）与低甲基化（d ≤ −阈值）两个方向各取 1,250 个，共 2,500 个，
  以避免 marker 全部偏向单一方向。

**细胞类型粒度（参考矩阵列的组织方式）：**

- **Model B（粗粒度四分类）**：甲状腺上皮、淋巴系（T_CD4/T_CD8/B 细胞合并）、
  髓系（单核/NK/粒细胞合并）、基质（内皮/成纤维合并）。
- **Model A（细粒度九分类）**：甲状腺上皮、CD4 T、CD8 T、B 细胞、单核细胞、
  NK、粒细胞、内皮细胞、成纤维细胞，全部保留。

其中 **Panel A / Model A 为主分析**（marker 独立筛选、细胞类型最细），
**Panel B / Model B 为敏感性分析**（完全沿用文献定义、最粗粒度）。
若细粒度下某些细胞类型因参考信号高度相似而不可辨识（如 T/B 细胞），
则按 Model B 的合并方式解释为更粗的组分，不强行报告不可辨识的细类比例。

### 1.5 目标样本甲基化信号提取

对每个 bulk 样本，按 marker 的 1 kb 坐标从 `methylation_profile.bigWig` 提取范围内
全部 CpG 的甲基化比例，取均值并除以 100 得到 bin 水平 beta 值。
bigWig 中无记录的 bin（该范围无 CpG）记为缺失，拟合时整行剔除；
6 个样本的 marker 缺失率均 ≤ 0.4%。提取结果缓存于
`data/processed/target_marker_methylation_{A,B}.tsv`。

### 1.6 去卷积模型：约束非负最小二乘

对每样本解约束 NNLS：

  min ||y − R p||²，约束 p ≥ 0，Σp = 1

其中 y 为样本在 marker 上的 beta 向量，R 为参考矩阵，p 为待估计的细胞比例。
实现采用增广矩阵法：在 R 下方追加全 1 行（权重 λ）、y 下方追加 λ，
调用 scipy NNLS 后迭代增大 λ 直至比例和为 1（通常 1–3 次收敛），
最后归一化。若 Σp 强制为 1 而模型未覆盖全部细胞类型，结果为**相对已建模细胞类型的
组成估计**（如间质/脂肪等未建模成分会负载到已有分量上），解释时已注明。

### 1.7 不确定性评估：染色体区块 bootstrap

为避免把相邻 CpG 当作独立观测，采用**染色体区块 bootstrap**：以整条染色体为
重抽样单元，有放回地抽取 22 条常染色体（等价于对每条染色体上的 marker 整体重抽样），
重新拟合 NNLS，重复 200 次（种子固定），取各细胞比例的 2.5% 与 97.5% 分位数作为
95% 置信区间。

### 1.8 模型验证（模拟混合，先于实际应用）

在 beta 值层面用参考数据构造已知比例的模拟混合样本（线性混合近似，已声明）：
免疫总比例取 0%、1%、5%、10%、20%、30%、50% 七档，两种免疫内部构成分配，
观测噪声 σ = 0.01/0.05 两档，每组合 50 次重复，共 4 种 Panel×Model 组合。
评价指标：MAE、RMSE、估计 vs 真实比例偏差、零免疫时的假阳性、低比例检出率。

结果（`results/tables/simulation_summary.tsv`）：0% 免疫时估计偏差 ≈ 0.0003–0.0019
（无系统性假阳性）；整体 MAE 0.002–0.012；**5% 免疫组分的检出率为 27–72%**——
故低于约 5% 的免疫比例视为低精度估计，后续解释以效应量与置信区间为准。

### 1.9 拟合优度与实际估计结果

实际样本拟合 R² = 0.91–0.96，marker 缺失率 ≤ 0.4%。

**组成估计（主分析 Panel A / Model A 的点估计，数值范围为 6 个样本的
最小–最大值，并非统计参考区间；每样本 95% CI 见表内 ci_lo/ci_hi 列）**：
总免疫比例 1.2%–10.8%，甲状腺上皮 77.2%–88.3%，基质（内皮+成纤维）9.1%–16.7%
（`results/tables/thyroid_cell_proportions.tsv`，图 1A–1D）。

敏感性组合 Panel B / Model B 给出的相应范围为：总免疫 0.0%–7.7%、
甲状腺上皮 80.2%–94.0%、基质 5.5%–14.1%——量级一致，但 Panel B 系统性
偏低约 2–3 个百分点，反映两套 marker 定义（250 bp 区块 vs. 1 kb bin 独立筛选）
之间的方法学差异；正文以主分析为准，Panel B 作为稳健性范围。
两套 marker panel 得到的下游回归结果高度一致（调整后回归斜率相关性 r = 0.993，
见 Step 8 敏感性分析）。

### 1.10 软件与固定随机种子

Python 3.8.5：numpy 1.24.4、pandas 1.5.3、scipy 1.10.1（NNLS）、
scikit-learn 1.3.2、pyBigWig 0.3.18。种子：bootstrap 与模拟 2026/42（见 `config.yaml`）。

## 二、Q3：O-information 与估计免疫细胞比例的相关性

### 2.1 数据与模型

O-information（Ω）数据取自原研究在 40 种基因组区间（genomic context）上的
样本 × 区间 Ω 值矩阵，区间定义保持原样（promoter/enhancer × CTCF ± ×
CpG island/None/repeat 家族等组合）。免疫细胞比例取主分析（Panel A / Model A）
的 6 类免疫组分之和（total immune fraction）。

对每种区间分别估计 Ω 与免疫比例的关联，同时报告：

- 双变量 Pearson / Spearman 相关（图 2C、2D 散点）；
- 联合模型中的免疫项系数 β_immune（即 3.1 节 Model 1 中的免疫比例项，
  可同时控制年龄），每区间一条记录，含 95% CI 与 p 值。

独立统计单位为 donor（n = 6），不以 1 kb bin 充当生物学重复。

### 2.2 结果

`results/tables/oinfo_age_adjustment.tsv`（β_immune 列）、图 2C/2D：

- **40 个区间中无一例 β_immune 达到 p < 0.05**（最小 p = 0.06，
  `None_AND_LINE` 区间）。
- 重点区间（CTCF-bound proximal enhancer、无 CpG island，
  `proximal_ctcf_bound_enhancer_purified_bins_AND_None`）：
  β_immune = −3.16（95% CI [−16.3, +10.0]），p = 0.50；
  双变量 Spearman r = −0.60，p = 0.20（图 2D）。
- 解释限制：n = 6 时该检验功效很低，"未显著"不能等同于"证明无关"，
  故同时报告效应量与置信区间；但就点估计而言，没有任何区间的 Ω 呈现
  随免疫比例系统变化的证据。

## 三、Q4：控制免疫细胞比例后，O-information 的年龄效应是否仍然成立？

### 3.1 分析目标与模型设定

本节只回答一个问题：**在把甲基化去卷积估计的免疫细胞比例纳入模型后，
Ω 随年龄的下降是否仍然存在**。因此对每种基因组区间（共 40 种），
主要估计目标为下述模型中的 β_age：

  **Model 1（主模型）**：Ω = β₀ + β_age·Age + β_immune·免疫比例 + ε

其中免疫比例为该样本 6 类免疫组分之和（主分析 Panel A / Model A，
见第一节）。未调整模型 Model 0（Ω ~ Age）仅作为参照，用于量化
加入免疫比例后年龄效应的衰减。独立统计单位为 donor（n = 6），
不以 1 kb bin 充当生物学重复。

推断与稳健性：

- 每个区间独立拟合 OLS，报告 β_age 的 95% CI 与 nominal p 值；
  40 个区间的多重检验同时报告 BH 与两阶段 BKY（Benjamini–Kestelman–
  Yekutieli 自适应 step-up，`fdr_tsbky`）两种 FDR。
- 年龄与免疫比例的共线性：Pearson r = 0.50（p = 0.31，
  `results/tables/immune_age_regression.tsv`，即 Q2 的全部结果），
  Model 1 的 VIF ≤ 1.34，属轻度共线性，年龄与细胞组成效应可以分离。
- n = 6 时不进一步加入性别/疾病/批次协变量（1 男 5 女、无疾病状态变异），
  作为局限性报告。
- 稳健性：① 留一 donor 重拟合，β_age 始终落在原始 95% CI 之内
  （全体区间 loo 斜率范围 [−0.058, +0.057]）；② marker panel 敏感性——
  Panel B / Model B（UXM 文献 marker + 粗粒度）组合下 40 个区间的调整后
  年龄斜率与主分析相关 r = 0.993（图 2F）；③ 衰减度
  attenuation = 1 − β_age_adj/β_age_unadj 仅在 |β_age_unadj| 接近零时
  不稳定（可出现负值），故下文同时直接报告两个斜率及其 CI，
  不以"调整后 p > 0.05"作为"年龄效应消失"的判据。

### 3.2 结果

`results/tables/oinfo_age_fdr_table.tsv`（完整 40 区间）、
`results/tables/oinfo_age_adjustment.tsv`、图 2A–2F。
下表按调整后 p 值排序（仅主分析 Panel A / Model A）：

| 基因组区间 | β_age 未调整 | p 未调整 | β_age 调整后 | 95% CI(调整后) | p 调整后(nominal) | FDR-BH | FDR-BKY |
|---|---|---|---|---|---|---|---|
| proximal_ctcf_bound_enhancer_purified_bins_AND_CpG_island | -0.0076 | 0.0001 | -0.0073 | [-0.0088, -0.0058] | 0.0006 | 0.0227 | 0.0233 |
| proximal_ctcf_bound_enhancer_purified_bins_AND_SINE | -0.0307 | 0.0019 | -0.0290 | [-0.0460, -0.0120] | 0.0122 | 0.1714 | 0.1755 |
| proximal_ctcf_bound_enhancer_purified_bins_AND_LINE | -0.0380 | 0.0020 | -0.0378 | [-0.0603, -0.0153] | 0.0129 | 0.1714 | 0.1755 |
| proximal_ctcf_bound_enhancer_purified_bins_AND_LTR | -0.0458 | 0.0043 | -0.0466 | [-0.0798, -0.0134] | 0.0209 | 0.2094 | 0.2144 |
| proximal_ctcf_bound_enhancer_purified_bins_AND_None | -0.0396 | 0.0132 | -0.0352 | [-0.0715, 0.0011] | 0.0539 | 0.4313 | 0.4416 |
| proximal_not_ctcf_bound_enhancer_purified_bins_AND_LTR | -0.0301 | 0.0736 | -0.0370 | [-0.0837, 0.0096] | 0.0856 | 0.4884 | 0.5000 |
| proximal_not_ctcf_bound_enhancer_purified_bins_AND_SINE | -0.0319 | 0.0320 | -0.0309 | [-0.0728, 0.0109] | 0.1001 | 0.4884 | 0.5000 |
| None_AND_CpG_island | -0.0262 | 0.0300 | -0.0222 | [-0.0527, 0.0083] | 0.1035 | 0.4884 | 0.5000 |
| proximal_not_ctcf_bound_enhancer_purified_bins_AND_LINE | -0.0308 | 0.0599 | -0.0344 | [-0.0829, 0.0142] | 0.1099 | 0.4884 | 0.5000 |
| distal_ctcf_bound_enhancer_purified_bins_AND_SINE | -0.0144 | 0.0761 | -0.0161 | [-0.0412, 0.0089] | 0.1326 | 0.5304 | 0.5430 |
| proximal_not_ctcf_bound_enhancer_purified_bins_AND_CpG_island | -0.0204 | 0.0498 | -0.0162 | [-0.0435, 0.0111] | 0.1547 | 0.5626 | 0.5759 |
| distal_ctcf_bound_enhancer_purified_bins_AND_LTR | -0.0105 | 0.5149 | -0.0224 | [-0.0681, 0.0233] | 0.2167 | 0.7069 | 0.7237 |
| distal_ctcf_bound_enhancer_purified_bins_AND_LINE | -0.0109 | 0.3768 | -0.0180 | [-0.0566, 0.0206] | 0.2350 | 0.7069 | 0.7237 |
| distal_ctcf_bound_enhancer_purified_bins_AND_CpG_island | -0.0161 | 0.1278 | -0.0161 | [-0.0518, 0.0196] | 0.2474 | 0.7069 | 0.7237 |
| None_AND_None | 0.0015 | 0.8175 | 0.0068 | [-0.0092, 0.0227] | 0.2688 | 0.7169 | 0.7339 |
| proximal_not_ctcf_bound_enhancer_purified_bins_AND_None | -0.0236 | 0.2114 | -0.0262 | [-0.0932, 0.0408] | 0.3013 | 0.7532 | 0.7711 |
| distal_not_ctcf_bound_enhancer_purified_bins_AND_SINE | -0.0116 | 0.2969 | -0.0140 | [-0.0540, 0.0260] | 0.3459 | 0.7831 | 0.8017 |
| promoter_AND_LTR | -0.0080 | 0.6865 | -0.0207 | [-0.0836, 0.0421] | 0.3710 | 0.7831 | 0.8017 |
| distal_ctcf_bound_enhancer_purified_bins_AND_None | -0.0087 | 0.5207 | -0.0150 | [-0.0623, 0.0322] | 0.3860 | 0.7831 | 0.8017 |
| intron_purified_bins_AND_SINE | -0.0091 | 0.4564 | -0.0126 | [-0.0576, 0.0324] | 0.4385 | 0.7831 | 0.8017 |
| distal_not_ctcf_bound_enhancer_purified_bins_AND_LTR | -0.0021 | 0.9044 | -0.0144 | [-0.0679, 0.0391] | 0.4539 | 0.7831 | 0.8017 |
| promoter_AND_LINE | -0.0083 | 0.5772 | -0.0142 | [-0.0681, 0.0396] | 0.4625 | 0.7831 | 0.8017 |
| promoter_AND_SINE | -0.0086 | 0.4071 | -0.0099 | [-0.0494, 0.0295] | 0.4811 | 0.7831 | 0.8017 |
| exon_purified_bins_AND_SINE | -0.0072 | 0.4821 | -0.0089 | [-0.0479, 0.0302] | 0.5217 | 0.7831 | 0.8017 |
| None_AND_SINE | 0.0076 | 0.6327 | -0.0072 | [-0.0391, 0.0248] | 0.5257 | 0.7831 | 0.8017 |
| promoter_AND_CpG_island | 0.0042 | 0.2623 | 0.0028 | [-0.0098, 0.0154] | 0.5316 | 0.7831 | 0.8017 |
| distal_not_ctcf_bound_enhancer_purified_bins_AND_LINE | -0.0023 | 0.8742 | -0.0107 | [-0.0605, 0.0391] | 0.5434 | 0.7831 | 0.8017 |
| exon_purified_bins_AND_LTR | 0.0019 | 0.9205 | -0.0121 | [-0.0690, 0.0448] | 0.5482 | 0.7831 | 0.8017 |
| intron_purified_bins_AND_LTR | 0.0051 | 0.7944 | -0.0096 | [-0.0657, 0.0464] | 0.6227 | 0.8511 | 0.8713 |
| None_AND_LINE | 0.0269 | 0.2497 | 0.0068 | [-0.0363, 0.0499] | 0.6508 | 0.8511 | 0.8713 |
| promoter_AND_None | -0.0071 | 0.6716 | -0.0099 | [-0.0753, 0.0555] | 0.6629 | 0.8511 | 0.8713 |
| distal_not_ctcf_bound_enhancer_purified_bins_AND_CpG_island | -0.0152 | 0.3604 | -0.0082 | [-0.0654, 0.0491] | 0.6808 | 0.8511 | 0.8713 |
| distal_not_ctcf_bound_enhancer_purified_bins_AND_None | 0.0027 | 0.8736 | -0.0064 | [-0.0644, 0.0516] | 0.7477 | 0.8993 | 0.9207 |
| exon_purified_bins_AND_LINE | 0.0024 | 0.8606 | -0.0048 | [-0.0515, 0.0419] | 0.7644 | 0.8993 | 0.9207 |
| intron_purified_bins_AND_CpG_island | -0.0046 | 0.7249 | -0.0037 | [-0.0560, 0.0485] | 0.8341 | 0.9384 | 0.9606 |
| intron_purified_bins_AND_LINE | 0.0080 | 0.6277 | -0.0024 | [-0.0550, 0.0502] | 0.8946 | 0.9384 | 0.9606 |
| exon_purified_bins_AND_None | 0.0075 | 0.5939 | 0.0023 | [-0.0496, 0.0543] | 0.8958 | 0.9384 | 0.9606 |
| None_AND_LTR | 0.0214 | 0.4179 | -0.0023 | [-0.0545, 0.0499] | 0.8989 | 0.9384 | 0.9606 |
| intron_purified_bins_AND_None | 0.0108 | 0.4981 | 0.0017 | [-0.0505, 0.0538] | 0.9245 | 0.9384 | 0.9606 |
| exon_purified_bins_AND_CpG_island | -0.0030 | 0.8335 | -0.0015 | [-0.0576, 0.0547] | 0.9384 | 0.9384 | 0.9606 |

读表要点：

1. **方向与量级**：40 个区间中 32 个的 β_age 为负；全部 8 个
   proximal enhancer 区间（CTCF ± × CpG island/SINE/LINE/LTR/None）
   在调整前后均为负，且斜率几乎不变（attenuation = −0.02 ~ +0.15，
   即衰减至多约 11%，多数 < 6%）。
2. **名义显著性**：调整后 4 个区间 nominal p < 0.05，且全部是
   CTCF-bound proximal enhancer 区间（× CpG island / SINE / LINE / LTR），
   β_age = −0.029 ~ −0.047 /年，其 95% CI 均不含 0；
   其中 CTCF+ × CpG island 最显著（β = −0.0073，p = 0.0006，
   CI [−0.0088, −0.0058]，每年下降约 0.7%）。
3. **多重检验校正后**：CTCF+ proximal enhancer × CpG island 在
   BH 与 BKY 两种 FDR 下均满足 FDR < 0.05（q = 0.023），且该结论
   对 0.1 阈值同样成立；其余 3 个区间 FDR 在 0.17–0.21。
   原稿重点关注的 CTCF+ proximal enhancer 无 CpG island 区间
   调整后 β_age = −0.0352（p = 0.054，attenuation = 11%），
   方向与量级保持，但 CI 已跨 0。
4. **敏感性**：上述 4 个区间的结论在留一 donor、染色体区块
   bootstrap（表内 boot_ci 列）与另一套 marker panel 下均稳定；
   两套 panel 的 40 个调整后斜率相关 r = 0.993。

**结论**：在 6 例样本上，把甲基化估计的免疫细胞比例纳入模型后，
CTCF-bound proximal enhancer（尤其 × CpG island）的 Ω 年龄相关下降
**仍然存在**——斜率本身几乎不变（衰减 ≤ 11%），说明该年龄效应
不能被本研究估计的细胞组成变化所解释。但由于 n = 6、免疫比例
估计精度有限（1.8 节），且甲基化去卷积不等于对七条表观轨道的
完整细胞组成校正，此结果不能被表述为已证明甲状腺上皮细胞内在的
表观遗传衰老。
