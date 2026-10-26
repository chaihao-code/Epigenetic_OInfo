# 男性差异基因在女性 CLL-R 样本中的表达一致性

结论：多数基因的变化方向可以延续到女性，但四位女性并非都完整复现男性 R 组的表达强度。

## 样本与检验范围

现有 RNA 数据包括男性 CLL-R 4 人、男性 CLL-S 5 人，以及女性 CLL-R 4 人：CEMT0005、CEMT0025、CEMT0028、CEMT0030。没有女性 CLL-S，因此没有计算女性内部的 R/S 差异 p 值或 FDR，也不能将下列结果称为女性中显著差异的独立复现。

固定采用 `results_rna_de_male_ighv` 的男性 FDR < 0.05 名单（38 个基因），另以 FDR < 0.10 的 60 个基因及去除性染色体后的 37 个基因作敏感性分析。女性数据不参与选基因、估计归一化参考、拟合协变量或定义男性参照中心。

## 主要方法

1. 对既有 tximport lengthScaledTPM counts 取整；按男性原分析口径预筛（至少 4 位男性 count >= 10），在九位男性均有正计数的基因上建立冻结的几何均值参考，以 median-of-ratios 风格的 size factor 归一化所有样本。
2. 在 log2(normalized counts + 1) 上，仅用男性拟合 `intercept + RNA method + IGHV + R_vs_S`，去除 RNA method 与 IGHV 对表达的估计贡献。
3. 使用男性的基因均值和标准差做 z-score，再将女性投射到固定的男性 S→R 中心连线。0 表示男性 S 中心，1 表示男性 R 中心；0.5 是等距离位置。这是描述性表达相似度，不是诊断分类器。此对数表达投射模型与原 DESeq2 差异检验不同；所有 38 个基因在两种男性模型中的 R/S 方向相同。

## 主要结果

| 女性样本（实际全部为 R） | 对男性 R/S 差异同向的基因数 | S→R 投射位置 |
|---|---:|---:|
| CEMT0005 | 33/38 | 0.488 |
| CEMT0025 | 35/38 | 0.561 |
| CEMT0028 | 29/38 | 0.438 |
| CEMT0030 | 35/38 | 0.819 |

这里的“同向”指女性相对调整后的男性 S 中心，方向与原男性 DESeq2 R/S log2FoldChange 相同，不表示该女性样本存在统计显著差异。

- 四位女性的平均表达相对男性 S：**37/38 个基因同向**，唯一均值反向的是 **LPL**（男性 R 高，女性平均低于调整后的男性 S 中心；只有 1/4 位女性同向）。
- **25/38 个基因在四位女性中全部同向**。
- CEMT0030 的整体模式最接近男性 R；CEMT0025 略偏 R；CEMT0005 接近两组中间；CEMT0028 略偏 S。两位接近中间的样本不能据此改分临床组。
- 扩展到 60 个基因时，女性均值有 **59/60 同向**，38 个基因在四位女性中都同向。

## 敏感性与限制

38 基因名单在校正 log-count、校正 log-TPM、去除 PNMA3（chrX）后，都有 2/4 位女性更接近男性 R。仅使用六位 stranded 男性作参考并调整 IGHV 时为 3/4。60 基因名单在不同参照下为 2–4/4。结果支持总体方向的延续，但个体距离判定存在方法敏感性。

女性只有四人且全部为 R，RNA 方法和 IGHV 在男性组间不平衡；因此不能验证女性中的组间特异性、显著性或分类准确率。原男性 DE 基因在选取它们的九位男性中区分明显，本身不是额外验证证据。

## 文件与重算

- `female_male_signature_concordance.pdf/.png`：38 基因的调整后表达热图及四位女性的 S→R 投射位置。
- `female_signature_concordance.tsv`：四位女性在所有基因名单/敏感性模型下的逐样本结果。
- `gene_direction_validation.tsv`：60 个基因的男性原始 DESeq2 结果、女性逐人方向与 TPM；`in_male_FDR05` 标记主要 38 基因。
- `all_sample_signature_scores.tsv`、`primary_adjusted_male_scaled_expression.tsv`：全部 13 个样本的分数和主要表达矩阵。
- `validation_summary.tsv`：不同模型与基因名单的汇总。
- `*_projection_coefficients.tsv`：用于冻结投射的男性模型系数。
- `sample_metadata.tsv`、`validation_config.json`、`verification.json`：样本、方法、输入哈希与数值核对。

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 project/aging/bcells/.rna_de_env/bin/python \
  project/aging/bcells/scripts/validate_male_de_in_females.py
```

结果写入本目录，未修改原男性 DE 结果。
