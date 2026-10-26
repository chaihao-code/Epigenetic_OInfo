# Thyroid local GSEA 作图重跑

两个原始入口，共572幅图、1144个新PDF/PNG文件。历史曲线清单为pathway 204幅、marginal 368幅。

| 原始入口 | 内容 | 状态 | 新图/文件数 | 最大ES误差 |
|---|---|---|---:|---:|
| `02_gsea_pathway.py` | 甲状腺局部O-information年龄效应：GO BP/Hallmark GSEA曲线 | success | 204/408 | 4.98e-08 |
| `03_gsea_marginal.py` | 甲状腺局部O-information年龄效应：各信号up/down/flat GSEA曲线 | success | 368/736 | 4.99e-08 |

所有40个注释区域、80个pathway区域/库及40个marginal区域均纳入检查。仅使用原历史PDF文件名选择应重画terms，没有读取或复制其图片内容。原始report/RNK/GMT的来源、目标和SHA256记录在[restored_inputs.json](restored_inputs.json)。全部输出位于`cases/aging/thyroid/local_O_info_regression/results/gsea_pathway`及`gsea_marginal`。

每条曲线从完整RNK和GMT计算weighted running ES，断言与保存ES差异≤2e-6；报告中保存的NES/P/FDR保持原值。使用matplotlib重新布局与渲染，不重跑回归、基因注释、信号提取或GSEA置换。所有PNG和PDF逐文件验证。

原`02_gsea_pathway.py`需要回归表及GENCODE GTF，并执行基因排名与置换；原`03_gsea_marginal.py`需要回归表及提取信号，执行bin分类、基因集构建及置换。因此这两个原入口不是纯绘图程序。本次专门提供只读中间结果的作图入口：

```bash
python validation/plot_rerun_20261026/gsea/run.py --jobs 3
```

逐图term、ES核验、输入组和结果SHA256在[report.json](report.json)。三个绘图进程各限制BLAS/OMP为1线程；每组执行证据在[groups_results](groups_results)。

[全部572幅缩略图](contact_sheet_all.jpg) · [24幅代表图预览](contact_sheet_examples.jpg) · [预览索引](preview_index.json)
