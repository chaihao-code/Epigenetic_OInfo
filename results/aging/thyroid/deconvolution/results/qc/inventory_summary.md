# Step 1/2 数据盘点摘要

- 样本: ['CEMT0040', 'CEMT0042', 'CEMT0044', 'CEMT0098', 'CEMT0130', 'CEMT0132'](原分析排除 ['CEMT0086', 'CEMT0088', 'CEMT0090'])
- bigWig 基因组版本(实测 chr21 长度): {'hg38'}
- 甲基化单位(实测): {'percent_0_100'}
- 缺失文件: 无
- 关键限制: methylation_profile bigWig 为 0-100 百分制且无 coverage, 无法运行 UXM_deconv; 采用 constrained NNLS beta 混合近似并单独报告。
- 零值语义: zero_frac 列见 input_inventory.tsv; bigWig 无值区间返回 NaN(缺失), 0 为实测信号零值; WGBS 未覆盖 CpG 处以 NaN/无记录表示。
