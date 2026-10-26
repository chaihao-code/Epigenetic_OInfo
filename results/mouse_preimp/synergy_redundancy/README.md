# 协同性和冗余性分析：5 bins，10000 次置换重算

按用户要求，原有 results_n_bins_5 下的 2240 个分区 O-information CSV 已删除，
使用共享 compute_hoi_oinfo.py 从原始信号重新计算，结果写回同一个目录。
calculation_results 下的既有熵和 mutual_information 下的数据、图表保留。

## 参数

7 个阶段：zygote、2cell、4cell、8cell、morula、ICM、TE；完整 9 信号系统。
5 bins，40 个 mm10 区间，每区间最多 100000 个位点，抽样 seed 42。
10000 次独立区间内置换，permutation seed 1234。
观测熵及 O-information 计算所有子系统；置换背景仅计算完整 9 信号系统。
没有启用删除 H2A.Z/H3K36me3 的约简系统，也没有在本任务中做 FDR 校正。

## 输出

每个阶段写入 results_n_bins_5/<stage>/：
- <stage>_entropy_all_subsets.csv：该次重算的全子集熵。
- <stage>_oinfo.csv：各区间、各阶子系统的观测 O-information、TC、DTC。
- <stage>_permutation_oinfo.csv：完整系统的逐次置换统计量。
- <stage>_permutation_summary.csv：完整系统观测值、背景统计量和经验 P 值。
- <stage>_metadata.json：运行参数。

新输出采用合并表结构，不再使用旧的逐分区 *_oinfo_missing_*.csv 文件结构。
hoi < 0 为 synergy，hoi > 0 为 redundancy；统计显著性需根据置换检验另行判定。

## 提交及追溯

从本目录运行：sbatch scripts/submit_oinfo.slurm。
任务数组包含 7 个任务，每任务 16 CPUs、24 小时，无并发上限。
Slurm 日志在 results_n_bins_5/logs/。
输入顺序在 scripts/oinfo_samples.tsv；参数、源代码 SHA256 和任务号在 rerun_10000_manifest.json。
