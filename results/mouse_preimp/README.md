# Mouse preimplantation analysis

目录分为三类；最初仅整理文件和路径；随后按用户要求删除并重新提交了 5 bins O-information，详见 synergy_redundancy/README.md。
7 个阶段：zygote、2cell、4cell、8cell、morula、ICM、TE。
9 个信号：ATAC、H2A.Z、H3K4me1、H3K4me3、H3K9me3、H3K27ac、H3K27me3、H3K36me3、CpG methylation。

| 分类 | 目录 | 内容 |
|---|---|---|
| 计算结果 | `calculation_results/` | 输入清单；3/5/10 bins 全子集熵；5 bins 合并熵表；计算提交脚本 |
| 互信息分析 | `mutual_information/` | 3/5/10 bins 自身熵、两两 MI、联合 MI；5 bins 相关性；绘图代码、图片和源数据 |
| 协同性和冗余性分析 | `synergy_redundancy/` | 5 bins 的 O-information、TC、DTC，包含完整系统和各子系统；复现入口 |

## 复现依赖与结果含义

`calculation_results/data_info/*.txt` 指向项目外的原始 1 kb 信号；原始信号没有移动或删除。
共有计算代码在仓库根目录的 `script/hoi_code/`；mm10 注释在 `data/genome_features/mm10/`。
计算环境为 `/lustre/home/2101110354/anaconda3/envs/hoi_env/bin/python3`。
主结果为 5 bins，3/10 bins 用于分箱敏感性分析。

O-information < 0 为 synergy，> 0 为 redundancy。5 bins O-information 现按 10000 次区间内置换重新计算，观测值覆盖所有子系统，
置换背景覆盖完整 9 信号系统；新任务没有自动进行 FDR 校正。熵和 MI 既有结果保留。
互信息绘图入口位于 `mutual_information/scripts/`，输出在 `mutual_information/plots/`。
历史日志保留原运行路径；可执行脚本已更新为当前路径。

## 整理记录

replacement 结果目录、压缩包、7 份配置、配置生成程序和专用提交/汇总脚本及专用日志已删除。
只清理本目录中的 replacement 内容，其他项目的独立分析未删除。
`organization_manifest.json` 记录旧路径、新路径、文件大小、数据/图片校验值以及校验结果。

## GitHub source export

See [UPLOAD_NOTES.md](UPLOAD_NOTES.md) for the uploaded scope, excluded legacy/TF analyses,
shared calculation snapshot, original HPC paths and external data requirements.
