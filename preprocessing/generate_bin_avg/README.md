# generate_bin_avg — 基因组信号分箱与格式转换工具集

## 功能概述

本目录包含两套基因组信号处理流水线：

1. **CpG 甲基化数据转换**：将 Bismark / NOMe-seq 输出的 per-CpG BED 文件转换为连续的 bedGraph 格式
2. **基因组信号分箱**：将 bigWig 或 bedGraph 格式的连续信号按固定窗口（bin）分箱，计算每个窗口内的平均信号值

## 依赖环境

| 工具 | 用途 | 加载方式 |
|:-----|:-----|:---------|
| `bedtools` | 生成 genome bins、bedGraph 映射平均 | `module load bedtools` |
| `bigWigAverageOverBed` | BigWig 分箱平均 | `module load deeptools` |
| `bedGraphToBigWig` | bedGraph → BigWig 转换 | UCSC 二进制工具 |
| `sbatch` | 集群任务提交 | SLURM |

## 文件说明

### 核心处理脚本

| 文件 | 功能 | 输入 | 输出 |
|:-----|:-----|:-----|:-----|
| `convert_cpg_to_bedgraph.sh` | 单文件：per-CpG BED → bedGraph | `.bed` / `.bed.gz` | `.bedGraph` |
| `compute_bin_average.sh` | 单文件：BigWig / bedGraph → 分箱平均 | `.bw` / `.bedGraph` + `chrom.sizes` | `.bin{SIZE}.step{STEP}.bp` |

### SLURM 单任务提交

| 文件 | 功能 | 调用核心脚本 |
|:-----|:-----|:-------------|
| `slurm_convert_cpg_to_bedgraph.sh` | 提交单个 CpG BED 转换任务 | 内联 `awk`（等价于 `convert_cpg_to_bedgraph.sh`） |
| `slurm_compute_bin_average.sh` | 提交单个 BigWig / bedGraph 分箱任务 | `compute_bin_average.sh` |

### 批量提交脚本

| 文件 | 功能 |
|:-----|:-----|
| `batch_convert_cpg_to_bedgraph.sh` | 批量提交 CpG BED 转换（遍历目录下所有 `.bed.gz`） |
| `batch_compute_bin_average.sh` | 批量提交 BigWig 分箱（遍历目录下所有 `.bw`） |
| `example_mouse_atac_binning.sh` | 使用示例：小鼠胚胎 ATAC 数据分箱 |

---

## 使用说明

### 1. CpG BED → bedGraph 转换

```bash
# 单文件直接运行
bash convert_cpg_to_bedgraph.sh \
    input.bed.gz \
    output.bedGraph \
    8          # coverage threshold（默认 8）

# 或批量提交到 SLURM
bash batch_convert_cpg_to_bedgraph.sh \
    /path/to/bed_dir \
    /path/to/output_dir \
    8
```

输入 BED 格式要求（Bismark 输出）：
- `$1-$3` : chr, start, end
- `$5`    : coverage（`slurm_convert_cpg_to_bedgraph.sh` 使用）
- `$10`   : coverage（`convert_cpg_to_bedgraph.sh` 使用）
- `$11`   : methylation fraction（0–1）

### 2. BigWig / bedGraph → 分箱平均

```bash
# 单文件直接运行
bash compute_bin_average.sh \
    input.bw \
    /path/to/chrom.sizes \
    /path/to/output_dir \
    1000       # bin size（默认 1000 bp）
    1000       # step size（默认等于 bin size，即无重叠）

# 或批量提交到 SLURM
bash batch_compute_bin_average.sh \
    /path/to/bw_dir \
    /path/to/output_dir \
    /path/to/chrom.sizes \
    1000 \
    1000
```

输入文件支持：`.bw` / `.bigwig` / `.bigWig` / `.bedGraph`

输出格式（4 列 BED）：
```
chr\tstart\tend\tmean_signal
```

---

## 输出文件命名规则

```
${BASENAME}.bin${BIN_SIZE}.step${STEP_SIZE}.bp
```

例如：`H3K27ac.bin1000.step1000.bp`

---

## 注意事项

1. `compute_bin_average.sh` 处理 bedGraph 时会产生 `.expanded` 和 `.sorted` 临时文件，脚本末尾会自动清理
2. `slurm_convert_cpg_to_bedgraph.sh` 与 `convert_cpg_to_bedgraph.sh` 的列号约定略有不同（前者用 `$5` 作 coverage，后者用 `$10`），根据实际 BED 格式选择
3. 所有 SLURM 脚本均提交到 `cpu1,cpu2,cpu_short` 队列，内存 8G，时限 1 小时
