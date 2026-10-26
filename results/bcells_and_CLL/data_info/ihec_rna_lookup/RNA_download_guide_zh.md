# IHEC / CEMT RNA 下载清单

核对日期：2026-10-12。只读取元数据、目录索引和一个定量文件的表头，没有下载完整 RNA 数据。

**13 个比较样本以及另外 4 个未分组样本，在最新版 EpiATLAS 中均有公开基因表达定量文件。建议优先下载 `genes.results` 和 `stats.json`。**

## 文件选择

- `genes.results`：基因级 RSEM 定量。已核实表头包含 `gene_id`、`expected_count`、`TPM`、`FPKM` 及后验估计列；`expected_count` 为估计计数，不是逐基因原始整数计数。后续可按定量方法选择适合的统计处理。
- `isoforms.results`：可选，转录本级定量。
- `stats.json` / `report.pdf`：RNA QC，建议保留。
- `*.md5`：官方校验和。
- BigWig 是覆盖轨迹；基因差异分析优先使用定量表。此处 portal 文件名 mRNA-Seq 与 CEMT 原始样本 total RNA 记录属于不同层次标签，不据此认定 polyA 富集协议。

## 样本和下载目录

|样本|当前分组|EpiATLAS ID|RNA 分析方式|定量文件目录|
|---|---|---|---|---|
|CEMT0004|CLL-R|IHECRE00000241.3|链特异|[下载目录](https://ihec-epigenomes.org/epiatlas/index/base_data/IHECRE00000241/mRNA-Seq/)|
|CEMT0005|CLL-R|IHECRE00000246.3|链特异|[下载目录](https://ihec-epigenomes.org/epiatlas/index/base_data/IHECRE00000246/mRNA-Seq/)|
|CEMT0006|CLL-S|IHECRE00000229.3|链特异|[下载目录](https://ihec-epigenomes.org/epiatlas/index/base_data/IHECRE00000229/mRNA-Seq/)|
|CEMT0025|CLL-R|IHECRE00000247.3|链特异|[下载目录](https://ihec-epigenomes.org/epiatlas/index/base_data/IHECRE00000247/mRNA-Seq/)|
|CEMT0026|CLL-R|IHECRE00000239.3|链特异|[下载目录](https://ihec-epigenomes.org/epiatlas/index/base_data/IHECRE00000239/mRNA-Seq/)|
|CEMT0027|CLL-R|IHECRE00000234.3|链特异|[下载目录](https://ihec-epigenomes.org/epiatlas/index/base_data/IHECRE00000234/mRNA-Seq/)|
|CEMT0028|CLL-R|IHECRE00000232.3|链特异|[下载目录](https://ihec-epigenomes.org/epiatlas/index/base_data/IHECRE00000232/mRNA-Seq/)|
|CEMT0029|CLL-R|IHECRE00000244.3|链特异|[下载目录](https://ihec-epigenomes.org/epiatlas/index/base_data/IHECRE00000244/mRNA-Seq/)|
|CEMT0030|CLL-R|IHECRE00000226.3|链特异|[下载目录](https://ihec-epigenomes.org/epiatlas/index/base_data/IHECRE00000226/mRNA-Seq/)|
|CEMT0094|CLL-S|IHECRE00003687.3|非链特异|[下载目录](https://ihec-epigenomes.org/epiatlas/index/base_data/IHECRE00003687/mRNA-Seq/)|
|CEMT0095|CLL-S|IHECRE00003660.4|非链特异|[下载目录](https://ihec-epigenomes.org/epiatlas/index/base_data/IHECRE00003660/mRNA-Seq/)|
|CEMT0096|CLL-S|IHECRE00003688.3|非链特异|[下载目录](https://ihec-epigenomes.org/epiatlas/index/base_data/IHECRE00003688/mRNA-Seq/)|
|CEMT0097|CLL-S|IHECRE00003689.3|链特异|[下载目录](https://ihec-epigenomes.org/epiatlas/index/base_data/IHECRE00003689/mRNA-Seq/)|
|CEMT0129|unannotated|IHECRE00003692.4|链特异|[下载目录](https://ihec-epigenomes.org/epiatlas/index/base_data/IHECRE00003692/mRNA-Seq/)|
|CEMT0134|unannotated|IHECRE00003821.3|链特异|[下载目录](https://ihec-epigenomes.org/epiatlas/index/base_data/IHECRE00003821/mRNA-Seq/)|
|CEMT0135|unannotated|IHECRE00003828.3|链特异|[下载目录](https://ihec-epigenomes.org/epiatlas/index/base_data/IHECRE00003828/mRNA-Seq/)|
|CEMT0138|unannotated|IHECRE00003836.3|链特异|[下载目录](https://ihec-epigenomes.org/epiatlas/index/base_data/IHECRE00003836/mRNA-Seq/)|

## 在另一台机器下载

复制本目录的 `download_rna_cll13.sh` 到有网络的机器，然后：

```bash
bash download_rna_cll13.sh ./bcells_rna_cll13
```

脚本下载 13 个基因级定量文件及 13 个 stats.json，输出按 CEMT 样本号命名。只下载上述 26 个文件，支持断点续传。脚本做过 bash 语法检查，未实际运行完整下载。

如果希望保留原始文件名，可复制 `rna_genes_cll13.urls.txt` 并运行：

```bash
mkdir -p bcells_rna_cll13
wget --continue --timeout=60 --tries=5 -P bcells_rna_cll13 -i rna_genes_cll13.urls.txt
```

可选清单：`rna_isoforms_cll13.urls.txt` 是转录本文件；`rna_qc_cll13.urls.txt` 是 stats、PDF 和 md5。`all17` 版本包括 4 个未分组样本。`rna_availability.tsv` 包含 CEMT、原始 RNA 样本 ERS accession、EpiRR ID、输入 EGAR run、uuid、供者年龄性别以及各类直接下载链接。

## 混杂因素

官方实验元数据中 `analyzed_as_stranded`：CEMT0094/0095/0096 为 False；其余 10 个比较样本为 True。说明这三个 S 样本在统一流程中按非链特异方式分析；这是明确的分析方式差异，但不单独证明所有甲基化/组蛋白差异来自批次。

S 组 5 个供者都是男性；R 组为 4 男/4 女。年龄和其它 QC 也应在后续 RNA 分析中核对。新增元数据未修改原有分组文件和已完成的甲基化结果。

旧版 IHEC 2020-10 hg38 的 CEMT 命名 mRNA 轨迹未列出 CEMT0094/0095/0096，不能据此说它们没有 RNA。最新版使用 IHECRE 命名，可以找到三者的轨迹与 genes.results。smRNA-Seq 是小 RNA 数据，不能用来替代这些蛋白编码基因的表达定量。

## 原始测序文件

IHEC 最新索引将这些样本链接至 [EGAD00001003963](https://ega-archive.org/datasets/EGAD00001003963)，所属研究 [EGAS00001000552](https://ega-archive.org/studies/EGAS00001000552)。这是包含多种实验的累计 BAM/FASTQ 数据集，不应下载全套后才筛选 RNA。EGA 数据集受控；官方政策要求先完成机构间数据访问协议，之后再向 DAC 申请。获得权限后按 RNA 样本、实验和输入 EGAR run 筛选，使用 EGA 下载工具。公开的 EpiATLAS genes.results 无需这一原始测序访问流程。

## 核对来源

- [EpiATLAS 下载门户](https://ihec-epigenomes.org/EpiATLAS/data/)
- [Observed 数据 README](https://ihec-epigenomes.org/epiatlas/index/readme/observed_data_readme.txt)
- [EpiATLAS 实验元数据](https://ihec-epigenomes.org/epiatlas/index/readme/EpiATLAS_experiment_metadata.csv)
- [RNA QC 汇总](https://ihec-epigenomes.org/epiatlas/index/readme/EpiATLAS_rnaseq_qc_summary.csv)
- [CEMT 供者及 RNA 样本元数据](https://www.bcgsc.ca/downloads/edcc/data/CEMT/temp/cemt_samples.html)
- IHEC 公开只读 API：`/cgi-bin/api/getDataHub.py?donor_id=...` 将 CEMT 供者映射到 EpiRR 实验；`/api/datahub/download?build=EpiATLAS_v1&assembly=4` 交叉核对最新 BigWig 索引。各样本的定量目录均已请求验证存在，并与实验 uuid 匹配。
