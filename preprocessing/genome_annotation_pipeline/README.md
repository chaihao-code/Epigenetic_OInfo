# Genome annotation → 1-kb bin masks → purified region masks

从原始注释文件到最终 `*_purified_bins.tsv` 的完整流程。同时支持 hg38 与 mm10。

## 流程

```
00_download_annotations.sh   下载 GTF / UCSC 注释 / ENCODE / EPDnew / chrom.sizes
01_make_bins.sh              bedtools makewindows → bins_1kb.bed
02_extract_gene_features.sh  GTF → exons/promoters/transcripts/introns + cpg/TE beds
03_ctcf_split.sh             enhancer ∩ CTCF → ctcf_bound / not_ctcf_bound
04_bed_to_bins.py            BED → 0/1 bin mask（通用）
05_run_bed_to_bins.sh        对每个特征调用 04，生成 *_bins.tsv
06_purify_masks.py           按优先级重建 *_purified_bins.tsv（互斥）
```

按编号顺序执行；`config.sh` 里统一配置路径。

## 输入文件与版本

| 文件 | 来源 | 用于 |
|---|---|---|
| gencode.v47.annotation.gtf [58] | GENCODE release 47 (hg38) | exon / promoter / intron |
| gencode.vM25.annotation.gtf [57] | GENCODE release M25 (mm10) | exon / intron（小鼠） |
| mm10.EPDnewv6.promoter.tsv | EPDnew v6 | promoter（小鼠） |
| hg38.cpgIslandExt / mm10.cpgIslandExt | UCSC | CpG island |
| hg38.RepeatMasker.rmsk / mm10.RepeatMasker.rmsk | UCSC/RepeatMasker | LINE / SINE / LTR |
| ENCFF936NSJ | ENCODE (proximal enhancer, hg38) | proximal enhancer ±CTCF |
| ENCFF535MKS | ENCODE (distal enhancer, hg38) | distal enhancer ±CTCF |
| ENCFF996ZSG | ENCODE (CTCF, hg38) | CTCF 划分 |
| *.chrom.sizes | UCSC / ENCODE | bin 生成 |

## 关键定义（重要）

- **promoter = 包含 TSS 的 1 个碱基对的区间**（`02_extract_gene_features.sh`），
  即"含 TSS 的 1-kb bin"。**不是** TSS±2kb——脚本注释中"TSS ± 2kb"是笔误。
  小鼠 promoter 用 EPDnew v6（TSS 附近 60 bp 区间），同样是窄定义。
- **层级（06_purify_masks.py）**：promoter > prox enh CTCF+ > prox enh CTCF− >
  distal enh CTCF+ > distal enh CTCF− > exon > intron，每层剔除所有上层，
  各层互斥。
- bin→基因映射（下游 GSEA 用）是 **TSS±2kb** 窗口，与 promoter 定义不同尺度，
  因此 exon/intron 区间中 TSS 近端（1–2 kb）的 bin 仍会被映射到基因，
  方法/局限性中需披露这一点。

## 已知问题 / 历史

- `06_purify_masks.py` 取代旧的 `code/genome_feature_to_bin/purify_genomic_regions.py`
  （旧版 distal 层不剔除 prox enh 重合、exon/intron 只剔除 promoter 且彼此不互斥）。
- 旧版 purified mask 的备份：同目录 `*_purified_bins.tsv.bak`。
- mm10 的 exon/intron 注释最初在集群上用 mm10 GTF 生成，本目录脚本以
  gencode.vM25 重生成，bin 数可能与历史文件略有出入，使用前应比对。

## 参考文献编号

- [57] GENCODE mouse release M25 (gencode.vM25.annotation.gtf, mm10)
- [58] GENCODE human release 47 (gencode.v47.annotation.gtf, hg38)

编号与论文参考文献列表一致。
