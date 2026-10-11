#!/usr/bin/env bash
# 统一路径配置，被 00–06 各脚本 source
set -euo pipefail

REPO=/media/dell/data5/sequence_feature/sequence_entropy
ANN=$REPO/data/genome_annotation          # 原始注释文件
FEAT=$REPO/result/genome_features         # 输出 *_bins.tsv / *_purified_bins.tsv

GENOME=${GENOME:-hg38}                    # hg38 | mm10
BIN_SIZE=1000

if [ "$GENOME" = "hg38" ]; then
    GTF=$ANN/hg38/gencode.v47.annotation.gtf
    CPG=$ANN/hg38/hg38.cpgIslandExt.tsv
    RMSK=$ANN/hg38/hg38.RepeatMasker.rmsk
    CHROM_SIZES=$ANN/hg38/hg38.chrom.sizes
    ANN_DIR=$ANN/hg38/hg38_annotations
    # ENCODE (hg38 only)
    ENH_PROX=$ANN_DIR/ENCFF936NSJ_hg38_proximal_enhancer.bed
    ENH_DIST=$ANN_DIR/ENCFF535MKS_hg38_distal_enhancer.bed
    CTCF=$ANN/ENCFF996ZSG.bed
else
    GTF=$ANN/mm10/gencode.vM25.annotation.gtf   # 小鼠 exon/intron
    EPD=$ANN/mm10.EPDnewv6.promoter.tsv         # 小鼠 promoter (TSS±30bp)
    CPG=$ANN/mm10.cpgIslandExt.tsv
    RMSK=$ANN/mm10.RepeatMasker.rmsk.tsv
    CHROM_SIZES=$ANN/mm10.chrom.sizes
    ANN_DIR=$ANN/mm10/mm10_annotations
fi

BINS=$ANN_DIR/bins_${BIN_SIZE}bp.bed
OUT_DIR=$FEAT/$GENOME
mkdir -p "$ANN_DIR" "$OUT_DIR"
