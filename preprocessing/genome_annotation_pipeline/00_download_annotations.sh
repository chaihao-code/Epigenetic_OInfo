#!/usr/bin/env bash
# 00 下载原始注释文件（如已存在则跳过）
# 用法: GENOME=hg38 bash 00_download_annotations.sh / GENOME=mm10 bash ...
set -euo pipefail
source "$(dirname "$0")/config.sh"

dl() {  # dl <url> <out>
    if [ -s "$2" ]; then echo "[skip] $2"; else
        echo "[GET] $1"; wget -O "$2" "$1"
    fi
}

if [ "$GENOME" = "hg38" ]; then
    # GENCODE v47
    dl https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz "$GTF.gz"
    gunzip -f "$GTF.gz"
    # UCSC 注释
    dl https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz "$ANN/hg38/hg38.cpgIslandExt.txt.gz"
    dl https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz "$ANN/hg38/hg38.RepeatMasker.rmsk.gz"
    dl https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes "$CHROM_SIZES"
    # ENCODE 增强子/CTCF（hg38）
    dl https://www.encodeproject.org/files/ENCFF936NSJ/@@download/ENCFF936NSJ.bed.gz "$ENH_PROX.gz"
    dl https://www.encodeproject.org/files/ENCFF535MKS/@@download/ENCFF535MKS.bed.gz "$ENH_DIST.gz"
    dl https://www.encodeproject.org/files/ENCFF996ZSG/@@download/ENCFF996ZSG.bed.gz "$CTCF.gz"
    gunzip -f "$ENH_PROX.gz" "$ENH_DIST.gz" "$CTCF.gz"
else
    # GENCODE mm10（exon/intron）
    dl https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz "$GTF.gz"
    gunzip -f "$GTF.gz"
    # EPDnew v6 启动子（小鼠 promoter 的 TSS 来源）
    dl https://epd.expasy.org/ftp/EPDnew_mouse/mm10/EPDnew_mm10.tsv "$EPD"
    # UCSC 注释
    dl https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/cpgIslandExt.txt.gz "$ANN/mm10.cpgIslandExt.txt.gz"
    dl https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/rmsk.txt.gz "$ANN/mm10.RepeatMasker.rmsk.tsv.gz"
    dl https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes "$CHROM_SIZES"
fi
echo "[done] $GENOME annotations in $ANN_DIR"
