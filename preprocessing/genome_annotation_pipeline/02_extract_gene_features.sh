#!/usr/bin/env bash
# 02 从 GTF/EPDnew/UCSC 提取基因特征与序列特征 BED
#
# 注意：promoter 这里是"包含 TSS 的 1-bp 区间"（即含 TSS 的 1-kb bin），
# 不是 TSS±2kb。该定义与下游 GSEA 的 TSS±2kb 基因映射窗口不同尺度，
# exon/intron 区间会保留 TSS 近端 1–2 kb 的 bin（见 README"关键定义"）。
set -euo pipefail
source "$(dirname "$0")/config.sh"

if [ "$GENOME" = "hg38" ]; then
    # exons
    awk '$3=="exon" {
        match($9, /gene_id "([^"]+)"/, gid);
        match($9, /transcript_id "([^"]+)"/, tid);
        print $1"\t"$4-1"\t"$5"\t"gid[1]"\t"tid[1]"\t"$7
    }' "$GTF" > "$ANN_DIR/exons.bed"

    # promoters = TSS 位点（1 bp，BED 半开区间）
    awk '$3=="transcript" {
        match($9, /gene_id "([^"]+)"/, gid);
        match($9, /transcript_id "([^"]+)"/, tid);
        if($7=="+") tss=$4; else tss=$5;
        print $1"\t"tss-1"\t"tss"\t"gid[1]"\t"tid[1]"\t"$7
    }' "$GTF" > "$ANN_DIR/promoters.bed"

    # transcripts → introns = transcript − exon
    awk '$3=="transcript" {
        match($9, /gene_id "([^"]+)"/, gid);
        match($9, /transcript_id "([^"]+)"/, tid);
        print $1"\t"$4-1"\t"$5"\t"gid[1]"\t"tid[1]"\t"$7
    }' "$GTF" > "$ANN_DIR/transcripts.bed"
    bedtools subtract -a "$ANN_DIR/transcripts.bed" -b "$ANN_DIR/exons.bed" \
        > "$ANN_DIR/introns.bed"
else
    # mm10 promoter: EPDnew v6，TSS 位于 60-bp 窗口的内侧边缘
    # (+ 链 TSS=chromEnd，- 链 TSS=chromStart，BED 约定)
    awk 'NR>1 {
        if($6=="+"){s=$3-1; e=$3} else {s=$2; e=$2+1}
        print $1"\t"s"\t"e"\t"$4"\t.\t"$6
    }' "$EPD" > "$ANN_DIR/promoters.bed"

    # mm10 exon/intron 与 hg38 相同方式从 GTF 提取
    awk '$3=="exon" {
        match($9, /gene_id "([^"]+)"/, gid);
        match($9, /transcript_id "([^"]+)"/, tid);
        print $1"\t"$4-1"\t"$5"\t"gid[1]"\t"tid[1]"\t"$7
    }' "$GTF" > "$ANN_DIR/exons.bed"
    awk '$3=="transcript" {
        match($9, /gene_id "([^"]+)"/, gid);
        match($9, /transcript_id "([^"]+)"/, tid);
        print $1"\t"$4-1"\t"$5"\t"gid[1]"\t"tid[1]"\t"$7
    }' "$GTF" > "$ANN_DIR/transcripts.bed"
    bedtools subtract -a "$ANN_DIR/transcripts.bed" -b "$ANN_DIR/exons.bed" \
        > "$ANN_DIR/introns.bed"
fi

# CpG island（UCSC cpgIslandExt 前三列）
tail -n +2 "$CPG" | awk 'BEGIN{OFS="\t"} {print $2,$3,$4}' > "$ANN_DIR/cpg_island.bed"

# 重复序列 LINE/SINE/LTR（RepeatMasker repClass 列）
if [ "$GENOME" = "hg38" ]; then RCOL=12; C1=6; C2=7; C3=8; SKIP=1; else RCOL=12; C1=6; C2=7; C3=8; SKIP=1; fi
for cls in LINE SINE LTR; do
    awk -v C="$cls" -v r="$RCOL" -v a="$C1" -v b="$C2" -v c="$C3" \
        'NR>'"$SKIP"' && $r==C {print $a,$b,$c}' OFS='\t' "$RMSK" \
        > "$ANN_DIR/${cls}.bed"
done
echo "[done] gene/sequence feature BEDs in $ANN_DIR"
