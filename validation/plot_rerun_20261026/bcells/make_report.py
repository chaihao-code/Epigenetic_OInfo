from pathlib import Path
import json,hashlib
from PIL import Image
ROOT=Path(__file__).resolve().parents[3]; OUT=Path(__file__).resolve().parent; CASE=ROOT/'cases/bcells_and_CLL'
runs=json.loads((OUT/'runs.json').read_text())
meta={
'plot_abc_triplet_classes':('ABC三元组：A/B/C每样本nominal与FDR热图、各组均值热图，共9幅','ABC_triplet_classes/plots','[ABC]_*heatmap*','triplet significance TSV + clinical/RNA sample metadata','中间结果+样本元数据；重新聚合计数'),
'plot_abc_pooled_combined':('ABC三类的R/S/pB组均值组合热图','ABC_triplet_classes/plots','ABC_pooled_average_heatmap.*','triplet significance TSV + generated sample labels and class counts','中间结果；独立验证23,800-test BH后聚合'),
'plot_abc_pooled_compact_preview':('ABC组均值紧凑版组合热图','ABC_triplet_classes/plots','ABC_pooled_average_heatmap_compact_preview.*','generated ABC_pooled_average_counts.tsv','中间结果；纯绘图'),
'plot_global_oinfo_and_ABC_triplets':('全系统O-information组间差异与ABC三元组组合图','global_oinfo_and_triplets/plots','global_oinfo_and_ABC_triplets_18x8p5cm.*','order-7 O-info CSV + permutation summary CSV + generated ABC counts','中间结果；重新做40区域Welch检验'),
'plot_de_male_ighv_volcano_clear':('男性RNA method/IGHV校正差异表达火山图','results_rna_de_male_ighv','RNA_IGHV_adjusted_volcano_nominalP_labeled10_6cm_clear.*','rna_method_ighv_adjusted_all_genes.tsv','中间结果；纯绘图'),
'plot_gsea_eight_panel_v2_18x10cm':('8条代表GO通路的GSEA富集曲线组合图','results_gsea_rna_male/male_RNA_method_IGHV_adjusted/eight_panel_GSEA_v2_18x10cm','GSEA_8pathways_18x10cm.*','BP/MF_all_terms.tsv + GMT gene sets + prerank_data.rnk','GSEA中间结果+冻结基因集；重建running ES曲线，无置换'),
'plot_gsea_representative_bars':('代表GO通路NES/FDR条形图','results_gsea_rna_male/male_RNA_method_IGHV_adjusted/representative_pathways','representative_GO_FDR_NES_12x6cm.*','BP/MF/CC_all_terms.tsv','中间结果；纯绘图'),
'plot_gsea_top20_bars':('Top20 GO通路NES/FDR条形图','results_gsea_rna_male/male_RNA_method_IGHV_adjusted/top20_pathways','top20_GO_FDR_NES_12x6cm.*','BP/MF/CC_all_terms.tsv','中间结果；纯绘图'),
'plot_dmr_volcano':('13名CLL供体1 kb差异甲基化火山图','dmr_repro_cemt/output','cemt_R_vs_S_dmr_volcano.*','all_bins_1kb_cemt_g2_g3_limma_results.csv','中间结果；纯绘图，全部2,697,259行'),
'plot_igh_tracks_corrected':('IGH区域甲基化轨迹、组差异和coverage条带','dmr_repro_cemt/output','chr14_105800000_106500000_methylation_tracks_corrected.*','prepared all_bins_1kb_cemt_g2_g3.csv + gencode.v47.annotation.gtf','已处理beta矩阵+GENCODE原始参考注释；不是只读分析结果表，现场计算组均值'),
'validate_male_de_in_females':('38个男性DE基因在4名女性CLL-R中的表达热图与投影图','results_rna_de_female_validation','female_male_signature_concordance.*','length_scaled_counts.tsv + gene_tpm_tximport.tsv + metadata + male DE FDR tables','已导入表达矩阵+中间DE表；混合分析/绘图，重新归一化并拟合男性投影模型'),
}
entries=[]; all_outputs=[]
for run in runs:
 stem=Path(run['script']).stem; title,folder,pattern,inputs,nature=meta[stem]
 paths=sorted(p for p in (CASE/folder).glob(pattern) if p.suffix in ('.pdf','.png','.svg'))
 outputs=[]
 for p in paths:
  if p.suffix=='.png':
   with Image.open(p) as im: im.verify()
  if p.suffix=='.pdf':
   b=p.read_bytes();assert b.startswith(b'%PDF-') and b'%%EOF' in b[-2048:]
  outputs.append(dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
 entry=dict(run,description_zh=title,inputs=inputs,input_nature_zh=nature,status='success' if run['returncode']==0 and len(outputs)>0 else 'failed',figure_count=len([p for p in paths if p.suffix=='.png']),file_count=len(outputs),outputs=outputs)
 entries.append(entry);all_outputs.extend(outputs)
report=dict(case='bcells_and_CLL',date='2026-10-26',entrypoint_count=len(entries),success_count=sum(x['status']=='success' for x in entries),figure_count=sum(x['figure_count'] for x in entries),file_count=len(all_outputs),entries=entries,restoration_manifest='validation/plot_rerun_20261026/bcells/restored_inputs.json',external_input_manifest='validation/plot_rerun_20261026/bcells/external_input_links.json',notes=['All figure files were newly rendered; no historical images copied.','65 missing CSV/TSV/GMT/RNK intermediate files restored from original sibling workspace; 448,053,695 bytes; source/destination/SHA256 recorded.','The prepared beta matrix and GENCODE GTF remain existing external input symlinks; all outputs are in the release.','Initial probing overlapped with input restoration: 8 initial attempts failed; 3 later initial probes already saw restored files. Formal rerun covers all 11 after restoration.','Historical 1,000-permutation triplet/GSEA inputs retain their original meaning; no new permutations or entropy calculations run.'],code_fixes=['plot_abc_pooled_combined.py: added missing os import; layout metadata reads already validated permutation count rather than undefined df.','plot_abc_pooled_compact_preview.py: added missing os import.','plot_igh_tracks_corrected.py: create output directory before saving.'])
(OUT/'report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
lines=['# B-cell/CLL 作图重跑记录','',f'11个入口全部成功，新生成19幅图（19 PNG + 19 PDF）。PNG逐个解码验证；PDF检查文件头/EOF。所有输出位于当前release目录。','', '原release缺少多数计算中间表。本次从原 sibling workspace 复制65个CSV/TSV/GMT/RNK文件（448,053,695 bytes），没有复制历史图片。[输入恢复清单](restored_inputs.json)包含原路径、目标路径和SHA256；[外部输入链接清单](external_input_links.json)标明仍经现有软链接读取的beta矩阵和GTF。','', '| 程序 | 新生成内容 | 幅数/文件数 | 输入性质 | 输出目录 |','|---|---|---:|---|---|']
for e in entries:
 folder=meta[Path(e['script']).stem][1]
 lines.append(f"| `{Path(e['script']).name}` | {e['description_zh']} | {e['figure_count']}/{e['file_count']} | {e['input_nature_zh']} | [{folder}]({CASE/folder}) |")
lines += ['', '## 输入核实','', '纯绘图程序没有读取FASTQ、原始bigWig或重新计算O-information。大部分输入是已有O-information、triplet显著性、DE、GSEA、DMR结果。以下入口需要明确区分：','', '- `plot_igh_tracks_corrected.py`读取全基因组prepared beta矩阵及GENCODE GTF参考注释，计算局部组均值和差异；GTF不属于分析中间文件。两者是现有外部软链接。','- `validate_male_de_in_females.py`从已导入counts/TPM重做归一化和男性模型拟合再画图，属于分析与作图混合入口。','- `plot_global_oinfo_and_ABC_triplets.py`从已有order-7数值重做40区域Welch检验；`plot_abc_pooled_combined.py`重新验证全局BH，均不运行基因组计算。','- GSEA 8面板从已有排名和GMT重建running ES曲线，并断言与已有ES一致，不运行GSEA置换。','', '## 运行证据与修复','', '最终逐项命令、退出码、时长、结果路径与SHA256在[report.json](report.json)；[run.py](run.py)可在release中复跑，先并行8项，再按依赖顺序画ABC组合、紧凑版和global组合图。日志在[logs](logs)。','', '实际执行发现并修复：ABC两个入口缺少`os`导入；ABC pooled布局元数据引用未定义`df`；IGH入口未自动创建输出目录。科学定义未修改。','', '初始探测与输入恢复有时间重叠：8项初始失败（日志保留），其余3项探测时已可见恢复后的文件。正式重跑是在恢复后完成，并对失败修复后再次执行。','', '历史triplet和GSEA输入中的1,000次置换仍保留原义；本次只重画和必要表层统计，没有将其标成新的10,000次计算。']
(OUT/'REPORT.md').write_text('\n'.join(lines)+'\n')
print(report['success_count'],report['figure_count'],report['file_count'])
