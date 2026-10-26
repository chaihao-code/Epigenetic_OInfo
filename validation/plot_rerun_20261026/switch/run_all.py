"""Rerun every H3K36me3 figure entrypoint from retained intermediate tables."""
import os,sys,json,time,subprocess
from pathlib import Path
root=Path(__file__).resolve().parents[3]
base=root/'cases/H3K36me3_developmental_switch'
report=root/'validation/plot_rerun_20261026/switch'
(report/'logs').mkdir(parents=True,exist_ok=True)
for sub in ['batch_effect_figures/plots','redundancy_core/plots','redundancy_core/full_signal_occurrence/plots','synergistic_triplets/plots']:(base/sub).mkdir(parents=True,exist_ok=True)
env=os.environ.copy(); env.update(EPINFO_ROOT=str(root),MPLBACKEND='Agg',MPLCONFIGDIR='/tmp/epinfo-switch-mpl')
labels={
'plot_batch_effect':'16种信号来源替换组合的同阶段标记参与比例',
'plot_support_without_h3k27me3_replacement':'去除H3K27me3替换后的8组合支持度分布（两个FDR阈值）',
'plot_triplet_support_distribution':'16组合的triplet支持度分布（两个FDR阈值）',
'plot_combined_trajectories':'triplet与冗余核心标记轨迹并排合图',
'plot_core_occurrence_frequency_e11.5':'历史发育后期冗余核心Top25频率、累计线和点阵',
'plot_core_occurrence_frequency_preimp':'历史植入前冗余核心Top25频率、累计线和点阵',
'plot_full_signal_combination_occurrence':'全信号triplet/冗余核心Top25频率（4单图和2×2合图）',
'plot_redundancy_core_trajectories':'8标记在冗余核心中的发育轨迹',
'plot_synergistic_triplet_trajectories':'8标记在协同triplet中的发育轨迹',
'plot_triplet_combination_occurrence_with_H3K9ac_H3K4me2':'含H3K9ac/H3K4me2全信号A/B/C组合频次点阵（仅符号筛选）',
'plot_mark_occurrence_per_tissue_by_stage':'标记参与次数：植入前、发育后期及合图',
'plot_original_development_heatmap':'原始8信号14阶段标记参与比例热图',
'plot_triplet_ABC_selected_regions_post_E11p5':'8信号选定区域A/B/C相对出现次数轨迹',
'plot_triplet_ABC_stacked_post_E11p5':'8信号发育后期A/B/C三层热图',
'plot_triplet_ABC_stacked_preimplantation':'8信号植入前A/B/C三层热图',
'plot_triplet_A_B_C_post_E11p5':'8信号4区域A/B/C相对出现次数轨迹'}
records=[]
for p in sorted(base.rglob('plot*.py')):
 if p.name in ('plot_batch_effect.py','plot_redundancy_core_trajectories.py','plot_synergistic_triplet_trajectories.py','plot_full_signal_combination_occurrence.py'):
  variants=[['--plots-only']]
 elif p.name=='plot_triplet_combination_occurrence_with_H3K9ac_H3K4me2.py':
  variants=[[cls,'--plots-only'] for cls in 'ABC']
 else:variants=[[]]
 cmds=[];logs=[];codes=[];outputs=[];errors=[]
 for args in variants:
  cmd=[sys.executable,str(p.relative_to(root)),*args];cmds.append(' '.join(cmd));stamp=time.time()
  result=subprocess.run(cmd,cwd=root,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
  suffix='_'+args[0] if args and args[0] in 'ABC' else ''
  log=report/'logs'/f'{p.stem}{suffix}.log';log.write_text(result.stdout);logs.append(str(log.relative_to(root)))
  changed=[str(x.relative_to(root)) for x in base.rglob('*') if x.is_file() and x.suffix in ('.png','.pdf','.svg') and x.stat().st_mtime>=stamp]
  outputs+=changed;codes.append(result.returncode)
  if result.returncode:errors.append(result.stdout.strip().splitlines()[-1])
 name=p.stem
 if 'core_occurrence_frequency_' in name:inputs=['historical all_supreme_cores_filtered.csv (missing)']
 elif name=='plot_batch_effect':inputs=['batch_effect_figures/same_stage_mark_summary.tsv','batch_effect_figures/data/late_common_tissue_means.tsv','batch_effect_figures/stability_by_stage*.tsv','batch_effect_figures/manifest_16_fdr.json']
 elif name=='plot_support_without_h3k27me3_replacement':inputs=['batch_effect_figures/data/all_16variant_triplets_bh56.tsv.gz','batch_effect_figures/manifest_16_fdr.json']
 elif name=='plot_triplet_support_distribution':inputs=['batch_effect_figures/fdr_lt_0p1/triplet_support_across_variants.tsv','batch_effect_figures/fdr_lt_0p1/threshold_comparison.tsv']
 elif name=='plot_combined_trajectories':inputs=['redundancy_core/results/*_mark_trajectories.tsv','redundancy_core/inputs/analysis_manifest.json']
 elif name in ('plot_redundancy_core_trajectories','plot_synergistic_triplet_trajectories'):inputs=['redundancy_core/results/*_mark_counts_by_sample.tsv','redundancy_core/inputs/analysis_manifest.json']
 elif name=='plot_full_signal_combination_occurrence':inputs=['redundancy_core/full_signal_occurrence/data/*_combination_occurrence_*_top25.tsv','redundancy_core/full_signal_occurrence/data/occurrence_summary.tsv','redundancy_core/full_signal_occurrence/inputs/analysis_manifest.json']
 elif name=='plot_triplet_combination_occurrence_with_H3K9ac_H3K4me2':inputs=['redundancy_core/full_signal_occurrence/data/synergistic_triplet_combination_occurrence_post_E11p5_all.tsv']
 elif name=='plot_mark_occurrence_per_tissue_by_stage':inputs=['synergistic_triplets/data/original_development_mark_participation*.tsv','synergistic_triplets/inputs/source_manifest.json']
 else:inputs=['batch_effect_figures/data/*triplets_bh56.tsv.gz','synergistic_triplets/inputs/source_manifest.json','synergistic_triplets/inputs/preimplantation_additional_stages_bh56.tsv.gz (early panels only)','synergistic_triplets/data/*summary.tsv (cross-checks)']
 row=dict(script=str(p.relative_to(root)),title=labels[name],status='success' if all(c==0 for c in codes) and outputs else 'blocked',command=' && '.join(cmds),commands=cmds,exit_codes=codes,outputs=sorted(set(outputs)),figure_count=sum(x.endswith('.png') for x in set(outputs)),logs=logs,inputs=inputs,input_type='intermediate results + metadata',inputs_only_intermediate=True,reads_raw_signal=False,error='; '.join(errors),notes='No O-information/permutation/core recomputation. Existing results retain original scientific settings.')
 records.append(row);(report/'results.json').write_text(json.dumps(records,ensure_ascii=False,indent=2))
 print(p.name,row['status'],row['figure_count'],row['error'],flush=True)
lines=['# H3K36me3作图复跑记录','',f'入口{len(records)}个，成功{sum(r["status"]=="success" for r in records)}个，缺输入{sum(r["status"]!="success" for r in records)}个。新画{sum(r["figure_count"] for r in records)}张图（每张PNG和PDF）。','', '全部入口均读取中间结果与元数据；未读取bigWig、FASTQ或原始信号，也未重跑置换或冗余核心分析。部分原入口会从已有统计全表重汇总并核对BH值。','', '复跑：`python validation/plot_rerun_20261026/switch/run_all.py`。日志和逐文件路径见`results.json`。','', '|入口|画了什么|状态|张数|输出目录|','|---|---|---|---:|---|']
for r in records:
 dirs=sorted(set(str(Path(x).parent) for x in r['outputs']))
 lines.append(f'|{Path(r["script"]).name}|{r["title"]}|{r["status"]}|{r["figure_count"]}|'+('<br>'.join(dirs) or r['error'])+'|')
lines+=['','修复记录：','', '- 从provenance指定的本机upstream备份和原工作目录恢复6个缺失的统计中间文件。SHA256、源路径和大小见`restored_inputs.json`。这些文件现已复制到release中，复跑不依赖备份目录。','- 为原脚本输出创建plots目录；trajectory和full-signal入口增加`--plots-only`，从已有汇总表重画。','- 旧triplet组合脚本默认D没有分类定义，原本退出0但不出图。本次保留ABC分类规则、显式跑A/B/C并读取现有符号筛选计数；默认修为合法A，非法类会报错。','- 植入前堆叠热图的Zygote标签超出画布；调整同一绘图函数左右布局后通过原有文字边界和单元格尺寸断言。后期共享同一布局函数。','- 两个historical core Top25入口缺原始的`all_supreme_cores_filtered.csv`；未用其他核心结果替代。','- `trajectory_utils.py`是无main的共享绘图模块，由两个trajectory入口调用，不单独计入口。','', '输入表可能记录历史统计设置；本次仅重画，不将旧结果标记为重新完成10000次置换。']
(report/'REPORT.md').write_text('\n'.join(lines)+'\n')
