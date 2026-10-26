#!/usr/bin/env python3
"""Re-render the two thyroid GSEA workflows solely from saved reports/rankings/GMTs.
The historical PDF names choose terms; historical PDF contents are not read.
No regression, enrichment permutations, GTF parsing, or signal extraction runs.
"""
import os
os.environ.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',MPLBACKEND='Agg')
os.environ.setdefault('MPLCONFIGDIR','/tmp/thyroid-gsea-rerun-mpl')
from pathlib import Path
import argparse,concurrent.futures,hashlib,json,time,textwrap,traceback
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).resolve().parent

def sha256(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for part in iter(lambda:f.read(8*1024*1024),b''):h.update(part)
 return h.hexdigest()

def run_group(group):
 start=time.time();folder=ROOT/group['input_folder'];gid=str(folder.relative_to(ROOT/'cases/aging/thyroid/local_O_info_regression/results')).replace('/','__')
 record=dict(group,status='success',figures=[],seconds=0.)
 try:
  if not group['historical_pdf_names']:
   record['reason']='This original region/library had no historical enrichment curves.'
  else:
   report=pd.read_csv(folder/'gseapy.gene_set.prerank.report.csv',usecols=['Term','ES','NES','NOM p-val','FDR q-val'])
   mapping={row['Term'].replace('/','-').replace(':','_')+'.pdf':row for row in report.to_dict('records')}
   assert len(mapping)==len(report),'Nonunique normalized term names'
   assert set(group['historical_pdf_names'])<=mapping.keys(),'Historical term not in report'
   selected={mapping[name]['Term'] for name in group['historical_pdf_names']}
   sets={}
   with (folder/'gene_sets.gmt').open() as f:
    for line in f:
     fields=line.rstrip('\n').split('\t')
     if fields[0] in selected:sets[fields[0]]=fields[2:]
   assert sets.keys()==selected,'Selected gene set absent from frozen GMT'
   rank=pd.read_csv(folder/'prerank_data.rnk',sep='\t',header=None,names=['gene','score'],dtype={'gene':str,'score':np.float64})
   assert rank.gene.is_unique and np.isfinite(rank.score.to_numpy()).all()
   genes=pd.Index(rank.gene);scores=rank.score.to_numpy();magnitude=np.abs(scores);x=np.arange(1,len(rank)+1)
   assert np.all(scores[:-1]>=scores[1:]),'Preserve descending saved order'
   target=folder/'prerank';target.mkdir(exist_ok=True)
   for name in group['historical_pdf_names']:
    row=mapping[name];term=row['Term'];positions=genes.get_indexer(sets[term]);positions=positions[positions>=0]
    hit=np.zeros(len(rank),dtype=bool);hit[positions]=True
    assert 0<int(hit.sum())<len(hit)
    hitweight=hit*magnitude;normalizer=hitweight.sum();assert normalizer>0
    running=np.cumsum(hitweight/normalizer-(~hit)/(~hit).sum());es=float(running[np.argmax(np.abs(running))]);error=abs(es-float(row['ES']))
    # Reports preserve seven significant digits and ranked scores were exported as text.
    assert error<=2e-6,f'ES mismatch {term}: saved={row["ES"]} recomputed={es}'
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    fig=plt.figure(figsize=(7.2,5.5))
    grid=fig.add_gridspec(3,1,height_ratios=[3.,.55,1.4],hspace=.10,left=.13,right=.96,bottom=.14,top=.82)
    a=fig.add_subplot(grid[0]);b=fig.add_subplot(grid[1],sharex=a);c=fig.add_subplot(grid[2],sharex=a)
    a.plot(x,running,color='#3b9950',lw=.9,rasterized=True)
    a.axhline(0,color='gray',lw=.5);a.set_ylabel('Enrichment score (ES)')
    a.text(.985,.96,f'NES = {float(row["NES"]):.3g}\nP = {float(row["NOM p-val"]):.3g}\nFDR q = {float(row["FDR q-val"]):.3g}',ha='right',va='top',transform=a.transAxes,fontsize=9,bbox={'facecolor':'white','alpha':.8,'edgecolor':'none'})
    # Full hit positions and full saved ranking are rendered; no locus subsampling.
    b.vlines(np.flatnonzero(hit)+1,0,1,color='black',lw=.3,rasterized=True)
    b.set_yticks([]);b.set_ylabel('Hits',fontsize=8);b.spines['left'].set_visible(False);b.spines['bottom'].set_visible(False)
    c.plot(x,scores,color='#66859a',lw=.65,rasterized=True);c.axhline(0,color='gray',lw=.5)
    c.set_ylabel('Ranked metric');c.set_xlabel('Rank in ordered gene/bin list')
    c.set_xlim(1,len(rank));a.tick_params(labelbottom=False);b.tick_params(labelbottom=False)
    fig.suptitle('\n'.join(textwrap.wrap(term,83)),fontsize=11,y=.97)
    region=folder.parent.name if group['kind']=='gsea_pathway' else folder.parent.name
    library=folder.name if group['kind']=='gsea_pathway' else 'signal-change gene sets'
    fig.text(.13,.855,region+'\n'+library,fontsize=8)
    fig.text(.13,.025,'Reconstructed from saved ranks and gene sets; NES/P/FDR retained from the original report.',fontsize=7)
    outputs=[]
    for ext in ['pdf','png']:
     p=target/(Path(name).stem+'.'+ext);fig.savefig(p,dpi=180)
     if ext=='png':
      with Image.open(p) as im:im.verify()
     else:
      content=p.read_bytes();assert content.startswith(b'%PDF-') and b'%%EOF' in content[-2048:]
     outputs.append(dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=sha256(p)))
    plt.close(fig)
    record['figures'].append(dict(term=term,historical_pdf_name=name,n_ranked=len(rank),n_hits=int(hit.sum()),saved_ES=float(row['ES']),recomputed_ES=es,absolute_ES_error=error,NES=float(row['NES']),nominal_P=float(row['NOM p-val']),FDR_q=float(row['FDR q-val']),outputs=outputs))
 except Exception:
  record['status']='failed';record['error']=traceback.format_exc()
 finally:
  record['seconds']=round(time.time()-start,2)
  (OUT/'groups_results').mkdir(exist_ok=True)
  (OUT/'groups_results'/(gid+'.json')).write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n')
 return record

def make_report(records):
 entries=[]
 for kind,script,title,nature in [
  ('gsea_pathway','02_gsea_pathway.py','甲状腺局部O-information年龄效应：GO BP/Hallmark GSEA曲线','原入口读取回归表和GENCODE GTF，进行TSS注释、基因排名、获取基因集及GSEA置换；不是纯作图。此次仅读冻结report/RNK/GMT。'),
  ('gsea_marginal','03_gsea_marginal.py','甲状腺局部O-information年龄效应：各信号up/down/flat GSEA曲线','原入口读取回归表及extracted_signals.tsv，重新按年龄相关性分类bins、构建基因集并做GSEA置换；不是纯作图。此次仅读冻结report/RNK/GMT。')]:
  groups=[r for r in records if r['kind']==kind];figs=[f for g in groups for f in g['figures']];expected_groups=[g for g in json.loads((OUT/'groups.json').read_text()) if g['kind']==kind];expected=sum(len(g['historical_pdf_names']) for g in expected_groups);outputs=[o for f in figs for o in f['outputs']]
  complete_groups={g['input_folder'] for g in groups}=={g['input_folder'] for g in expected_groups}
  entries.append(dict(script='cases/aging/thyroid/local_O_info_regression/scripts/'+script,description_zh=title,status='success' if complete_groups and all(g['status']=='success' for g in groups) and len(figs)==expected else 'failed',command=f'python validation/plot_rerun_20261026/gsea/run.py --kind {kind} --jobs 3',inputs='restored gseapy.gene_set.prerank.report.csv, prerank_data.rnk, gene_sets.gmt; historical PDF names only select terms',input_nature_zh=nature,figure_count=len(figs),file_count=len(outputs),expected_historical_figure_count=expected,region_count=40,region_library_count=len(groups),max_absolute_ES_error=max((f['absolute_ES_error'] for f in figs),default=0),outputs=outputs,groups=groups))
 summary=dict(entries=entries,entrypoint_count=2,success_count=sum(e['status']=='success' for e in entries),figure_count=sum(e['figure_count'] for e in entries),file_count=sum(e['file_count'] for e in entries),restoration_manifest='validation/plot_rerun_20261026/gsea/restored_inputs.json',checks=['Every historical PDF name mapped to one saved report term.','Every running ES recomputed from full saved ranks/GMT with error <= 2e-6.','Saved NES, P and FDR retained; no new null distribution/permutations.','Every PNG decoded and every PDF header/EOF verified.'],notes=['Only CSV/RNK/GMT intermediate files restored. Historical images were not copied or read.','Custom matplotlib rendering replaces library styling while preserving full curves, hits, ranked metrics and saved statistics.','Enrichment plots were originally conditional on library plotting criteria; all 120 region/library outputs are audited, including those with zero historical curves.'])
 (OUT/'report.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n')
 lines=['# Thyroid local GSEA 作图重跑','',f'两个原始入口，共{summary["figure_count"]}幅图、{summary["file_count"]}个新PDF/PNG文件。历史曲线清单为pathway 204幅、marginal 368幅。','', '| 原始入口 | 内容 | 状态 | 新图/文件数 | 最大ES误差 |','|---|---|---|---:|---:|']
 for e in entries:lines.append(f'| `{Path(e["script"]).name}` | {e["description_zh"]} | {e["status"]} | {e["figure_count"]}/{e["file_count"]} | {e["max_absolute_ES_error"]:.3g} |')
 lines += ['', '所有40个注释区域、80个pathway区域/库及40个marginal区域均纳入检查。仅使用原历史PDF文件名选择应重画terms，没有读取或复制其图片内容。原始report/RNK/GMT的来源、目标和SHA256记录在[restored_inputs.json](restored_inputs.json)。全部输出位于`cases/aging/thyroid/local_O_info_regression/results/gsea_pathway`及`gsea_marginal`。','', '每条曲线从完整RNK和GMT计算weighted running ES，断言与保存ES差异≤2e-6；报告中保存的NES/P/FDR保持原值。使用matplotlib重新布局与渲染，不重跑回归、基因注释、信号提取或GSEA置换。所有PNG和PDF逐文件验证。','', '原`02_gsea_pathway.py`需要回归表及GENCODE GTF，并执行基因排名与置换；原`03_gsea_marginal.py`需要回归表及提取信号，执行bin分类、基因集构建及置换。因此这两个原入口不是纯绘图程序。本次专门提供只读中间结果的作图入口：','', '```bash','python validation/plot_rerun_20261026/gsea/run.py --jobs 3','```','', '逐图term、ES核验、输入组和结果SHA256在[report.json](report.json)。三个绘图进程各限制BLAS/OMP为1线程；每组执行证据在[groups_results](groups_results)。']
 (OUT/'REPORT.md').write_text('\n'.join(lines)+'\n')
 print(json.dumps({k:summary[k] for k in ['success_count','figure_count','file_count']}) ,flush=True)

if __name__=='__main__':
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--jobs',type=int,default=3);parser.add_argument('--kind',choices=['gsea_pathway','gsea_marginal']);parser.add_argument('--report-only',action='store_true');args=parser.parse_args()
 groups=json.loads((OUT/'groups.json').read_text());selected=[g for g in groups if args.kind is None or g['kind']==args.kind]
 if not args.report_only:
  with concurrent.futures.ProcessPoolExecutor(max_workers=args.jobs) as pool:
   for i,r in enumerate(pool.map(run_group,selected),1):print(i,len(selected),r['kind'],Path(r['input_folder']).parent.name,r['status'],len(r['figures']),r.get('error','')[:150],flush=True)
 records=[json.loads(p.read_text()) for p in sorted((OUT/'groups_results').glob('*.json'))]
 make_report(records)
