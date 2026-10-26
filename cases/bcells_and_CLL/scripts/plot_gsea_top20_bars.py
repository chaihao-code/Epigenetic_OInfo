#!/usr/bin/env python3
from pathlib import Path
import os
import csv
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Patch
P=Path(__file__).resolve().parents[1]
S=P/'results_gsea_rna_male/male_RNA_method_IGHV_adjusted'
O=S/'top20_pathways';O.mkdir(exist_ok=True)
# Selected representatives of the functional themes reported in the text.
selection=[
 ('oxidative phosphorylation','Oxidative\nphosphorylation','Energy metabolism'),
 ('mitochondrial ATP synthesis coupled electron transport','ATP synthesis /\nelectron transport','Energy metabolism'),
 ('cytoplasmic translation','Cytoplasmic\ntranslation','Protein synthesis / turnover'),
 ('mitochondrial translation','Mitochondrial\ntranslation','Protein synthesis / turnover'),
 ('proteasome complex','Proteasome\ncomplex','Protein synthesis / turnover'),
 ('chemokine receptor binding','Chemokine\nreceptor binding','Immune functions'),
 ('antigen processing and presentation of peptide antigen via MHC class I','MHC-I antigen\npresentation','Immune functions'),
 ('positive regulation of cytokine production','Positive regulation\nof cytokine production','Immune functions'),
 ('lymphocyte homeostasis','Lymphocyte\nhomeostasis','Immune functions'),
 ('regulation of G1/S transition of mitotic cell cycle','G1/S transition\nregulation','Cell cycle / survival'),
 ('regulation of signal transduction by p53 class mediator','p53 signaling\nregulation','Cell cycle / survival'),
 ('negative regulation of apoptotic process','Negative regulation\nof apoptosis','Cell cycle / survival'),
 ('extracellular matrix structural constituent','ECM structural\nconstituent','ECM / receptors'),
 ('G protein-coupled receptor activity','GPCR\nactivity','ECM / receptors')]
lookup={}
for ont in ['BP','MF','CC']:
 with (S/(ont+'_all_terms.tsv')).open() as h:
  for r in csv.DictReader(h,delimiter='\t'):
   r['ontology']=ont;lookup[r['GO_term']]=r
# Rank using ontology-specific FDR; ties by |NES| and GO ID.
top=sorted(lookup.values(),key=lambda r:(float(r['FDR q-val']),-abs(float(r['NES'])),r['GO_ID']))[:20]
import textwrap
selection=[(r['GO_term'],'\n'.join(textwrap.wrap(r['GO_term'],width=23)),r['ontology']) for r in top]
rows=[]
for term,label,category in selection:
 r=dict(lookup[term]);assert float(r['FDR q-val'])<.1
 r['display_label']=label.replace('\n',' ');r['category']=category
 # Visualization cap only: does not change stored statistics or claim a P-value bound.
 r['display_minus_log10_FDR']=min(4.,-np.log10(float(r['FDR q-val']))) if float(r['FDR q-val'])>0 else 4.
 r['empirical_FDR_zero']=float(r['FDR q-val'])==0
 rows.append(r)
pd.DataFrame(rows).to_csv(O/'top20_pathway_statistics.tsv',sep='\t',index=False)
font = Path(os.environ.get('BCELLS_ARIAL_FONT', str(Path.home() / '.local/share/fonts/arial.ttf')))
if font.exists():
    font_manager.fontManager.addfont(str(font))
plt.rcParams.update({'font.family':'Arial' if font.exists() else 'DejaVu Sans','font.size':5,'axes.labelsize':5,'axes.titlesize':5,'xtick.labelsize':5,'ytick.labelsize':5,'legend.fontsize':5,'pdf.fonttype':42})
colors={'BP':'#4aab77','MF':'#4d92bc','CC':'#e39c35'}
fig,ax=plt.subplots(figsize=(12/2.54,6/2.54))
x=np.arange(len(rows));heights=[r['display_minus_log10_FDR'] for r in rows]
ax.bar(x,heights,width=.68,color=[colors[r['category']] for r in rows])
for i,r in enumerate(rows):
 suffix='*' if r['empirical_FDR_zero'] else ''
 ax.text(i,heights[i]+.08,f"{float(r['NES']):+.2f}{suffix}",ha='center',va='bottom',fontsize=5)
ax.axhline(1,color='.5',ls='--',lw=.45)
ax.text(.995,1.04,'FDR = 0.10',transform=ax.get_yaxis_transform(),ha='right',fontsize=5)
ax.set_xticks(x);ax.set_xticklabels([label for term,label,cat in selection],rotation=65,ha='right',rotation_mode='anchor')
ax.set_ylabel('-log10(FDR)')
ax.set_title('Top 20 GO gene sets by ontology-specific FDR',pad=4)
ax.set_ylim(0,4.6);ax.set_xlim(-.6,len(rows)-.4)
ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
ax.tick_params(length=2,width=.4,pad=2)
fig.legend(handles=[Patch(color=c,label=k) for k,c in colors.items()],loc='lower center',bbox_to_anchor=(.5,.075),ncol=3,frameon=False,columnspacing=1.1,handlelength=1)
fig.text(.12,.025,'Labels: signed NES (+ R / - S). * Empirical FDR = 0; bar shown at plotting cap of 4.',fontsize=5)
fig.subplots_adjust(left=.10,right=.985,top=.88,bottom=.52)
for ext in ['pdf','png']:
 fig.savefig(O/('top20_GO_FDR_NES_12x6cm.'+ext),dpi=400)
(O/'caption_zh.txt').write_text('RNA方法及IGHV校正模型的GO FDR前20条目。BP/MF/CC合并后按各类别内FDR升序，同值按绝对NES降序、GO编号升序排序。柱高为-log10(FDR)，颜色表示GO类别，柱顶数字为带符号NES（正值R方向、负值S方向）。横线为FDR=0.10。经验FDR为0的条目以星号标记，柱高统一显示为4，仅作绘图截断，不表示真实FDR等于0.0001。GO类别内分别计算FDR；条目之间仍可能共享基因，且并非统一检验集合的全局FDR排名。\n')
print('Saved',O,flush=True)
