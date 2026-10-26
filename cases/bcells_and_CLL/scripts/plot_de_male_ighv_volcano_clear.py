#!/usr/bin/env python3
from pathlib import Path
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
P=Path(__file__).resolve().parents[1]
O=P/'results_rna_de_male_ighv'
font = Path(os.environ.get('BCELLS_ARIAL_FONT', str(Path.home() / '.local/share/fonts/arial.ttf')))
if font.exists():
    font_manager.fontManager.addfont(str(font))
plt.rcParams.update({'font.family':'Arial' if font.exists() else 'DejaVu Sans','font.size':5,'axes.labelsize':5,'axes.titlesize':5,'xtick.labelsize':5,'ytick.labelsize':5,'legend.fontsize':5,'pdf.fonttype':42})
d=pd.read_csv(O/'rna_method_ighv_adjusted_all_genes.tsv',sep='\t')
d=d[np.isfinite(d.pvalue)&np.isfinite(d.log2FoldChange)].copy()
d['y']=-np.log10(d.pvalue.clip(lower=1e-300))
size_cm=float(os.environ.get('VOLCANO_SIZE_CM','6'))
fig,ax=plt.subplots(figsize=(size_cm/2.54,size_cm/2.54))
sig=d.padj.lt(.1)
for mask,color,label in [(~sig,'#b9b9b9','Other tested genes'),(sig&d.log2FoldChange.gt(0),'#bc4749','R higher (25)'),(sig&d.log2FoldChange.lt(0),'#457b9d','S higher (35)')]:
 ax.scatter(d.loc[mask,'log2FoldChange'],d.loc[mask,'y'],s=2,alpha=.6,color=color,label=label,rasterized=True)
cutoff=float(d.loc[sig,'pvalue'].max())
assert int(sig.sum())==60
ax.axhline(-np.log10(cutoff),color='.4',ls='--',lw=.4)
ax.text(.98,-np.log10(cutoff),'FDR < 0.1',transform=ax.get_yaxis_transform(),ha='right',va='bottom',fontsize=5)
print('Nominal P threshold:',cutoff,flush=True)
ax.axvline(0,color='.7',lw=.4)
label_genes=['CCL3','CCL4','LPL','CLEC4C','CLEC2B','PTPRS','GDPD5','RGS12','ADHFE1','D2HGDH']
selected=d[d.gene_name.isin(label_genes)].copy()
assert len(selected)==10 and selected.padj.lt(.1).all()
selected.to_csv(O/'volcano_labeled10_statistics.tsv',sep='\t',index=False)
xlo,xhi=ax.get_xlim()
for positive in [False,True]:
 group=selected[selected.log2FoldChange.gt(0).eq(positive)].sort_values('y',ascending=False)
 text_x=.96 if positive else .04
 label_y=np.linspace(.74,.51,len(group))
 for row,ypos in zip(group.itertuples(),label_y):
  ax.scatter(row.log2FoldChange,row.y,s=6,facecolors='none',edgecolors='#333333',linewidths=.4,zorder=4)
  ax.annotate(row.gene_name,(row.log2FoldChange,row.y),xytext=(text_x,ypos),textcoords='axes fraction',
   fontsize=5,ha='right' if positive else 'left',va='center',
   arrowprops=dict(arrowstyle='-',color='.55',lw=.25),
   zorder=5)
ax.set_xlabel('Adjusted log2 fold change (R / S)')
ax.set_ylabel('-log10(nominal P value)')
ax.set_title('Differential gene expression\nbetween CLL-R and CLL-S cells',fontsize=5,pad=4)
ax.legend(frameon=False,loc='upper left',markerscale=2)
ax.tick_params(width=.4,length=2)
fig.subplots_adjust(left=.17,right=.985,bottom=.15,top=.89)
fig.canvas.draw()
renderer=fig.canvas.get_renderer()
# All gene-label boxes must remain inside the axes; the title inside the figure.
for text in ax.texts:
 if text.get_text() in label_genes:
  box=text.get_window_extent(renderer)
  assert ax.bbox.contains(box.x0,box.y0) and ax.bbox.contains(box.x1,box.y1), text.get_text()
box=ax.title.get_window_extent(renderer)
assert fig.bbox.contains(box.x0,box.y0) and fig.bbox.contains(box.x1,box.y1)
for ext in ['pdf','png']:
 fig.savefig(O/(f'RNA_IGHV_adjusted_volcano_nominalP_labeled10_{size_cm:g}cm_clear.'+ext),dpi=300)
(O/f'RNA_IGHV_adjusted_volcano_nominalP_labeled10_{size_cm:g}cm_clear_caption_zh.txt').write_text(f'横轴为 RNA方法及IGHV校正模型的 R/S log2倍数变化；纵轴为未进行多重检验校正的 nominal P 值负对数。红/蓝色分别为实际 padj<0.10 的 R/S 高表达基因。横线对应这些显著基因中最大的 nominal P 值 {cutoff:.8g}，仅适用于本批数据。灰色包括不显著及因独立过滤而无有效 padj 的基因；后者不依横线判定显著。\n')
