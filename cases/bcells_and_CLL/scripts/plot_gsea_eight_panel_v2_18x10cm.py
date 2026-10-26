#!/usr/bin/env python3
from pathlib import Path
import os
import csv
import textwrap
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
P=Path(__file__).resolve().parents[1]
S=P/'results_gsea_rna_male/male_RNA_method_IGHV_adjusted'
O=S/'eight_panel_GSEA_v2_18x10cm';O.mkdir(exist_ok=True)
PANELS=[('BP','GO:0006119','Oxidative phosphorylation'),('BP','GO:0042775','ATP synthesis/electron transport'),('BP','GO:0002181','Cytoplasmic translation'),('MF','GO:0042379','Chemokine receptor binding'),('BP','GO:0002474','MHC-I antigen presentation'),('BP','GO:2000045','G1/S transition regulation'),('BP','GO:1901796','p53 signaling regulation'),('MF','GO:0005201','ECM structural constituent')]
font = Path(os.environ.get('BCELLS_ARIAL_FONT', str(Path.home() / '.local/share/fonts/arial.ttf')))
if font.exists():
    font_manager.fontManager.addfont(str(font))
plt.rcParams.update({'font.family':'Arial' if font.exists() else 'DejaVu Sans','font.size':5,'axes.labelsize':5,'xtick.labelsize':5,'ytick.labelsize':5,'axes.titlesize':6,'pdf.fonttype':42,'axes.linewidth':.35})
lookup={};sets={};ranks={}
for ont in ['BP','MF']:
 with (S/(ont+'_all_terms.tsv')).open() as h:
  for r in csv.DictReader(h,delimiter='\t'):lookup[r['GO_ID']]=r
 with (S/(ont+'.gmt')).open() as h:
  for line in h:
   a=line.rstrip().split('\t');sets[a[0].split('|')[0]]=set(a[2:])
 with (S/ont/'prerank_data.rnk').open() as h:
  ranks[ont]=[(l.split('\t')[0],float(l.split('\t')[1])) for l in h]
with (O/'panel_statistics.tsv').open('w') as h:
 fields=['panel','display_name','GO_ID','GO_term','NES','FDR q-val','NOM p-val']
 w=csv.DictWriter(h,fieldnames=fields,delimiter='\t');w.writeheader()
 for i,(ont,gid,label) in enumerate(PANELS,1):
  r=lookup[gid];w.writerow(dict(panel=i,display_name=label,**{k:r[k] for k in fields[2:]}))
for width,height in [(18,10)]:
 fig=plt.figure(figsize=(width/2.54,height/2.54))
 outer=fig.add_gridspec(2,4,left=.065,right=.985,bottom=.13,top=.82,wspace=.20,hspace=.55)
 for i,(ont,gid,label) in enumerate(PANELS):
  r=lookup[gid];rank=ranks[ont];scores=np.array([v for g,v in rank]);hit=np.array([g in sets[gid] for g,v in rank])
  weights=np.abs(scores)*hit
  run=np.cumsum(weights/weights.sum()-(~hit)/(~hit).sum())
  assert np.isclose(run[np.argmax(np.abs(run))],float(r['ES']),atol=1e-7)
  inner=outer[i//4,i%4].subgridspec(3,1,height_ratios=[3,.55,1.3],hspace=.10)
  a=fig.add_subplot(inner[0]);b=fig.add_subplot(inner[1],sharex=a);c=fig.add_subplot(inner[2],sharex=a)
  x=np.arange(1,len(rank)+1)
  a.plot(x,run,color='#38964a',lw=.6);a.axhline(0,color='.6',lw=.3)
  a.set_title(gid+'\n'+'\n'.join(textwrap.wrap(label,25)),fontsize=6,pad=3)
  q=float(r['FDR q-val']);fdr='0*' if q==0 else f'{q:.2g}'
  a.text(.98,.97,f"NES={float(r['NES']):.2f}\nFDR={fdr}",transform=a.transAxes,ha='right',va='top',fontsize=5)
  b.vlines(x[hit],.22,1,lw=.25,color='black');b.imshow(scores.reshape(1,-1),aspect='auto',cmap='RdBu_r',extent=[1,len(rank),0,.18],vmin=-max(abs(scores)),vmax=max(abs(scores)))
  b.set_ylim(0,1);b.set_yticks([])
  c.fill_between(x,scores,color='#c6d1db',lw=0);c.axhline(0,color='.6',lw=.3)
  c.set_xlim(1,len(rank));c.set_xticks([1,9000,17787]);c.set_xticklabels(['0','9k','18k'])
  for ax in [a,b,c]:ax.tick_params(length=1.5,width=.35,pad=1)
  a.tick_params(labelbottom=False);b.tick_params(labelbottom=False)
  if i%4==0:a.set_ylabel('Running ES');c.set_ylabel('Wald stat.')
  c.set_xlabel('Gene rank (R → S)',labelpad=1)
 fig.text(.065,.025,'RNA method + IGHV adjusted. Positive NES: R; negative NES: S. * Empirical FDR = 0 (1,000 permutations).',fontsize=5)
 for ext in ['pdf','png']:fig.savefig(O/f'GSEA_8pathways_{width}x{height}cm.{ext}',dpi=400)
 plt.close(fig)
print('Saved eight-panel figures:',O,flush=True)
