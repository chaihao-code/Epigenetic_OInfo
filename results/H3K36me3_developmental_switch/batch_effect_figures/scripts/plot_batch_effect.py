#!/usr/bin/env python3
"""Four-mark 2^4 source replacement, BH56 participation and strict set stability."""
import argparse
import csv
import gzip
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MARKS = ['ATAC','DNA_methylation','H3K9me3','H3K27me3','H3K36me3','H3K4me1','H3K4me3','H3K27ac']
STAGES = ['2cell','8cell','ICM']
LATE_STAGES = ['E11.5','E12.5','E13.5','E14.5','E15.5','E16.5','P0']


def canonical(name):
    match = re.search(r'H3K\d+(?:me\d|ac)|ATAC|CpG_methyl|DNA_methylation|WGBS',name)
    if not match: raise ValueError(name)
    return 'DNA_methylation' if match[0] in {'CpG_methyl','WGBS'} else match[0]


def bh(p):
    q=[1.]*len(p);running=1.
    order=sorted(range(len(p)),key=p.__getitem__)
    for rank in range(len(p),0,-1):
        i=order[rank-1];running=min(running,p[i]*len(p)/rank);q[i]=running
    return q


def table(name, rows, fields=None):
    path=ROOT/name;opener=gzip.open if path.suffix=='.gz' else open
    with opener(path,'wt',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),delimiter='\t')
        w.writeheader();w.writerows(rows)


def load(folder,name):
    meta=json.loads((folder/(name+'_metadata.json')).read_text())
    assert (meta['n_features'],meta['n_bins'],meta['n_permutations'],meta['permutation_level'],meta['permutation_scope'])==(8,5,10000,'all','block')
    p=folder/(name+'_null_summary.csv');groups=defaultdict(list)
    with p.open() as f:
        for r in csv.DictReader(f):
            if int(r['order'])!=3:continue
            signals=frozenset(canonical(x) for x in r['subset_names'].split('/'))
            assert len(signals)==3 and signals<=set(MARKS) and int(r['n_perm'])==10000
            obs,pval=float(r['observed_hoi']),float(r['p_down'])
            assert math.isfinite(obs) and 1/10001-1e-15<=pval<=1
            groups[r['block']].append(dict(block=r['block'],triplet='/'.join(sorted(signals)),observed_hoi=obs,p_down=pval,n_perm=10000))
    assert len(groups)==40
    records=[]
    for block,rows in sorted(groups.items()):
        assert len(rows)==56 and len({r['triplet'] for r in rows})==56
        for r,q in zip(rows,bh([r['p_down'] for r in rows])):
            r['q_down_block_bh56']=q;r['significant']=r['observed_hoi']<0 and q<=.05
        records.extend(rows)
    return records,dict(sample=name,path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())


def ratio(sets):
    intersection=set.intersection(*sets);union=set.union(*sets)
    return dict(n_intersection=len(intersection),n_union=len(union),intersection_over_union=len(intersection)/len(union) if union else '')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--validate-inputs',action='store_true')
    parser.add_argument('--plots-only',action='store_true',help='Redraw validated summaries without repeating FDR calculations')
    args=parser.parse_args()
    manifest=json.loads((ROOT/'manifest_16_fdr.json').read_text())
    assert len(manifest['tasks'])==48
    if args.plots_only:
        with (ROOT/'same_stage_mark_summary.tsv').open() as f:
            summaries=[{**r,'mask':int(r['mask']),'mean_fraction':float(r['mean_fraction'])} for r in csv.DictReader(f,delimiter='\t')]
        with (ROOT/'data/late_common_tissue_means.tsv').open() as f:
            late=[{**r,'mean_fraction':float(r['mean_fraction'])} for r in csv.DictReader(f,delimiter='\t')]
        def stability_table(name):
            with (ROOT/name).open() as f:
                return [{**r,'n_intersection':int(r['n_intersection']),'n_union':int(r['n_union']),
                         'intersection_over_union':float(r['intersection_over_union']) if r['intersection_over_union'] else ''}
                        for r in csv.DictReader(f,delimiter='\t')]
        plot(summaries,late,stability_table('stability_by_stage_mark.tsv'),stability_table('stability_by_stage.tsv'),manifest)
        print('Redrew the primary all16 participation bar figure.',flush=True)
        return
    if args.validate_inputs:
        for stage in STAGES:
            tasks=[t for t in manifest['tasks'] if t['stage']==stage]
            assert {t['mask'] for t in tasks}==set(range(16))
            for t in tasks:assert len(t['tracks'])==8 and all(Path(p).is_file() for p in t['tracks'])
        print('All 48 input panels and the 16-mask design validated.');return
    allrows=[];sources=[];sets_by_stage=defaultdict(dict);observed_by_stage=defaultdict(dict);participation=[]
    for t in manifest['tasks']:
        name=t['name'];rows,source=load(ROOT/'results'/name,name);sources.append(source)
        for r in rows:r.update(sample=name,stage=t['stage'],mask=t['mask'],variant=t['label'])
        allrows.extend(rows)
        sets_by_stage[t['stage']][t['mask']]={(r['block'],r['triplet']) for r in rows if r['significant']}
        observed_by_stage[t['stage']][t['mask']]={(r['block'],r['triplet']):r['observed_hoi'] for r in rows}
        for mark in MARKS:
            groups=defaultdict(list)
            for r in rows:
                if mark in r['triplet'].split('/'):groups[r['block']].append(r)
            for block,z in sorted(groups.items()):
                assert len(z)==21
                n=sum(r['significant'] for r in z)
                participation.append(dict(stage=t['stage'],sample=name,mask=t['mask'],variant=t['label'],mark=mark,block=block,n_significant=n,n_possible=21,fraction=n/21))
        print('Validated',name,flush=True)
    assert len(allrows)==48*40*56
    table('data/all_16variant_triplets_bh56.tsv.gz',allrows)
    table('data/mark_region_participation.tsv',participation)
    blocks=sorted({r['block'] for r in allrows});assert len(blocks)==40
    # Numerical validation for triplets whose signal files were not replaced.
    unaffected=[]
    for t in manifest['tasks']:
        original=observed_by_stage[t['stage']][0];current=observed_by_stage[t['stage']][t['mask']]
        differences=[abs(v-original[k]) for k,v in current.items() if not set(k[1].split('/')) & set(t['replaced_marks'])]
        diff=max(differences,default=0.)
        assert diff<1e-9,(t['name'],diff)
        unaffected.append(dict(sample=t['name'],n_unaffected=len(differences),max_abs_observed_difference=diff))
    table('data/unaffected_triplet_validation.tsv',unaffected)
    totals=[];byregion=[];bymark=[];support=[]
    for stage in STAGES:
        sets=[sets_by_stage[stage][mask] for mask in range(16)]
        totals.append(dict(stage=stage,n_variants=16,**ratio(sets)))
        for block in blocks:
            byregion.append(dict(stage=stage,block=block,n_variants=16,**ratio([{t for b,t in s if b==block} for s in sets])))
        for mark in MARKS:
            bymark.append(dict(stage=stage,mark=mark,n_variants=16,**ratio([{(b,t) for b,t in s if mark in t.split('/')} for s in sets])))
        for block,triplet in sorted(observed_by_stage[stage][0]):
            n=sum((block,triplet) in s for s in sets)
            support.append(dict(stage=stage,block=block,triplet=triplet,n_significant_variants=n,in_union=n>0,in_intersection=n==16))
    table('stability_by_stage.tsv',totals)
    table('stability_by_stage_region.tsv',byregion)
    table('stability_by_stage_mark.tsv',bymark)
    table('data/triplet_support_across_variants.tsv',support)
    # Current 10000-permutation late reference, balanced over the five tissues present at every stage.
    late=[];late_sources=[];late_root=Path(manifest['late_source'])
    for stage in LATE_STAGES:
        for tissue in manifest['late_reference_tissues']:
            name=stage+'_'+tissue;rows,source=load(late_root/name,name);late_sources.append(source)
            for r in rows:r.update(sample=name,stage=stage,tissue=tissue)
            late.extend(rows)
    assert len(late)==35*40*56
    table('data/late_reference_triplets_bh56.tsv.gz',late)
    late_cells=defaultdict(list)
    for r in late:
        for mark in r['triplet'].split('/'):
            late_cells[(mark,r['tissue'],r['stage'],r['block'])].append(r['significant'])
    late_tissue=[]
    for mark in MARKS:
        for tissue in manifest['late_reference_tissues']:
            values=[]
            for stage in LATE_STAGES:
                region=[]
                for block in blocks:
                    z=late_cells[(mark,tissue,stage,block)];assert len(z)==21
                    region.append(sum(z)/21)
                values.append(sum(region)/40)
            late_tissue.append(dict(mark=mark,tissue=tissue,mean_fraction=sum(values)/7,n_stages=7))
    table('data/late_common_tissue_means.tsv',late_tissue)
    summaries=[]
    for t in manifest['tasks']:
        for mark in MARKS:
            z=[r['fraction'] for r in participation if r['sample']==t['name'] and r['mark']==mark]
            assert len(z)==40
            summaries.append(dict(stage=t['stage'],mask=t['mask'],variant=t['label'],mark=mark,mean_fraction=sum(z)/40))
    table('same_stage_mark_summary.tsv',summaries)
    plot(summaries,late_tissue,bymark,totals,manifest)
    validation=dict(status='complete',n_runs=48,n_variants_per_stage=16,n_permutations=10000,BH_family_size=56,
                    n_variant_tests=len(allrows),n_late_tests=len(late),stability=totals,
                    n_significant=sum(r['significant'] for r in allrows),
                    max_unaffected_observed_difference=max(r['max_abs_observed_difference'] for r in unaffected),
                    sources=sources,late_sources=late_sources)
    (ROOT/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    print('STABILITY',json.dumps(totals),flush=True)
    print('Completed the primary BH56 participation figure and numeric stability summaries.',flush=True)


def plot(summaries,late_tissue,bymark,totals,manifest):
    import os
    os.environ.setdefault('MPLCONFIGDIR','/tmp/bulk16-support-mpl')
    os.environ.setdefault('XDG_CACHE_HOME','/tmp/bulk16-support-cache')
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    font_path=font_manager.findfont('Arial',fallback_to_default=False)
    plt.rcParams.update({'font.family':'Arial','font.size':8,
                         'axes.titlesize':8,'axes.labelsize':8,'xtick.labelsize':8,'ytick.labelsize':8,
                         'legend.fontsize':8,'pdf.fonttype':42,'ps.fonttype':42,'axes.linewidth':.5})
    variants=sorted(manifest['variants'],key=lambda v:(v['n_bulk'],v['mask']))
    short={'H3K36me3':'K36','H3K9me3':'K9','H3K27me3':'K27','H3K4me3':'K4'}
    labels=['Original' if not v['replaced_marks'] else '+'.join(short[m] for m in v['replaced_marks']) for v in variants]
    # Light, moderately saturated colors without the fluorescent turbo midpoint.
    colors=['#999999','#758EC4','#819DD0','#69A9D3','#63B5C5','#53B5B3',
            '#79BFAC','#94C19F','#B0C38C','#CDD08D','#E8CA79','#F2BD61',
            '#EBA16C','#E59876','#E58078','#CD7F91']
    late_mean=[np.mean([r['mean_fraction'] for r in late_tissue if r['mark']==m]) for m in MARKS]
    late_sd=[np.std([r['mean_fraction'] for r in late_tissue if r['mark']==m],ddof=1) for m in MARKS]
    peak=max(max(r['mean_fraction'] for r in summaries),max(np.array(late_mean)+np.array(late_sd)))
    upper=min(1.,max(.1,np.ceil(peak*1.1*10)/10))
    # Each plotting rectangle is exactly 3 cm high; additional space is for labels/legend.
    height_cm=15
    fig=plt.figure(figsize=(18/2.54,height_cm/2.54))
    axes=[]
    for bottom_cm in [11.4,7.9,4.4]:
        axes.append(fig.add_axes([.10,bottom_cm/height_cm,.89,3/height_cm],
                                sharex=axes[0] if axes else None,
                                sharey=axes[0] if axes else None))
    width=.8/17
    for ax,stage in zip(axes,STAGES):
        for j,v in enumerate(variants):
            values=[next(r['mean_fraction'] for r in summaries if r['stage']==stage and r['mask']==v['mask'] and r['mark']==m) for m in MARKS]
            ax.bar(np.arange(8)+(j-8)*width,values,width=width*.93,color=colors[j],label=labels[j])
        xpos=np.arange(8)+8*width
        ax.bar(xpos,late_mean,width=width*.93,color='black',label='Late (5 tissues)')
        ax.errorbar(xpos,late_mean,yerr=late_sd,fmt='none',ecolor='black',elinewidth=.5,capsize=1)
        ax.set_title({'2cell':'2-cell','8cell':'8-cell','ICM':'ICM'}[stage],loc='left',pad=3)
        ax.set_ylim(0,upper);ax.spines[['top','right']].set_visible(False)
        ax.tick_params(width=.5,length=2)
    for ax in axes[:-1]:ax.tick_params(labelbottom=False)
    axes[-1].set_xticks(np.arange(8),['ATAC','DNA\nmethylation']+MARKS[2:],rotation=45,ha='right')
    fig.text(.025,.62,'Significant / possible triplets',rotation=90,va='center',fontsize=8)
    fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.54,.035),
               ncol=4,frameon=False,handlelength=.8,handletextpad=.3,columnspacing=.55,labelspacing=.35)
    fig.text(.10,.009,'BH56 FDR ≤ 0.05; Ω < 0. Late: mean ± SD of 5 tissues.',fontsize=8)
    fig.canvas.draw()
    assert all(abs(ax.get_position().height*height_cm-3)<1e-10 for ax in axes)
    for ext in ['pdf','png']:fig.savefig(ROOT/'plots'/('same_stage_bulk_replacement_all_marks.'+ext),dpi=300)
    plt.close(fig)
    (ROOT/'data/bulk_replacement_plot_format.json').write_text(json.dumps(dict(
        width_cm=18,total_height_cm=height_cm,subplot_height_cm=3,n_subplots=3,legend_columns=4,
        font_family='Arial',font_file=font_path,font_size_pt=8,legend='bottom',
        variant_colors=dict(zip(labels,colors)),late_color='#000000',
        only_primary_figure=True),indent=2)+'\n')


if __name__=='__main__':main()
