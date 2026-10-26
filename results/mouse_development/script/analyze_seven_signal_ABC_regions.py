import csv,gzip,collections,statistics,math,json
from pathlib import Path
base=Path(__file__).resolve().parents[1] / 'synergistic_triplets'
stages=['zygote','2cell','4cell','8cell','morula','ICM','TE','E6.5','E7.0','E7.5','E11.5','E12.5','E13.5','E14.5','E15.5','E16.5','P0']
maxima={'A':9,'B':16,'C':10}; samples=collections.defaultdict(set); groups=collections.defaultdict(list)
with gzip.open(base/'all_triplets_fdr_10000.tsv.gz','rt') as f:
 for r in csv.DictReader(f,delimiter='\t'):
  marks=r['triplet'].split('/'); assert 'H3K36me3' not in marks
  repressive=bool(set(marks)&{'H3K9me3','H3K27me3'})
  cls='A' if 'CpG_methyl' in marks and repressive else 'B' if repressive else 'C'
  assert cls==r['triplet_class']
  groups[r['sample'],r['block'],cls].append(r);samples[r['stage']].add(r['sample'])
blocks=sorted({k[1] for k in groups});assert len(blocks)==40
counts={}
for (sample,b,c),rr in groups.items():
 assert len(rr)==maxima[c]
 counts[sample,b,c]=sum(float(r['observed_hoi'])<0 and float(r['q_down_sample_block_bh'])<=.05 for r in rr)/maxima[c]
lookup={}
for s in stages:
 for b in blocks:
  for c in 'ABC':lookup[s,b,c]=statistics.mean(counts[x,b,c] for x in samples[s])
# Match the retained final selected-region figure.
with open(base/'triplet_ABC_selected_region_trajectories_7signals_summary.tsv') as f:
 for r in csv.DictReader(f,delimiter='\t'): assert abs(lookup[r['stage'],r['block'],r['triplet_class']]-float(r['mean_relative']))<1e-12
out=[]
for label,ss in [('all',stages),('preimplantation',stages[:7]),('gastrulation',stages[7:10]),('late',stages[10:])]:
 for b in blocks:
  means={c:statistics.mean(lookup[s,b,c] for s in ss) for c in 'ABC'}
  d={c:statistics.mean(abs(lookup[s,b,c]-lookup[s,b,cc]) for s in ss for cc in 'ABC' if cc!=c) for c in 'ABC'}
  pairgaps=[abs(means[a]-means[c]) for a,c in [('A','B'),('A','C'),('B','C')]]
  out.append(dict(period=label,block=b,**means,range=max(means.values())-min(means.values()),min_pair_gap=min(pairgaps),mean_timewise_pair_gap=statistics.mean(abs(lookup[s,b,a]-lookup[s,b,c]) for s in ss for a,c in [('A','B'),('A','C'),('B','C')]),order='>'.join(sorted('ABC',key=means.get,reverse=True)),A_top_stages=sum(lookup[s,b,'A']>max(lookup[s,b,'B'],lookup[s,b,'C']) for s in ss),B_top_stages=sum(lookup[s,b,'B']>max(lookup[s,b,'A'],lookup[s,b,'C']) for s in ss),C_top_stages=sum(lookup[s,b,'C']>max(lookup[s,b,'A'],lookup[s,b,'B']) for s in ss)))
for period in ['all','late']:
 print('\nPERIOD',period,'TOP minimum pair gap')
 for r in sorted([r for r in out if r['period']==period],key=lambda r:r['min_pair_gap'],reverse=True)[:10]:print(r['block'],*[round(r[c],3) for c in 'ABC'],'gap',round(r['min_pair_gap'],3),r['order'])
 print('B or C ABOVE BOTH, period mean')
 for r in out:
  if r['period']==period and (r['B']>max(r['A'],r['C']) or r['C']>max(r['A'],r['B'])):print(r['block'],*[round(r[c],3) for c in 'ABC'],r['order'])
selected=['exon_purified_bins_AND_SINE','distal_not_ctcf_bound_enhancer_purified_bins_AND_SINE','proximal_not_ctcf_bound_enhancer_purified_bins_AND_CpG_island','distal_not_ctcf_bound_enhancer_purified_bins_AND_CpG_island','exon_purified_bins_AND_LINE','intron_purified_bins_AND_SINE']
print('\nCANDIDATES BY PERIOD')
for b in selected:
 print(b)
 for r in out:
  if r['block']==b:print(r['period'],*[round(r[c],3) for c in 'ABC'],'top stages',*[r[c+'_top_stages'] for c in 'ABC'])
print('\nDISTAL CpG timecourse')
for s in stages: print(s,*[round(lookup[s,selected[3],c],3) for c in 'ABC'])
path=base/'seven_signal_ABC_region_separation.tsv'
with path.open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
print('saved',path)
metadata=dict(source=str(base/'all_triplets_fdr_10000.tsv.gz'), n_signals=7,
              repressive_marks=['H3K9me3','H3K27me3'], H3K36me3_present=False,
              theoretical_max=maxima, n_samples=sum(len(z) for z in samples.values()),
              samples_per_stage={s:len(samples[s]) for s in stages},
              criterion='observed_hoi < 0 and BH35 q <= 0.05',
              aggregation='Equal mean across available samples within each stage, then equal mean across stages within each period',
              minimum_pair_gap='Minimum of the three absolute pairwise differences between period-mean A/B/C proportions',
              scope='Descriptive exploration across all 40 blocks; region selection is not a significance test',
              final_selected_region_trajectory_means_matched=True)
(base/'seven_signal_ABC_region_separation_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
