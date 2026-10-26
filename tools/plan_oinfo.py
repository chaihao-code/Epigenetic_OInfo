#!/usr/bin/env python3
"""Preview the single supported O-information command; no analysis is executed."""
import argparse,json,shlex,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def main():
    settings=json.loads((ROOT/'config/cases.json').read_text())
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--case',required=True,choices=settings)
    p.add_argument('--sample',required=True)
    p.add_argument('--jobs',type=int,default=4)
    p.add_argument('--n-permutations',type=int,default=10000)
    p.add_argument('--n-bins',type=int,default=5)
    args=p.parse_args();cfg=settings[args.case]
    if args.case == 'H3K36me3_developmental_switch':
        p.error('Follow the case README: this cohort reuses its specific pipeline/manifest.')
    manifest=Path(cfg['input_dir'])/(args.sample+'.txt')
    if not (ROOT/manifest).is_file():p.error(f'No manifest: {manifest}')
    signals=[x.strip() for x in (ROOT/manifest).read_text().splitlines() if x.strip()]
    missing=[x for x in signals if not (ROOT/x).is_file()]
    f=['promoter_bins.tsv','proximal_ctcf_bound_enhancer_purified_bins.tsv','proximal_not_ctcf_bound_enhancer_purified_bins.tsv','distal_ctcf_bound_enhancer_purified_bins.tsv','distal_not_ctcf_bound_enhancer_purified_bins.tsv','exon_purified_bins.tsv','intron_purified_bins.tsv']
    if cfg['genome']=='hg38' and args.case!='aging/thyroid':f[0]='promoter_purified_bins.tsv'
    r=['CpG_island_bins.tsv','LINE_bins.tsv','SINE_bins.tsv','LTR_bins.tsv']
    if args.case=='bcells_and_CLL':f.append('None_func_bins.tsv');r.append('None_rep_bins.tsv')
    prefix=Path('data/genome_features')/cfg['genome']
    for name in f+r:
        if not (ROOT/prefix/name).is_file():missing.append(str(prefix/name))
    groups=';'.join(str(prefix/n) for n in f)+'::BLOCKGROUP::'+';'.join(str(prefix/n) for n in r)
    command=[sys.executable,'common/02_oinformation/compute_hoi_oinfo.py','--input-file',str(manifest),
        '--output-dir',f'cases/{args.case}/outputs/oinfo/{args.sample}','--sample-name',args.sample,
        '--n-bins',str(args.n_bins),'--n-max','100000','--random-seed','42','--block-groups',groups,
        '--drop-signals','','--n-permutations',str(args.n_permutations),'--permutation-level',cfg.get('permutation_level','all'),
        '--permutation-min-order','3','--permutation-max-order','3','--permutation-seed','1234',
        '--fdr-family',cfg['fdr_family'],'--n-jobs',str(args.jobs)]
    print('Working directory:',ROOT);print(shlex.join(command))
    print(f'Input paths: {len(signals)} signals; {len(missing)} missing signal/mask paths.')
    for x in missing:print('MISSING',x)
    return bool(missing)
if __name__=='__main__':raise SystemExit(main())
