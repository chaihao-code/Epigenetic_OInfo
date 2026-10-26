#!/usr/bin/env python3
"""Complete thyroid donor: unchanged whole-genome grid, original 40 classes, shared global/local Omega."""
from pathlib import Path
import argparse,gzip,hashlib,json,os,shutil,subprocess,sys,time
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'demo/scripts'))
from prepare_from_bigwig import prepare

def compare_input(new, old):
    count=missing=changed=0;maximum=0.
    names=['chrom','start','end','value']
    aa=pd.read_csv(new,sep='\t',header=None,names=names,chunksize=250000)
    bb=pd.read_csv(old,sep='\t',header=None,names=names,chunksize=250000)
    for a,b in zip(aa,bb,strict=True):
        if not a.iloc[:,:3].equals(b.iloc[:,:3]):raise ValueError('Original fixed-grid coordinates changed')
        x=a.value.to_numpy(float);y=b.value.to_numpy(float)
        if not np.array_equal(np.isnan(x),np.isnan(y)):raise ValueError('NA coverage differs from the original input')
        good=np.isfinite(x)&np.isfinite(y);delta=np.abs(x[good]-y[good])
        maximum=max(maximum,float(delta.max(initial=0)));changed+=int(np.count_nonzero(delta));missing+=int((~good).sum());count+=len(a)
    return dict(n_rows=count,n_NA=missing,n_different_signal_values=changed,max_abs_signal_difference=maximum)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--from-bigwig',action='store_true')
    p.add_argument('--output',type=Path)
    p.add_argument('--n-permutations',type=int,default=10000)
    p.add_argument('--jobs',type=int,default=32)
    p.add_argument('--reuse-raw-binned',action='store_true',help='Reuse a completed fixed-grid raw binning stage after checking its metadata')
    p.add_argument('--triplets',action='store_true',help='Optional generic triplet/FDR demonstration; default is the thyroid full-system null')
    args=p.parse_args();cfg=json.loads((ROOT/'demo/config/demo.json').read_text())
    mode='from_bigwig' if args.from_bigwig else 'from_bins'
    out=(args.output or ROOT/'demo/results'/mode).resolve();out.mkdir(parents=True,exist_ok=True)
    tracks=[x for x in cfg['signals'] if x['mark']!='Input']
    # Restore only compressed copies of the exact original prepared tracks when local links are unavailable.
    for entry in tracks:
        source=ROOT/entry['binned_file']
        if not source.exists():
            if source.is_symlink():source.unlink()
            with gzip.open(str(source)+'.gz','rb') as inp,source.open('wb') as dest:shutil.copyfileobj(inp,dest)
    started=time.time()
    manifest=prepare(out/'preprocessing',min(args.jobs,8),args.reuse_raw_binned) if args.from_bigwig else ROOT/'demo/input/signals.txt'
    actual=[Path(x) if Path(x).is_absolute() else ROOT/x for x in manifest.read_text().splitlines() if x.strip()]
    comparisons=[]
    for path,entry in zip(actual,tracks,strict=True):
        if args.from_bigwig:comparisons.append(dict(mark=entry['mark'],**compare_input(path,ROOT/entry['binned_file'])))
        else:
            digest=hashlib.sha256(path.read_bytes()).hexdigest()
            if digest!=entry['binned_sha256']:raise ValueError('Original prepared-input checksum mismatch')
    maskgroups=';'.join(str(ROOT/x) for x in cfg['functional_masks'])+'::BLOCKGROUP::'+';'.join(str(ROOT/x) for x in cfg['repeat_masks'])
    command=[sys.executable,str(ROOT/'common/02_oinformation/compute_hoi_oinfo.py'),
       '--input-file',str(manifest),'--output-dir',str(out/'oinfo'),'--sample-name',cfg['sample'],
       '--block-groups',maskgroups,'--n-bins','5','--n-max','100000','--random-seed','42','--drop-signals','',
       '--n-permutations',str(args.n_permutations),'--permutation-seed','1234','--permutation-scope','block',
       '--permutation-level','all' if args.triplets else 'system',
       '--permutation-min-order','3','--permutation-max-order','3','--fdr-family','sample-block',
       '--save-local-oinfo','--save-sampled-bins','--n-jobs',str(args.jobs)]
    with (out/'run.log').open('w') as log:subprocess.run(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,check=True)
    base=out/'oinfo'/cfg['sample'];oi=pd.read_csv(str(base)+'_oinfo.csv')
    identity=pd.read_csv(str(base)+'_local_global_check.tsv',sep='\t')
    full=oi[oi.n_missing.eq(0)]
    assert len(oi)==4800 and len(full)==40 and len(identity)==40
    assert full.order.eq(7).all() and identity.absolute_difference.max()<1e-8
    assert np.max(np.abs(oi.hoi-(oi.tc-oi.dtc)))<1e-10
    checks=dict(n_classes=40,n_features=7,n_observed_rows=len(oi),n_sampled_loci=int(identity.n_loci.sum()),
                maximum_local_global_difference=float(identity.absolute_difference.max()))
    if args.triplets:
        sig=pd.read_csv(str(base)+'_triplet_significance.tsv',sep='\t')
        assert len(sig)==1400 and sig.n_perm.eq(args.n_permutations).all() and sig.fdr_family_size.eq(35).all()
        assert sig.significant_synergy.eq(sig.observed_hoi.lt(0)&sig.q_down_BH.le(.05)).all()
        checks['n_triplet_tests']=len(sig)
    else:
        null=pd.read_csv(str(base)+'_permutation_summary.csv')
        assert len(null)==40 and null.n_perm.eq(args.n_permutations).all()
        assert null.p_down.between(1/(args.n_permutations+1),1).all()
        checks['n_full_system_permutation_tests']=40
    reference=ROOT/'cases/aging/thyroid/deconvolution/data/processed/oinfo_long.tsv'
    ref=pd.read_csv(reference,sep='\t');ref=ref[ref.sample_id.eq(cfg['donor'])]
    matched=full.merge(ref,left_on='block',right_on='region',suffixes=('_new','_reference'))
    checks['original_result_classes_matched']=len(matched)
    if len(matched):checks['maximum_difference_from_original_oinfo']=float(np.abs(matched.hoi_new-matched.hoi_reference).max())
    report=dict(status='passed',sample=cfg['sample'],mode=mode,slurm_job_id=os.environ.get('SLURM_JOB_ID'),
        fixed_grid_sha256=cfg['fixed_grid_sha256'],n_genome_bins=cfg['n_genome_bins'],n_autosomal_bins=cfg['n_autosomal_bins'],
        n_permutations=args.n_permutations,permutation_level='all-triplets' if args.triplets else 'full-system',
        input_comparisons=comparisons,elapsed_seconds=round(time.time()-started,3),
        python=sys.version,python_executable=sys.executable,numpy=np.__version__,pandas=pd.__version__,
        engine_sha256=hashlib.sha256((ROOT/'common/02_oinformation/compute_hoi_oinfo.py').read_bytes()).hexdigest(),
        command=command,checks=checks)
    (out/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    full.to_csv(out/'global_oinfo_by_class.tsv',sep='\t',index=False)
    print(json.dumps(dict(status=report['status'],sample=report['sample'],mode=mode,**checks),indent=2))
if __name__=='__main__':main()
