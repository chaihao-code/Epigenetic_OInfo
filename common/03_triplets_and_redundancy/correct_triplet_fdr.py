#!/usr/bin/env python3
"""BH-correct existing triplet null summaries; never recompute observed/permuted Omega."""
from pathlib import Path
import argparse
import json
import sys
import pandas as pd
ROOT = next(p for p in Path(__file__).resolve().parents if (p / '.epigenetic_oinfo_root').is_file())
sys.path.insert(0, str(ROOT / 'common/02_oinformation'))
from compute_hoi_oinfo import correct_triplet_fdr

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',nargs='+',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path)
    p.add_argument('--family',choices=['sample-block','sample','global'],default='sample-block')
    p.add_argument('--alpha',type=float,default=.05)
    args=p.parse_args();frames=[]
    for source in args.input:
        df=pd.read_csv(source,sep='\t' if '.tsv' in source.name else ',')
        if 'sample' not in df:
            sample=source.name.removesuffix('_null_summary.csv').removesuffix('_triplet_significance.tsv')
            df.insert(0,'sample',sample)
        frames.append(df)
    frame=pd.concat(frames,ignore_index=True)
    if frame.duplicated(['sample','block','subset_names']).any():
        raise ValueError('Duplicate sample/class/triplet rows; provide each complete testing family once')
    result=correct_triplet_fdr(frame,family=args.family,alpha=args.alpha)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    result.to_csv(args.output,sep='\t',index=False)
    args.output.with_suffix('.metadata.json').write_text(json.dumps(dict(
        family=args.family,alpha=args.alpha,n_tests=len(result),
        definition='BH on every lower-tail triplet P in the chosen family, then observed_hoi < 0 and q <= alpha',
        input_files=[str(x) for x in args.input]),indent=2)+'\n')
    print(f'Wrote {args.output}: {len(result)} triplets')
if __name__=='__main__':main()
