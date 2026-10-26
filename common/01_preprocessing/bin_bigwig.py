#!/usr/bin/env python3
"""Bin selected bigWig regions into four-column .bp tracks, with covered-base means."""
from pathlib import Path
import argparse
import math
import pyBigWig

def bin_bigwig(source, regions, output, bin_size=1000):
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    with pyBigWig.open(str(source)) as bw, Path(regions).open() as inp, output.open('w') as out:
        for line in inp:
            if not line.strip() or line.startswith('#'):continue
            chrom,start,end=line.split()[:3];start,end=int(start),int(end)
            if start < 0 or end <= start or end > bw.chroms(chrom):
                raise ValueError(f'Invalid region {chrom}:{start}-{end}')
            for left in range(start,end,bin_size):
                right=min(left+bin_size,end)
                mean=bw.stats(chrom,left,right,type='mean',exact=True)[0]
                value=0.0 if mean is None or not math.isfinite(mean) else float(mean)
                out.write(f'{chrom}\t{left}\t{right}\t{value:.17g}\n')
    return output

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',required=True,type=Path)
    p.add_argument('--regions',required=True,type=Path,help='BED3 regions; zero-based, half-open')
    p.add_argument('--output',required=True,type=Path)
    p.add_argument('--bin-size',type=int,default=1000)
    args=p.parse_args()
    if args.bin_size <= 0:p.error('--bin-size must be positive')
    bin_bigwig(args.input,args.regions,args.output,args.bin_size)
if __name__=='__main__':main()
