#!/usr/bin/env python3
"""Use the original immutable hg38 windows and the original binning/normalization arithmetic."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import argparse,gzip,hashlib,json,os,shutil,subprocess
ROOT=Path(__file__).resolve().parents[2]
def prepare(output, workers=8, reuse_raw=False):
    cfg=json.loads((ROOT/'demo/config/demo.json').read_text())
    output=Path(output);raw=output/'raw_binned';normalized=output/'input_normalized'
    raw.mkdir(parents=True,exist_ok=True);normalized.mkdir(exist_ok=True)
    grid=output/'fixed_windows.bed'
    digest=hashlib.sha256()
    with gzip.open(ROOT/cfg['fixed_grid'],'rb') as inp,grid.open('wb') as out:
        for block in iter(lambda:inp.read(1024*1024),b''):out.write(block);digest.update(block)
    if digest.hexdigest()!=cfg['fixed_grid_sha256']:raise ValueError('The original fixed genomic grid changed')
    bed4=output/'fixed_windows.named.bed'
    with grid.open() as inp,bed4.open('w') as out:
        for line in inp:
            c,a,b=line.rstrip().split('\t')[:3];out.write(f'{c}\t{a}\t{b}\t{c}:{a}-{b}\n')
    binary=shutil.which('bigWigAverageOverBed')
    if not binary:raise RuntimeError('Install the UCSC bigWigAverageOverBed utility and put it on PATH')
    def one(entry):
        source=ROOT/entry['bigwig_file']
        if not source.is_file():raise FileNotFoundError(f'Restore {source}; see input/bigwig_sources.tsv')
        table=raw/(entry['raw_binned_name']+'.average.tsv')
        with (raw/(entry['mark']+'.log')).open('w') as log:
            subprocess.run([binary,str(source),str(bed4),str(table)],stdout=log,stderr=subprocess.STDOUT,check=True)
        dest=raw/entry['raw_binned_name']
        with table.open() as inp,dest.open('w') as out:
            for line in inp:
                fields=line.split();chrom,coords=fields[0].split(':');a,b=coords.split('-')
                out.write(f'{chrom}\t{a}\t{b}\t{fields[5]}\n')
        table.unlink()
        return entry['mark'],dest
    reuse_ok = reuse_raw and (output/'preprocessing.json').exists()
    if reuse_ok:
        previous=json.loads((output/'preprocessing.json').read_text())
        reuse_ok=previous['fixed_grid_sha256']==cfg['fixed_grid_sha256'] and all((raw/x['raw_binned_name']).is_file() for x in cfg['signals'])
    if reuse_ok:
        paths={x['mark']:raw/x['raw_binned_name'] for x in cfg['signals']}
    else:
        with ThreadPoolExecutor(max_workers=workers) as pool:paths=dict(pool.map(one,cfg['signals']))
    input_path=paths['Input'];outputs=[]
    # The retained thyroid submission uses signal/(Input+0.01), not the old generic NA template.
    program=r"""BEGIN{OFS="\t"} {
      if ((getline other < inp) <= 0) {print "Input ended early" > "/dev/stderr"; exit 2}
      split(other,a,"\t");
      if ($1!=a[1] || $2!=a[2] || $3!=a[3]) {print "Coordinate mismatch" > "/dev/stderr"; exit 3}
      ratio=$4/(a[4]+0.01);
      print $1,$2,$3,ratio
    }"""
    for entry in cfg['signals']:
        if entry['mark']=='Input':continue
        if entry['mark']=='WGB-Seq':
            dest=paths[entry['mark']]
            ordered=dest.with_suffix('.sorted.bp')
            with ordered.open('w') as out:
                subprocess.run(['sort','-k1,1','-k2,2n',str(dest)],stdout=out,check=True,env={**os.environ,'LC_ALL':'C'})
            ordered.replace(dest)
        else:
            dest=normalized/Path(entry['binned_file']).name
            with dest.open('w') as out:
                subprocess.run(['awk','-v','inp='+str(input_path),program,str(paths[entry['mark']])],stdout=out,check=True)
        outputs.append(str(dest.resolve()))
    manifest=output/'signals.txt';manifest.write_text('\n'.join(outputs)+'\n')
    (output/'preprocessing.json').write_text(json.dumps(dict(
        fixed_grid_sha256=digest.hexdigest(),binning_executable=binary,n_tracks_binned=8,
        format='chr<TAB>start<TAB>end<TAB>signal; no header; original coordinate grid and order',
        histone_normalization='ChIP mean / (matching Input mean + 0.01); original thyroid awk formula',
        raw_binning_reused=bool(reuse_ok),
        methylation='Original supplied derived methylation bigWig, mean on the same windows; no Input division'),indent=2)+'\n')
    return manifest
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path)
    p.add_argument('--workers',type=int,default=8);a=p.parse_args();print(prepare(a.output,a.workers))
