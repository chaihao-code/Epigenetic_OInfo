#!/usr/bin/env python3
"""Prepare separate portable source and supplementary-data tarballs; no analyses."""
import argparse
import hashlib
import json
import tarfile
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
def kind(p):
    r = p.relative_to(ROOT)
    if p.is_symlink() or any(x in r.parts for x in ['release_archives', '__pycache__']):
        return None
    if r.parts[0] == 'demo':
        if 'results' in r.parts:
            if 'preprocessing' in r.parts and p.suffix in ['.bp', '.bed']:
                return None  # reproducible intermediates; original full input copies are supplied
            return 'supplementary'
        if 'input' in r.parts and p.name.endswith('.gz'):
            return 'supplementary'
        return 'source'
    if r.parts[0] == 'supplementary_data':
        return 'supplementary'
    if r.parts[0] == 'data':
        if p.name.endswith('.gz'):
            return 'supplementary'
        return 'source' if p.name == 'README.md' else None
    if 'tables' in r.parts and any(x in r.parts for x in ['per_sample', 'per_class', 'local_per_class']):
        return 'supplementary'
    if r.parts[0] == 'external' or (r.parts[0] == 'tools' and 'bin' in r.parts):
        return None
    if r.parts[0] == 'cases' and p.suffix in ['.csv', '.tsv', '.bed', '.rds'] and p.stat().st_size > 1_000_000:
        return 'supplementary'
    return 'source'

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true',help='Write archives; otherwise print destinations only')
    args=parser.parse_args();out=ROOT/'release_archives'
    targets={k:out/f'Epigenetic_OInfo_{k}.tar.gz' for k in ['source','supplementary']}
    for k,p in targets.items():print(k,p)
    if not args.write:return
    out.mkdir(exist_ok=True)
    items={k:[] for k in targets}
    for p in sorted(ROOT.rglob('*')):
        if p.is_file():
            k=kind(p)
            if k:items[k].append(p)
    checks=[]
    for k,target in targets.items():
        if target.exists():raise SystemExit(f'Refusing to overwrite {target}')
        with tarfile.open(target,'w:gz',compresslevel=1,dereference=False) as archive:
            for p in items[k]:archive.add(p,arcname=str(p.relative_to(ROOT)),recursive=False)
        digest=hashlib.sha256()
        with target.open('rb') as h:
            for chunk in iter(lambda:h.read(1024*1024),b''):digest.update(chunk)
        checks.append(f'{digest.hexdigest()}  {target.name}')
        print(k,'files',len(items[k]),'bytes',target.stat().st_size)
    (out/'SHA256SUMS').write_text('\n'.join(checks)+'\n')
if __name__=='__main__':main()
