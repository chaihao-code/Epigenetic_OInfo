#!/usr/bin/env python3
"""Decompress supplied existing masks. Does not run annotation or analysis."""
from pathlib import Path
import gzip
import shutil
ROOT = Path(__file__).resolve().parents[1]
def main():
    n = 0
    for source in sorted((ROOT / 'data/genome_features').glob('*/*.tsv.gz')):
        dest = source.with_suffix('')
        if dest.exists():
            continue
        if dest.is_symlink():
            dest.unlink()
        temp = dest.with_name(dest.name + '.unpacking')
        with gzip.open(source, 'rb') as inp, temp.open('wb') as out:
            shutil.copyfileobj(inp, out)
        temp.replace(dest)
        n += 1
    print(f'Restored {n} masks; existing files were retained.')
if __name__ == '__main__': main()
