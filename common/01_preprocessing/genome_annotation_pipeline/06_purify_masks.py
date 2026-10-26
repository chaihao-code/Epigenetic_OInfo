#!/usr/bin/env python3
"""06 按优先级重建 *_purified_bins.tsv（各层互斥）。

层级（每层剔除所有上层；promoter 保持原始"TSS bin"定义不变）：

    promoter > proximal enh. CTCF+ > proximal enh. CTCF-
    > distal enh. CTCF+ > distal enh. CTCF- > exon > intron

与旧版 purify_genomic_regions.py 的差异：
  - distal 两层额外剔除 proximal enhancer 重合；
  - exon 剔除全部 enhancer 层（旧版只剔除 promoter）；
  - intron 剔除全部 enhancer 层及 exon（旧版只剔除 promoter，且与 exon 可重叠）。

用法:
    python3 06_purify_masks.py --genome hg38 [--dry-run]
    python3 06_purify_masks.py --genome mm10 [--dry-run]

旧文件首次运行时备份为 *_purified_bins.tsv.bak。
"""

# Portable repository paths; no analysis is executed by this helper.
from pathlib import Path as _ReleasePath
_RELEASE_ROOT = next(p for p in _ReleasePath(__file__).resolve().parents
                     if (p / ".epigenetic_oinfo_root").is_file())
def _release_path(value):
    return value.replace("@ROOT@", str(_RELEASE_ROOT))
def _release_config(value):
    if isinstance(value, dict):
        return {k: _release_config(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_release_config(v) for v in value]
    return _release_path(value) if isinstance(value, str) else value

import argparse
import itertools
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(_release_path('@ROOT@'))
FEATURE_DIR = REPO / 'outputs/genome_features'

TIERS = [
    'promoter',
    'proximal_ctcf_bound_enhancer',
    'proximal_not_ctcf_bound_enhancer',
    'distal_ctcf_bound_enhancer',
    'distal_not_ctcf_bound_enhancer',
    'exon',
    'intron',
]


def load_mask(fdir: Path, name: str, index: pd.Index) -> np.ndarray:
    df = pd.read_csv(fdir / f'{name}_bins.tsv', sep='\t')
    col = df.columns[1]
    m = pd.Series(0, index=index, dtype=int)
    m.loc[df.loc[df[col] == 1, 'bin_id']] = 1
    return m.values.astype(bool)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--genome', required=True, choices=['hg38', 'mm10'])
    ap.add_argument('--dry-run', action='store_true')
    args = ap.parse_args()

    fdir = FEATURE_DIR / args.genome
    if not fdir.exists():
        sys.exit(f'missing feature dir: {fdir}')

    # bin universe（所有 mask 共用同一套 bin_id）
    ref = pd.read_csv(fdir / 'exon_bins.tsv', sep='\t', usecols=['bin_id'])
    bin_index = pd.Index(ref['bin_id'])
    print(f'{args.genome}: {len(bin_index)} bins in universe')

    masks = {}
    for t in TIERS:
        p = fdir / f'{t}_bins.tsv'
        if not p.exists():
            sys.exit(f'missing input mask: {p}')
        masks[t] = load_mask(fdir, t, bin_index)

    purified, taken = {}, np.zeros(len(bin_index), dtype=bool)
    for t in TIERS:
        purified[t] = masks[t] & ~taken
        taken |= purified[t]

    print(f"\n{'tier':40s} {'raw':>10s} {'purified':>10s}")
    for t in TIERS:
        print(f'{t:40s} {int(masks[t].sum()):>10d} {int(purified[t].sum()):>10d}')

    bad = [(a, b, int((purified[a] & purified[b]).sum()))
           for a, b in itertools.combinations(TIERS, 2)
           if (purified[a] & purified[b]).any()]
    if bad:
        for a, b, n in bad:
            print(f'validation FAILED: {a} ∩ {b} = {n}')
        sys.exit(1)
    print('validation: all purified tiers mutually exclusive')

    if args.dry_run:
        print('\n(dry-run, not writing)')
        return

    for t in TIERS:
        out = fdir / f'{t}_purified_bins.tsv'
        bak = fdir / f'{t}_purified_bins.tsv.bak'
        n_old = None
        if out.exists():
            old = pd.read_csv(out, sep='\t')
            n_old = int(old[old.columns[1]].sum())
            if not bak.exists():
                shutil.copy2(out, bak)
        col = f'{t}_purified_bins'
        pd.DataFrame({'bin_id': bin_index, col: purified[t].astype(int)}
                     ).to_csv(out, sep='\t', index=False)
        print(f'wrote {out.name}: {n_old} -> {int(purified[t].sum())} bins')


if __name__ == '__main__':
    main()
