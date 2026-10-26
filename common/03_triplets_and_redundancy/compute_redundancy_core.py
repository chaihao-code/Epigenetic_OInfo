#!/usr/bin/env python3
"""Identify redundancy cores using same-order single-signal replacement backgrounds.

For each block and order, test only the maximum-O-information subset. Require
positive O-information, modified z > threshold, and a one-component GMM selected
by BIC or replacement-identity ANOVA eta > threshold for a two-component fit.
Promote an accepted (N-1)-core to N if the full set has higher O-information;
remove accepted cores strictly contained in another accepted core.
This is a replacement-background criterion, not a permutation p-value/FDR test.
"""
import argparse
import csv
import hashlib
import json
import math
import re
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture
import statsmodels.api as sm
import statsmodels.formula.api as smf

ETA_THRESH = 0.5


def canonical_signal(token):
    match = re.search(r'H3K\d+(?:me\d|ac)|H2A\.Z|ATAC|DNA_methylation|CpG_methyl|WGBS', token)
    name = match[0] if match else token.strip()
    return {'CpG_methyl': 'DNA_methylation', 'WGBS': 'DNA_methylation'}.get(name, name)


def load_oinfo(path, signals, min_order):
    blocks = defaultdict(dict)
    with path.open() as f:
        for row in csv.DictReader(f):
            names = row.get('remaining_names', row.get('subset_names', ''))
            subset = frozenset(canonical_signal(x) for x in names.split('/') if x.strip())
            if not subset or not subset <= signals or len(subset) < min_order:
                continue
            if len(subset) != int(row['order']):
                raise ValueError(f'Duplicate/unrecognized signal names: {names}')
            value = float(row.get('hoi', row.get('observed_hoi', 'nan')))
            if not math.isfinite(value):
                raise ValueError(f'Non-finite O-information: {names}')
            block = row['block']
            if subset in blocks[block] and not math.isclose(blocks[block][subset], value, abs_tol=1e-12):
                raise ValueError(f'Inconsistent duplicate subset in {block}: {names}')
            blocks[block][subset] = value
    if not blocks:
        raise ValueError('No matching subsets found')
    expected = {frozenset(c) for k in range(min_order, len(signals) + 1)
                for c in combinations(sorted(signals), k)}
    for block, lookup in blocks.items():
        if set(lookup) != expected:
            raise ValueError(f'{block}: incomplete subsets, missing {len(expected - set(lookup))}')
    return blocks


def compute_bic_k1k2(values):
    values = np.asarray(values, dtype=float).reshape(-1, 1)
    n = len(values)
    if n < 2 or np.std(values, ddof=1) < 1e-12:
        return -1e12, 1e12
    try:
        gmm1 = GaussianMixture(n_components=1, random_state=0, max_iter=200, n_init=5)
        gmm1.fit(values)
        bic1 = gmm1.bic(values)
    except Exception:
        bic1 = 1e12
    try:
        gmm2 = GaussianMixture(n_components=2, random_state=0, max_iter=200, n_init=5)
        gmm2.fit(values)
        bic2 = gmm2.bic(values)
    except Exception:
        bic2 = 1e12
    return bic1, bic2


def partition_variance(df):
    try:
        if df['replaced'].nunique() < 2 or df['new'].nunique() < 2:
            return {'new': 0.0, 'replaced': 0.0}
        model = smf.ols('hoi ~ C(replaced) + C(new)', data=df).fit()
        if not np.isfinite(model.rsquared) or model.df_resid <= 0:
            return {'new': 0.0, 'replaced': 0.0}
        aov = sm.stats.anova_lm(model, typ=2)
        ss = aov['sum_sq']
        total = ss.sum()
        if total <= 0:
            return {'new': 0.0, 'replaced': 0.0}
        eta_new = ss.get('C(new)', 0.0) / total
        eta_replaced = ss.get('C(replaced)', 0.0) / total
        return {'new': float(eta_new), 'replaced': float(eta_replaced)}
    except Exception:
        return {'new': 0.0, 'replaced': 0.0}


def check_modality(repl_df, eta_thresh=ETA_THRESH):
    vals = repl_df['hoi'].values
    bic1, bic2 = compute_bic_k1k2(vals)
    best_k = 1 if bic1 <= bic2 else 2
    delta_bic = bic2 - bic1
    if best_k == 1:
        return True, {'n_components': 1, 'delta_bic': float(delta_bic), 'eta_new': np.nan}
    eta = partition_variance(repl_df)
    eta_new = eta.get('new', 0.0)
    pass_cluster = eta_new > eta_thresh
    return pass_cluster, {'n_components': 2, 'delta_bic': float(delta_bic), 'eta_new': float(eta_new)}


def modified_z_score(target_value, background):
    background = np.asarray(background, dtype=float)
    med = np.median(background)
    mad = np.median(np.abs(background - med)) + 1e-12
    return 0.6745 * (target_value - med) / mad


def replacement_background(target, full_signals, lookup):
    rows = []
    for old in sorted(target):
        for new in sorted(full_signals - target):
            rows.append({'replaced': old, 'new': new, 'hoi': lookup[target - {old} | {new}]})
    return pd.DataFrame(rows)


def find_redundancy_cores(lookup, signals, min_order=3, z_threshold=2.5, eta_threshold=.5):
    candidates, diagnostics = [], []
    n = len(signals)
    for order in range(min_order, n):
        subsets = sorted((s for s in lookup if len(s) == order), key=lambda s: tuple(sorted(s)))
        target = max(subsets, key=lambda s: lookup[s])
        value = lookup[target]
        background = replacement_background(target, signals, lookup)
        assert len(background) == order * (n - order)
        z = float(modified_z_score(value, background['hoi'].to_numpy()))
        passes_modality, info = check_modality(background, eta_thresh=eta_threshold)
        selected = value > 0 and z > z_threshold and passes_modality
        record = dict(order=order, signals='/'.join(sorted(target)), hoi=value, mod_z=z,
                      n_background=len(background), n_components=info['n_components'],
                      delta_bic=info['delta_bic'], eta_new=info['eta_new'],
                      positive=value > 0, passes_z=z > z_threshold,
                      passes_modality=passes_modality, selected_before_pruning=selected,
                      promoted=False)
        diagnostics.append(record.copy())
        if selected:
            candidates.append(record)
    full_value = lookup[signals]
    if any(c['order'] == n - 1 and full_value > c['hoi'] for c in candidates):
        candidates.append(dict(order=n, signals='/'.join(sorted(signals)), hoi=full_value,
                               mod_z=float('nan'), n_background=0, n_components=1,
                               delta_bic=float('nan'), eta_new=float('nan'), positive=True,
                               passes_z=True, passes_modality=True,
                               selected_before_pruning=True, promoted=True))
    kept = []
    for candidate in sorted(candidates, key=lambda c: -c['order']):
        subset = set(candidate['signals'].split('/'))
        if not any(subset < set(k['signals'].split('/')) for k in kept):
            kept.append(candidate)
    return sorted(kept, key=lambda c: c['order']), diagnostics


FIELDS = ['sample', 'block', 'order', 'signals', 'hoi', 'mod_z', 'n_background',
          'n_components', 'delta_bic', 'eta_new', 'positive', 'passes_z',
          'passes_modality', 'selected_before_pruning', 'promoted']


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True, type=Path, help='Combined *_oinfo.csv for one sample')
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--signals', required=True, nargs='+', help='Canonical signal names defining the comparison universe')
    parser.add_argument('--sample', help='Sample identifier (default: input stem without _oinfo)')
    parser.add_argument('--min-order', type=int, default=3)
    parser.add_argument('--z-threshold', type=float, default=2.5)
    parser.add_argument('--eta-threshold', type=float, default=.5)
    parser.add_argument('--expected-blocks', type=int)
    args = parser.parse_args()
    signals = frozenset(canonical_signal(x) for x in args.signals)
    if len(signals) != len(args.signals) or not 3 <= args.min_order < len(signals):
        parser.error('Signals must be unique and 3 <= min-order < number of signals')
    sample = args.sample or args.input.stem.removesuffix('_oinfo')
    blocks = load_oinfo(args.input, signals, args.min_order)
    if args.expected_blocks is not None and len(blocks) != args.expected_blocks:
        raise ValueError(f'Expected {args.expected_blocks} blocks, found {len(blocks)}')
    cores, diagnostics, block_counts = [], [], []
    for block, lookup in sorted(blocks.items()):
        selected, tested = find_redundancy_cores(lookup, signals, args.min_order,
                                                args.z_threshold, args.eta_threshold)
        cores.extend(dict(sample=sample, block=block, **r) for r in selected)
        diagnostics.extend(dict(sample=sample, block=block, **r) for r in tested)
        block_counts.append(dict(block=block, n_cores=len(selected)))
        print(f'{sample}: {block}: {len(selected)} core(s)', flush=True)
    write_csv(args.output, cores)
    write_csv(args.output.with_name(args.output.stem + '_diagnostics.csv'), diagnostics)
    import sklearn, statsmodels
    metadata = dict(status='complete', sample=sample, input=str(args.input.resolve()),
                    input_sha256=hashlib.sha256(args.input.read_bytes()).hexdigest(),
                    program_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    signals=sorted(signals), n_blocks=len(blocks), n_cores=len(cores),
                    min_order=args.min_order, z_threshold=args.z_threshold, eta_threshold=args.eta_threshold,
                    algorithm='same-order replacement background; GMM BIC; type-II ANOVA; maximal-set pruning',
                    core_FDR_applied=False, GMM_seed=0, GMM_n_init=5,
                    sklearn_version=sklearn.__version__, statsmodels_version=statsmodels.__version__,
                    block_counts=block_counts)
    args.output.with_suffix('.metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print(f'Completed {sample}: {len(cores)} redundancy cores in {len(blocks)} blocks', flush=True)


if __name__ == '__main__':
    main()
