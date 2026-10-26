#!/usr/bin/env python3
"""Recalculate BH35 for the complete 10,000-permutation run."""
import csv
import gzip
import hashlib
import json
import math
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
OUT = BASE / 'synergistic_triplets'
PCOLS = ('p_down', 'p_up', 'p_dev')


def bh(p):
    q = [1.] * len(p)
    running = 1.
    for rank, i in reversed(list(enumerate(sorted(range(len(p)), key=p.__getitem__), 1))):
        running = min(running, p[i] * len(p) / rank)
        q[i] = running
    return q


def write(path, rows, columns=None):
    columns = columns or list(rows[0])
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'wt', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def key(row):
    return row['sample'], row['block'], row['triplet']


def load_run(nperm, require_all=False):
    root = BASE / ('result_oinfo_n_bins_5_triplet_perm' + str(nperm))
    manifest = list(csv.DictReader((root / 'input_manifest.tsv').open(), delimiter='\t'))
    assert len(manifest) == 82
    results, inputs, missing = [], [], []
    for entry in manifest:
        sample = entry['sample']
        source = root / sample / (sample + '_null_summary.csv')
        metadata = root / sample / (sample + '_metadata.json')
        if not source.exists() or not metadata.exists():
            missing.append(sample)
            continue
        assert json.loads(metadata.read_text())['n_permutations'] == nperm
        payload = source.read_bytes()
        rows = list(csv.DictReader(payload.decode().splitlines()))
        assert len(rows) == 1400, sample
        groups = defaultdict(list)
        for row in rows:
            assert int(row['order']) == 3 and int(row['n_perm']) == nperm
            for column in (*PCOLS, 'observed_hoi', 'null_mean', 'null_sd', 'null_median'):
                row[column] = float(row[column])
                assert math.isfinite(row[column]), (sample, column)
            assert all(0 <= row[c] <= 1 for c in PCOLS)
            signals = []
            for name in row['subset_names'].split('/'):
                match = re.search(r'ATAC|H3K27ac|H3K27me3|H3K4me1|H3K4me3|H3K9me3|CpG_methyl|WGBS|K27ac|K4M3|K4me1', name)
                assert match, name
                aliases = {'WGBS': 'CpG_methyl', 'K27ac': 'H3K27ac', 'K4M3': 'H3K4me3', 'K4me1': 'H3K4me1'}
                signals.append(aliases.get(match[0], match[0]))
            assert len(set(signals)) == 3
            row['triplet'] = '/'.join(sorted(signals))
            row['sample'] = sample
            row['stage'] = sample.split('_')[0]
            methyl = 'CpG_methyl' in signals
            repressive = any(x in signals for x in ('H3K9me3', 'H3K27me3', 'H3K36me3'))
            row['triplet_class'] = 'A' if methyl and repressive else 'B' if repressive else 'C'
            groups[row['block']].append(row)
        assert len(groups) == 40
        for group in groups.values():
            assert len(group) == 35 and len({x['triplet'] for x in group}) == 35
            for row, q in zip(group, bh([x['p_down'] for x in group])):
                row['q_down_sample_block_bh'] = q
                row['raw_significant_negative'] = row['observed_hoi'] < 0 and row['p_down'] <= .05
                row['significant_negative'] = row['observed_hoi'] < 0 and q <= .05
                results.append(row)
        inputs.append(dict(n_permutations=nperm, sample=sample, source=str(source),
                           sha256=hashlib.sha256(payload).hexdigest()))
    if require_all and missing:
        raise RuntimeError('Incomplete samples: ' + ', '.join(missing))
    assert len({key(r) for r in results}) == len(results)
    return results, inputs, missing


def counts(rows):
    return dict(n_tests=len(rows), n_negative=sum(r['observed_hoi'] < 0 for r in rows),
                n_raw_significant_negative=sum(r['raw_significant_negative'] for r in rows),
                n_bh35_significant_negative=sum(r['significant_negative'] for r in rows))


def main():
    ten, sources10, missing10 = load_run(10000, require_all=True)
    OUT.mkdir(exist_ok=True)
    ten.sort(key=lambda r: (r['q_down_sample_block_bh'], r['p_down'], key(r)))
    write(OUT / 'all_triplets_fdr_10000.tsv.gz', ten)
    selected = [r for r in ten if r['significant_negative']]
    write(OUT / 'significant_negative_triplets_10000.tsv', selected, list(ten[0]))
    for fields in [('sample',), ('stage',), ('block',), ('triplet',), ('triplet_class',), ('sample', 'block')]:
        groups = defaultdict(list)
        for row in ten:
            groups[tuple(row[f] for f in fields)].append(row)
        summaries = [dict(zip(fields, k), **counts(v)) for k, v in groups.items()]
        summaries.sort(key=lambda x: -x['n_bh35_significant_negative'])
        write(OUT / ('summary_by_' + '_'.join(fields) + '.tsv'), summaries)
    write(OUT / 'input_snapshot.tsv', sources10)
    metadata = dict(created_at_utc=datetime.now(timezone.utc).isoformat(),
                    method='Benjamini-Hochberg', family='35 p_down values per sample x annotation block',
                    alpha=.05, selection='observed_hoi < 0 and q_down_sample_block_bh <= .05',
                    n_permutations=10000, n_samples=len(sources10), missing_10k_samples=missing10,
                    summary=counts(ten))
    (OUT / 'analysis_metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
    print(json.dumps(metadata, indent=2), flush=True)


if __name__ == '__main__':
    main()
