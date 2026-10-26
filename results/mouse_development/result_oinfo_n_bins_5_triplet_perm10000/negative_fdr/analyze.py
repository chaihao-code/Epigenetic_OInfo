#!/usr/bin/env python3
"""Snapshot completed triplet summaries and BH-correct every lower-tail test."""
import csv
import gzip
import hashlib
import json
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


def bh(values):
    n = len(values)
    order = sorted(range(n), key=values.__getitem__)
    adjusted = [1.0] * n
    running = 1.0
    for rank in range(n, 0, -1):
        idx = order[rank - 1]
        running = min(running, values[idx] * n / rank)
        adjusted[idx] = running
    return adjusted


def write_table(path, rows, columns):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'wt', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def main():
    output = Path(__file__).resolve().parent
    root = output.parent
    manifest = list(csv.DictReader((root / 'input_manifest.tsv').open(), delimiter='\t'))
    rows, sources, pending = [], [], []
    for entry in manifest:
        sample = entry['sample']
        source = root / sample / f'{sample}_null_summary.csv'
        metadata = root / sample / f'{sample}_metadata.json'
        if not source.exists() or not metadata.exists():
            pending.append(sample)
            continue
        payload = source.read_bytes()
        frame = list(csv.DictReader(payload.decode().splitlines()))
        meta = json.loads(metadata.read_text())
        assert meta['n_permutations'] == 10000, sample
        assert len(frame) == 1400, sample
        assert len({(x['block'], x['subset_idx']) for x in frame}) == 1400, sample
        assert len({x['block'] for x in frame}) == 40, sample
        for row in frame:
            assert int(row['order']) == 3 and int(row['n_perm']) == 10000, sample
            for key in ('p_down', 'p_up', 'p_dev', 'observed_hoi', 'null_mean'):
                row[key] = float(row[key])
            assert all(0 <= row[k] <= 1 for k in ('p_down', 'p_up', 'p_dev'))
            row['sample'] = sample
            signals = []
            for name in row['subset_names'].split('/'):
                match = re.search(r'ATAC|H3K27ac|H3K27me3|H3K4me1|H3K4me3|H3K9me3|CpG_methyl|WGBS', name)
                assert match, name
                signals.append(match[0].replace('WGBS', 'CpG_methyl'))
            row['triplet'] = '/'.join(signals)
            row['delta_from_null_mean'] = row['observed_hoi'] - row['null_mean']
            row['negative_raw_p05'] = row['observed_hoi'] < 0 and row['p_down'] <= .05
            rows.append(row)
        sources.append({'sample': sample, 'source': str(source), 'sha256': hashlib.sha256(payload).hexdigest(), 'n_tests': len(frame)})
    assert rows, 'No completed samples'
    families = defaultdict(list)
    for i, row in enumerate(rows):
        families[(row['sample'], row['block'])].append(i)
    for name, groups in (('q_down_sample_block_bh', list(families.values())),):
        for indices in groups:
            assert len(indices) == 35
            qvalues = bh([rows[i]['p_down'] for i in indices])
            for i, q in zip(indices, qvalues):
                rows[i][name] = q
    qcolumns = [k for k in rows[0] if k.startswith('q_down_')]
    rows.sort(key=lambda r: (r['q_down_sample_block_bh'], r['p_down'], r['observed_hoi']))
    columns = list(rows[0])
    write_table(output / 'all_triplets_fdr.tsv.gz', rows, columns)
    negative = [r for r in rows if r['observed_hoi'] < 0]
    write_table(output / 'negative_raw_p05.tsv', [r for r in negative if r['p_down'] <= .05], columns)
    for qcolumn in qcolumns:
        write_table(output / f'negative_{qcolumn}_05.tsv', [r for r in negative if r[qcolumn] <= .05], columns)

    def summarize(subset):
        neg = [r for r in subset if r['observed_hoi'] < 0]
        return dict(n_tests=len(subset), n_negative=len(neg), n_negative_raw_p05=sum(r['p_down'] <= .05 for r in neg),
                    **{k.replace('q_down_', 'n_negative_') + '_05': sum(r[k] <= .05 for r in neg) for k in qcolumns})

    summary = summarize(rows)
    for group in ('sample', 'block', 'triplet', 'sample_block'):
        grouped = defaultdict(list)
        for row in rows:
            key = (row['sample'], row['block']) if group == 'sample_block' else row[group]
            grouped[key].append(row)
        results = [dict(**({'sample': key[0], 'block': key[1]} if group == 'sample_block' else {group: key}), **summarize(value)) for key, value in grouped.items()]
        results.sort(key=lambda r: -r['n_negative_sample_block_bh_05'])
        write_table(output / f'summary_by_{group}.tsv', results, list(results[0]))
    write_table(output / 'input_snapshot.tsv', sources, list(sources[0]))
    metadata = dict(created_at_utc=datetime.now(timezone.utc).isoformat(), n_samples_complete=len(sources),
                    n_samples_expected=len(manifest), pending_samples=pending, alpha=.05,
                    method='Benjamini-Hochberg', pvalue='p_down = (1 + count(null <= observed))/(10000 + 1)',
                    signal_labels='WGBS and CpG_methyl grouped as CpG_methyl in triplet column; original subset_names retained.',
                    selection='Correct all lower-tail p-values first, then require observed_hoi < 0 and q <= 0.05.',
                    families={'sample_block': '35 lower-tail tests per sample/block; no pooling across samples, blocks, or tails'},
                    dependence='Ordinary BH; overlapping triplets do not have an arbitrary-dependence guarantee.',
                    summary=summary)
    (output / 'analysis_metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print(json.dumps(metadata, indent=2), flush=True)
    print('TOP SAMPLE/BLOCK BH35 SURVIVORS', flush=True)
    for row in [r for r in negative if r['q_down_sample_block_bh'] <= .05][:12]:
        print({k: row[k] for k in ('sample', 'block', 'triplet', 'observed_hoi', 'null_mean', 'p_down', 'q_down_sample_block_bh')}, flush=True)


if __name__ == '__main__':
    main()
