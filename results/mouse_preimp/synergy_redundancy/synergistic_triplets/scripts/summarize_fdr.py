#!/usr/bin/env python3
"""BH84 lower-tail FDR per stage and block, using 10000 permutations."""
import ast
import csv
import gzip
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STAGES = ['zygote', '2cell', '4cell', '8cell', 'morula', 'ICM', 'TE']


def bh(pvalues):
    order = sorted(range(len(pvalues)), key=pvalues.__getitem__)
    qvalues = [1.] * len(order)
    running = 1.
    for rank in range(len(order), 0, -1):
        i = order[rank - 1]
        running = min(running, pvalues[i] * len(order) / rank)
        qvalues[i] = running
    return qvalues


def write_table(path, rows, columns):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'wt', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def counts(rows):
    return dict(n_tests=len(rows), n_negative=sum(r['observed_hoi'] < 0 for r in rows),
                n_nominal_synergistic=sum(r['nominal_synergistic'] for r in rows),
                n_fdr_synergistic=sum(r['fdr_synergistic'] for r in rows),
                n_removed_by_fdr=sum(r['nominal_synergistic'] and not r['fdr_synergistic'] for r in rows))


def main():
    all_rows, summaries, stage_summaries, sources = [], [], [], []
    max_obs_diff = 0.
    for stage in STAGES:
        folder = ROOT / 'permutation_results' / stage
        meta = json.loads((folder / (stage + '_metadata.json')).read_text())
        assert meta['n_features'] == 9 and meta['n_permutations'] == 10000
        assert meta['permutation_level'] == 'all' and meta['permutation_scope'] == 'block'
        source = folder / (stage + '_null_summary.csv')
        with source.open() as handle:
            data = [r for r in csv.DictReader(handle) if int(r['order']) == 3]
        assert len(data) == 40 * 84, (stage, len(data))
        # Confirm the new permutation run uses the same observed triplet values
        # as the already completed nine-signal calculation.
        reference = ROOT.parent / 'results_n_bins_5' / stage / (stage + '_oinfo.csv')
        observed = {}
        with reference.open() as handle:
            for r in csv.DictReader(handle):
                if int(r['order']) == 3:
                    observed[(r['block'], tuple(ast.literal_eval(r['remaining_indices'])))] = float(r['hoi'])
        with (folder / (stage + '_subsets.csv')).open() as handle:
            indices = {int(r['subset_idx']): tuple(ast.literal_eval(r['subset_indices']))
                       for r in csv.DictReader(handle)}
        groups = defaultdict(list)
        for r in data:
            assert int(r['n_perm']) == 10000
            row = dict(stage=stage, block=r['block'],
                       triplet='/'.join(sorted(name.split('.', 1)[-1] for name in r['subset_names'].split('/'))),
                       subset_names=r['subset_names'], n_perm=10000,
                       observed_hoi=float(r['observed_hoi']), p_down=float(r['p_down']))
            for key in ('observed_hoi', 'p_down'):
                assert math.isfinite(row[key])
            assert 1 / 10001 - 1e-15 <= row['p_down'] <= 1
            diff = abs(row['observed_hoi'] - observed[(r['block'], indices[int(r['subset_idx'])])])
            max_obs_diff = max(max_obs_diff, diff)
            assert diff < 1e-9, (stage, r['block'], diff)
            groups[row['block']].append(row)
        assert len(groups) == 40
        stage_rows = []
        for block, rows in sorted(groups.items()):
            assert len(rows) == 84 and len({r['triplet'] for r in rows}) == 84
            for r, q in zip(rows, bh([r['p_down'] for r in rows])):
                r['q_down_stage_block_bh84'] = q
                r['nominal_synergistic'] = r['observed_hoi'] < 0 and r['p_down'] <= .05
                r['fdr_synergistic'] = r['observed_hoi'] < 0 and q <= .05
            summaries.append(dict(stage=stage, block=block, **counts(rows)))
            stage_rows.extend(rows)
        all_rows.extend(stage_rows)
        stage_summaries.append(dict(stage=stage, **counts(stage_rows)))
        sources.append(dict(stage=stage, path=str(source), sha256=hashlib.sha256(source.read_bytes()).hexdigest()))
    columns = list(all_rows[0])
    write_table(ROOT / 'all_triplets_fdr.tsv.gz', all_rows, columns)
    hits = [r for r in all_rows if r['fdr_synergistic']]
    write_table(ROOT / 'synergistic_triplets_fdr05.tsv', hits, columns)
    write_table(ROOT / 'summary_by_stage.tsv', stage_summaries, list(stage_summaries[0]))
    write_table(ROOT / 'summary_by_stage_region.tsv', summaries, list(summaries[0]))
    report = dict(status='complete', n_samples=7, n_blocks_per_sample=40,
                  triplets_per_BH_family=84, n_permutations=10000, alpha=.05,
                  criterion='observed_hoi < 0 and lower-tail BH q <= 0.05',
                  BH_scope='All 84 triplets within each stage x genomic block; no sign prefilter',
                  empirical_p='(1 + number(null <= observed)) / 10001',
                  observed_max_difference_vs_completed_run=max_obs_diff,
                  totals=counts(all_rows), sources=sources)
    (ROOT / 'fdr_validation.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
