#!/usr/bin/env python3
"""Compare strict 16-variant stability at existing FDR <= .05 and FDR < .1."""
import csv
import gzip
import json
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'fdr_lt_0p1'
STAGES = ['2cell', '8cell', 'ICM']


def write(name, rows):
    with (OUT / name).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def main():
    OUT.mkdir(exist_ok=True)
    groups = defaultdict(list)
    with gzip.open(ROOT / 'data/all_16variant_triplets_bh56.tsv.gz', 'rt') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            assert row['stage'] in STAGES
            row['mask'] = int(row['mask'])
            row['q'] = float(row['q_down_block_bh56'])
            row['observed'] = float(row['observed_hoi'])
            groups[row['stage'], row['block'], row['triplet']].append(row)
    assert len(groups) == 3 * 40 * 56
    sets = {s: {threshold: [set() for _ in range(16)] for threshold in ['0.05', '0.1']} for s in STAGES}
    support = []
    for (stage, block, triplet), rows in sorted(groups.items()):
        assert len(rows) == 16 and {r['mask'] for r in rows} == set(range(16))
        counts = {}
        for threshold in ['0.05', '0.1']:
            selected = [r for r in rows if r['observed'] < 0 and
                        (r['q'] <= .05 if threshold == '0.05' else r['q'] < .1)]
            counts[threshold] = len(selected)
            for r in selected:
                sets[stage][threshold][r['mask']].add((block, triplet))
        support.append(dict(stage=stage, block=block, triplet=triplet,
                            n_variants_fdr05=counts['0.05'], n_variants_fdr10=counts['0.1'],
                            newly_in_intersection=counts['0.1'] == 16 and counts['0.05'] < 16,
                            q_max=max(r['q'] for r in rows), observed_hoi_max=max(r['observed'] for r in rows)))
    comparison = []
    for stage in STAGES:
        for threshold in ['0.05', '0.1']:
            variants = sets[stage][threshold]
            intersection, union = set.intersection(*variants), set.union(*variants)
            if threshold == '0.1':
                assert all(a <= b for a, b in zip(sets[stage]['0.05'], variants))
            comparison.append(dict(stage=stage, criterion='FDR <= 0.05' if threshold == '0.05' else 'FDR < 0.1',
                                   n_intersection=len(intersection), n_union=len(union),
                                   intersection_over_union=len(intersection) / len(union),
                                   n_distinct_triplets=len({t for b, t in intersection}),
                                   n_with_DNA_methylation=sum('DNA_methylation' in t.split('/') for b, t in intersection)))
    stable = [r for r in support if r['n_variants_fdr10'] == 16]
    by_triplet = defaultdict(list)
    for r in stable:
        by_triplet[r['stage'], r['triplet']].append(r)
    summary = [dict(stage=s, triplet=t, n_regions=len(rows),
                    n_new_regions=sum(r['newly_in_intersection'] for r in rows),
                    regions=';'.join(r['block'] for r in rows))
               for (s, t), rows in by_triplet.items()]
    summary.sort(key=lambda r: (STAGES.index(r['stage']), -r['n_regions'], r['triplet']))
    write('threshold_comparison.tsv', comparison)
    write('stable_intersection_triplets.tsv', stable)
    write('new_intersection_triplets.tsv', [r for r in stable if r['newly_in_intersection']])
    write('stable_intersection_triplet_summary.tsv', summary)
    write('triplet_support_across_variants.tsv', support)
    with (ROOT / 'stability_by_stage.tsv').open() as handle:
        for old in csv.DictReader(handle, delimiter='\t'):
            current = next(r for r in comparison if r['stage'] == old['stage'] and r['criterion'] == 'FDR <= 0.05')
            assert current['n_intersection'] == int(old['n_intersection'])
            assert current['n_union'] == int(old['n_union'])
    (OUT / 'README.md').write_text('Threshold sensitivity analysis using existing 10,000-permutation BH56 q-values.\n'
        'BH family: 56 triplets per stage × variant × genomic annotation block.\n'
        'Significance also requires observed O-information < 0. New threshold is strictly FDR < 0.1.\n'
        'Intersection and union use all 16 source variants, including original. Late is excluded.\n'
        'Original FDR <= 0.05 results and plots are preserved. No new permutations are needed.\n')
    print(json.dumps(comparison, indent=2))
    for stage in STAGES:
        print(stage, 'TOP', json.dumps([r for r in summary if r['stage'] == stage][:5]))


if __name__ == '__main__':
    main()
