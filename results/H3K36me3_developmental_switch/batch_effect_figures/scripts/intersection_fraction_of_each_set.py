#!/usr/bin/env python3
"""Average |intersection of all 16 variants| / |each variant's significant set|."""
import csv
import gzip
import json
from pathlib import Path
from statistics import mean

ROOT = Path(__file__).resolve().parents[1]
STAGES = ['2cell', '8cell', 'ICM']


def write(path, rows):
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def main():
    sets = {(s, threshold): [set() for _ in range(16)]
            for s in STAGES for threshold in ['FDR <= 0.05', 'FDR < 0.1']}
    labels = {}
    n_rows = 0
    with gzip.open(ROOT / 'data/all_16variant_triplets_bh56.tsv.gz', 'rt') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            n_rows += 1
            stage, mask = r['stage'], int(r['mask'])
            assert stage in STAGES and 0 <= mask < 16
            labels[stage, mask] = r['variant']
            if float(r['observed_hoi']) >= 0:
                continue
            q = float(r['q_down_block_bh56'])
            key = r['block'], r['triplet']
            if q <= .05:
                sets[stage, 'FDR <= 0.05'][mask].add(key)
            if q < .1:
                sets[stage, 'FDR < 0.1'][mask].add(key)
    assert n_rows == 48 * 40 * 56
    details, summary = [], []
    for (stage, criterion), variants in sets.items():
        intersection = set.intersection(*variants)
        assert all(variants), 'Empty significant set: ratio undefined'
        fractions = [len(intersection) / len(s) for s in variants]
        for mask, (s, fraction) in enumerate(zip(variants, fractions)):
            details.append(dict(stage=stage, criterion=criterion, mask=mask,
                                variant=labels[stage, mask], n_intersection=len(intersection),
                                n_set=len(s), intersection_fraction=fraction))
        summary.append(dict(stage=stage, criterion=criterion, n_variants=16,
                            n_intersection=len(intersection),
                            mean_intersection_fraction=mean(fractions),
                            min_intersection_fraction=min(fractions),
                            max_intersection_fraction=max(fractions),
                            min_set_size=min(map(len, variants)), max_set_size=max(map(len, variants))))
    write(ROOT / 'data/intersection_fraction_each_variant.tsv', details)
    write(ROOT / 'data/intersection_fraction_summary.tsv', summary)
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
