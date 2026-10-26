#!/usr/bin/env python3
"""Plot counts of region–triplet records significant in 0–16 source variants."""
import csv
import json
import os
from collections import Counter
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/bulk16-support-mpl')
os.environ.setdefault('XDG_CACHE_HOME', '/tmp/bulk16-support-cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
STAGES = ['2cell', '8cell', 'ICM']


def main():
    with (ROOT / 'fdr_lt_0p1/triplet_support_across_variants.tsv').open() as f:
        records = list(csv.DictReader(f, delimiter='\t'))
    assert len(records) == 3 * 40 * 56
    assert len({(r['stage'], r['block'], r['triplet']) for r in records}) == len(records)
    with (ROOT / 'fdr_lt_0p1/threshold_comparison.tsv').open() as f:
        reference = {(r['stage'], r['criterion']): r for r in csv.DictReader(f, delimiter='\t')}
    plt.rcParams.update({'font.size': 5, 'axes.labelsize': 5, 'axes.titlesize': 5,
                         'xtick.labelsize': 5, 'ytick.labelsize': 5, 'pdf.fonttype': 42,
                         'ps.fonttype': 42, 'axes.linewidth': .5})
    output = []
    for field, label, suffix in [('n_variants_fdr05', 'FDR ≤ 0.05', 'fdr05'),
                                 ('n_variants_fdr10', 'FDR < 0.1', 'fdr10')]:
        criterion = 'FDR <= 0.05' if suffix == 'fdr05' else 'FDR < 0.1'
        fig, axes = plt.subplots(1, 3, figsize=(18 / 2.54, 8 / 2.54), sharey=True)
        max_count = 0
        for ax, stage in zip(axes, STAGES):
            stage_records = [r for r in records if r['stage'] == stage]
            assert len(stage_records) == 2240
            counts = Counter(int(r[field]) for r in stage_records)
            assert all(0 <= x <= 16 for x in counts)
            ref = reference[stage, criterion]
            assert counts[16] == int(ref['n_intersection'])
            assert 2240 - counts[0] == int(ref['n_union'])
            values = [counts[x] for x in range(17)]
            max_count = max(max_count, max(values))
            bars = ax.bar(range(17), values, width=.78,
                          color=['#bdbdbd'] + ['#4c78a8'] * 15 + ['#2a9d6f'])
            ax.bar_label(bars, labels=[str(v) for v in values], padding=2, fontsize=5)
            ax.set_title(stage)
            ax.set_xticks(range(17))
            ax.set_xlim(-.8, 16.8)
            ax.set_xlabel('Number of significant variants (out of 16)')
            ax.spines[['top', 'right']].set_visible(False)
            ax.tick_params(width=.5, length=2)
            for x in range(17):
                output.append(dict(stage=stage, criterion=criterion, n_significant_variants=x,
                                   n_region_triplet_records=counts[x], fraction=counts[x] / 2240))
        axes[0].set_ylabel('Number of region × triplet records')
        axes[0].set_ylim(0, max_count * 1.14)
        fig.suptitle(f'Support across 16 source combinations | {label}', y=.97)
        fig.subplots_adjust(left=.08, right=.99, bottom=.19, top=.85, wspace=.15)
        for ext in ['pdf', 'png']:
            fig.savefig(ROOT / 'plots' / f'triplet_support_distribution_{suffix}.{ext}', dpi=300)
        plt.close(fig)
    with (ROOT / 'data/triplet_support_distribution.tsv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(output[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(output)
    (ROOT / 'data/triplet_support_distribution_validation.json').write_text(json.dumps(dict(
        status='complete', n_records_per_stage=2240, n_genomic_annotation_blocks=40,
        n_possible_triplets=56, n_variants=16, late_included=False,
        unit='stage × genomic annotation block × canonical triplet',
        checks=['unique keys', 'support range 0–16', '0 bin agrees with union',
                '16 bin agrees with intersection']), indent=2) + '\n')
    print('Saved two PDF/PNG figures and the 102-row distribution table.')


if __name__ == '__main__':
    main()
