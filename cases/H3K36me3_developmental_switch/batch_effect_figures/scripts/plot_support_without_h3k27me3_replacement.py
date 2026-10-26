#!/usr/bin/env python3
"""Drop variants using replacement H3K27me3, retain original H3K27me3 and BH56."""

# Resolve repository-relative paths in retained operational JSON manifests.
_RELEASE_ROOT = next(p for p in __import__('pathlib').Path(__file__).resolve().parents
                     if (p / '.epigenetic_oinfo_root').is_file())
def _release_manifest(value):
    if isinstance(value, dict):
        return {k: _release_manifest(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_release_manifest(v) for v in value]
    if isinstance(value, str) and value.startswith(('cases/', 'common/', 'data/')):
        return str(_RELEASE_ROOT / value)
    return value
import csv
import gzip
import json
import os
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean

os.environ.setdefault('MPLCONFIGDIR', '/tmp/bulk16-support-mpl')
os.environ.setdefault('XDG_CACHE_HOME', '/tmp/bulk16-support-cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'without_h3k27me3_replacement'
STAGES = ['2cell', '8cell', 'ICM']


def table(name, rows):
    with (OUT / name).open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def main():
    OUT.mkdir(exist_ok=True)
    manifest = _release_manifest(json.loads((ROOT / 'manifest_16_fdr.json').read_text()))
    retained = {s: {t['mask'] for t in manifest['tasks'] if t['stage'] == s
                    and 'H3K27me3' not in t['replaced_marks']} for s in STAGES}
    assert all(masks == {0, 1, 2, 3, 8, 9, 10, 11} for masks in retained.values())
    rows = []
    with gzip.open(ROOT / 'data/all_16variant_triplets_bh56.tsv.gz', 'rt') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            r['mask'] = int(r['mask'])
            r['q'] = float(r['q_down_block_bh56'])
            r['omega'] = float(r['observed_hoi'])
            rows.append(r)
    assert len(rows) == 48 * 40 * 56
    plt.rcParams.update({'font.size': 5, 'axes.labelsize': 5, 'axes.titlesize': 5,
                         'xtick.labelsize': 5, 'ytick.labelsize': 5, 'pdf.fonttype': 42,
                         'ps.fonttype': 42, 'axes.linewidth': .5})
    distributions, support, statistics, details = [], [], [], []
    for criterion, suffix in [('FDR <= 0.05', 'fdr05'), ('FDR < 0.1', 'fdr10')]:
        fig, axes = plt.subplots(1, 3, figsize=(18 / 2.54, 8 / 2.54), sharey=True)
        max_count = 0
        for stage, ax in zip(STAGES, axes):
            all_sets = {mask: set() for mask in range(16)}
            keys = set()
            for r in rows:
                if r['stage'] != stage:
                    continue
                key = r['block'], r['triplet']
                keys.add(key)
                if r['omega'] < 0 and (r['q'] <= .05 if suffix == 'fdr05' else r['q'] < .1):
                    all_sets[r['mask']].add(key)
            assert len(keys) == 2240
            subset = [all_sets[m] for m in sorted(retained[stage])]
            common, union = set.intersection(*subset), set.union(*subset)
            for n, collection in [(16, list(all_sets.values())), (8, subset)]:
                i, u = set.intersection(*collection), set.union(*collection)
                assert all(collection)
                statistics.append(dict(stage=stage, criterion=criterion, n_variants=n,
                                       n_intersection=len(i), n_union=len(u),
                                       intersection_over_union=len(i) / len(u),
                                       mean_intersection_fraction=mean(len(i) / len(s) for s in collection)))
            assert set.intersection(*all_sets.values()) <= common
            counts = Counter()
            for block, triplet in sorted(keys):
                count = sum((block, triplet) in s for s in subset)
                counts[count] += 1
                support.append(dict(stage=stage, criterion=criterion, block=block,
                                    triplet=triplet, n_significant_variants=count))
            assert counts[8] == len(common) and counts[0] == 2240 - len(union)
            for mask in sorted(retained[stage]):
                details.append(dict(stage=stage, criterion=criterion, mask=mask,
                                    n_set=len(all_sets[mask]), n_intersection=len(common),
                                    intersection_fraction=len(common) / len(all_sets[mask])))
            values = [counts[x] for x in range(9)]
            max_count = max(max_count, max(values))
            bars = ax.bar(range(9), values, width=.75,
                          color=['#bdbdbd'] + ['#4c78a8'] * 7 + ['#2a9d6f'])
            ax.bar_label(bars, labels=[str(v) for v in values], padding=2, fontsize=5)
            ax.set_title(stage)
            ax.set_xticks(range(9))
            ax.set_xlabel('Number of significant variants (out of 8)')
            ax.spines[['top', 'right']].set_visible(False)
            ax.tick_params(width=.5, length=2)
            for x in range(9):
                distributions.append(dict(stage=stage, criterion=criterion,
                                          n_significant_variants=x, n_region_triplet_records=counts[x],
                                          fraction=counts[x] / 2240))
        axes[0].set_ylim(0, max_count * 1.14)
        axes[0].set_ylabel('Number of region × triplet records')
        fig.suptitle('8 source combinations; original H3K27me3 retained | ' + criterion.replace('<=', '≤'), y=.97)
        fig.subplots_adjust(left=.08, right=.99, bottom=.19, top=.85, wspace=.15)
        for ext in ['pdf', 'png']:
            fig.savefig(OUT / f'triplet_support_distribution_{suffix}.{ext}', dpi=300)
        plt.close(fig)
    table('support_distribution.tsv', distributions)
    table('triplet_support_across_variants.tsv', support)
    table('stability_comparison.tsv', statistics)
    table('intersection_fraction_each_variant.tsv', details)
    (OUT / 'analysis_manifest.json').write_text(json.dumps(dict(
        excluded_source='replacement H3K27me3 (GSE73952)',
        original_H3K27me3_retained=True, retained_masks={s: sorted(m) for s, m in retained.items()},
        late_included=False, BH_family_size=56, n_permutations=10000,
        unit='stage × genomic annotation block × triplet',
        q_values='Existing within-variant, within-block BH56 values; no recalculation needed',
        caveat='Reducing 16 variants to 8 makes strict intersection easier; this alone does not establish data quality.',
        status='complete'), indent=2) + '\n')
    print(json.dumps(statistics, indent=2))


if __name__ == '__main__':
    main()
