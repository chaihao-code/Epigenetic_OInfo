#!/usr/bin/env python3
"""Plot eight-signal A/B/C triplet proportions for the five common late tissues."""

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
import hashlib
import itertools
import json
import math
import os
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean, pstdev

os.environ.setdefault('MPLCONFIGDIR', '/tmp/bulk16-support-mpl')
os.environ.setdefault('XDG_CACHE_HOME', '/tmp/bulk16-support-cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT.parent / 'batch_effect_figures/data/late_reference_triplets_bh56.tsv.gz'
STAGES = ['E11.5', 'E12.5', 'E13.5', 'E14.5', 'E15.5', 'E16.5', 'P0']
MARKS = ['ATAC', 'DNA_methylation', 'H3K9me3', 'H3K27me3', 'H3K36me3',
         'H3K4me1', 'H3K4me3', 'H3K27ac']
REPRESSIVE = {'H3K9me3', 'H3K27me3', 'H3K36me3'}
PANELS = [
    ('exon / SINE', 'exon_purified_bins_AND_SINE'),
    ('proximal enh. CTCF− / SINE', 'proximal_not_ctcf_bound_enhancer_purified_bins_AND_SINE'),
    ('proximal enh. CTCF− / CpG island', 'proximal_not_ctcf_bound_enhancer_purified_bins_AND_CpG_island'),
    ('distal enh. CTCF− / LINE', 'distal_not_ctcf_bound_enhancer_purified_bins_AND_LINE'),
]
COLORS = {'A': '#ca2a28', 'B': '#349939', 'C': '#ef7c21'}
PREFIX = 'triplet_A_B_C_counts_trajectories_post_E11p5_n5'


def classify(marks):
    repressive = bool(set(marks) & REPRESSIVE)
    return 'A' if 'DNA_methylation' in marks and repressive else 'B' if repressive else 'C'


def write_tsv(path, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def main():
    tissues = _release_manifest(json.loads((ROOT / 'inputs/source_manifest.json').read_text()))['late_reference_tissues']
    universe = {frozenset(t) for t in itertools.combinations(MARKS, 3)}
    maxima = dict(Counter(classify(t) for t in universe))
    assert maxima == {'A': 15, 'B': 31, 'C': 10}
    with gzip.open(SOURCE, 'rt') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    tests = defaultdict(list)
    for r in rows:
        assert r['stage'] in STAGES and r['tissue'] in tissues
        assert r['sample'] == r['stage'] + '_' + r['tissue']
        assert int(r['n_perm']) == 10000
        assert frozenset(r['triplet'].split('/')) in universe
        tests[r['sample'], r['block']].append(r)
    assert len(tests) == 35 * 40 and len(rows) == 35 * 40 * 56
    # Verify complete eight-signal families and the preserved BH56 correction.
    for family in tests.values():
        assert len(family) == 56
        assert {frozenset(r['triplet'].split('/')) for r in family} == universe
        ordered = sorted(family, key=lambda r: float(r['p_down']))
        expected = np.minimum(1., np.minimum.accumulate(
            (np.array([float(r['p_down']) for r in ordered]) * 56 / np.arange(1, 57))[::-1])[::-1])
        assert np.allclose(expected, [float(r['q_down_block_bh56']) for r in ordered], rtol=0, atol=1e-12)
    per_sample = []
    for stage in STAGES:
        for tissue in tissues:
            sample = stage + '_' + tissue
            for _, block in PANELS:
                family = tests[sample, block]
                for cls in 'ABC':
                    selected = [r for r in family if classify(r['triplet'].split('/')) == cls]
                    assert len(selected) == maxima[cls]
                    count = sum(float(r['observed_hoi']) < 0 and float(r['q_down_block_bh56']) <= .05
                                for r in selected)
                    per_sample.append(dict(sample=sample, stage=stage, tissue=tissue, block=block,
                                           triplet_class=cls, count=count, theoretical_max=maxima[cls],
                                           relative=count / maxima[cls]))
    summary, lookup = [], {}
    positions = {s: 19.5 if s == 'P0' else float(s[1:]) for s in STAGES}
    for cls in 'ABC':
        for _, block in PANELS:
            for stage in STAGES:
                values = [r for r in per_sample if (r['triplet_class'], r['block'], r['stage']) == (cls, block, stage)]
                assert len(values) == 5 and {r['tissue'] for r in values} == set(tissues)
                proportions = [r['relative'] for r in values]
                std = pstdev(proportions)
                record = dict(triplet_class=cls, block=block, stage=stage, x_pos=positions[stage],
                              mean_count=mean(r['count'] for r in values), theoretical_max=maxima[cls],
                              mean_relative=mean(proportions), std_relative=std,
                              sem_relative=std / math.sqrt(5), n=5)
                assert 0 <= record['mean_relative'] <= 1
                summary.append(record)
                lookup[cls, block, stage] = record
    write_tsv(ROOT / 'data' / (PREFIX + '_by_sample.tsv'), per_sample)
    write_tsv(ROOT / 'data' / (PREFIX + '_summary.tsv'), summary)
    plt.rcParams.update({'font.family': 'Arial', 'font.size': 7,
                         'pdf.fonttype': 42, 'ps.fonttype': 42})
    fig, axes = plt.subplots(3, 4, figsize=(25 / 2.54, 16 / 2.54), sharex=True, sharey=True)
    fig.subplots_adjust(left=.075, right=.99, top=.86, bottom=.12, wspace=.10, hspace=.19)
    for row, cls in enumerate('ABC'):
        color = COLORS[cls]
        for col, (title, block) in enumerate(PANELS):
            ax = axes[row, col]
            values = [lookup[cls, block, stage] for stage in STAGES]
            ax.errorbar([v['x_pos'] for v in values], [v['mean_relative'] for v in values],
                        yerr=[v['sem_relative'] for v in values], color=color, marker='o',
                        markersize=2.4, markeredgewidth=.4, linewidth=.8, elinewidth=.5, capsize=1.3)
            ax.set_ylim(-.04, 1.06)
            ax.set_yticks([0, .25, .5, .75, 1])
            ax.set_xticks(list(positions.values()), STAGES, rotation=60, ha='right')
            ax.set_xlim(11.25, 19.75)
            ax.tick_params(length=2, pad=1, labelsize=6)
            ax.axhline(0, color='gray', linewidth=.4, alpha=.3)
            ax.spines[['top', 'right']].set_visible(False)
            if row == 0:
                ax.set_title(title, fontsize=7.5, pad=7)
            if col == 0:
                ax.set_ylabel('Triplet ' + cls + '\noccurrence proportion', color=color, fontsize=7)
                ax.tick_params(axis='y', colors=color)
            if row == 2:
                ax.set_xlabel('Stage', fontsize=7)
    fig.suptitle('Significant synergistic triplets: E11.5–P0 relative-count trajectories\n'
                 '10,000 permutations; BH56 q ≤ 0.05 and O-info < 0; mean ± SEM across five common tissues\n'
                 'Eight-signal panel including H3K36me3: A / 15, B / 31, C / 10', fontsize=9, y=.98)
    for ext in ('pdf', 'png'):
        fig.savefig(ROOT / 'plots' / (PREFIX + '.' + ext), dpi=300, bbox_inches='tight', pad_inches=.03)
    plt.close(fig)
    metadata = dict(status='complete', source=str(SOURCE), source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
                    signals=MARKS, repressive_marks=sorted(REPRESSIVE),
                    class_definitions={'A': 'DNA methylation and at least one repressive mark',
                                       'B': 'No DNA methylation and at least one repressive mark',
                                       'C': 'No repressive mark'},
                    theoretical_max=maxima, stages=STAGES, tissues=tissues, n_samples=35,
                    n_permutations=10000, criterion='observed_hoi < 0 and BH56 q <= 0.05',
                    panels=[dict(title=t, block=b) for t, b in PANELS],
                    aggregation='Equal tissue mean of count / class maximum; population SD / sqrt(5) SEM, matching reference figure',
                    P0_plot_day=19.5, complete_BH56_families_checked=True, BH56_values_recomputed_and_matched=True,
                    n_summary_rows=len(summary), n_sample_rows=len(per_sample))
    (ROOT / 'data' / (PREFIX + '_validation.json')).write_text(json.dumps(metadata, indent=2) + '\n')
    print(json.dumps(metadata, indent=2))


if __name__ == '__main__':
    main()
