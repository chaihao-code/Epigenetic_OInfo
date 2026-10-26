#!/usr/bin/env python3
"""Original-only developmental participation heatmap, with BH56 FDR <= .05."""

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
from collections import defaultdict
from pathlib import Path
from statistics import mean

os.environ.setdefault('MPLCONFIGDIR', '/tmp/bulk16-support-mpl')
os.environ.setdefault('XDG_CACHE_HOME', '/tmp/bulk16-support-cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import PercentFormatter

ROOT = Path(__file__).resolve().parents[1]
BATCH = ROOT.parent / 'batch_effect_figures'
EARLY = ['zygote', '2cell', '4cell', '8cell', 'morula', 'ICM', 'TE']
LATE = ['E11.5', 'E12.5', 'E13.5', 'E14.5', 'E15.5', 'E16.5', 'P0']
MARKS = ['ATAC', 'DNA_methylation', 'H3K9me3', 'H3K27me3', 'H3K36me3', 'H3K4me1', 'H3K4me3', 'H3K27ac']


def read(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write(path, rows):
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def main():
    manifest = _release_manifest(json.loads((ROOT / 'inputs/source_manifest.json').read_text()))
    tissues = manifest['late_reference_tissues']
    # Match the three original bars exactly, including their permutation realizations.
    rows = [r for r in read(BATCH / 'data/all_16variant_triplets_bh56.tsv.gz') if int(r['mask']) == 0]
    for r in rows:
        r['q'] = float(r['q_down_block_bh56'])
    present = {r['stage'] for r in rows}
    for r in read(ROOT / 'inputs/preimplantation_additional_stages_bh56.tsv.gz'):
        if r['stage'] in EARLY and r['stage'] not in present:
            r['q'] = float(r['q_down_sample_block_bh56'])
            rows.append(r)
    for r in read(BATCH / 'data/late_reference_triplets_bh56.tsv.gz'):
        r['q'] = float(r['q_down_block_bh56'])
        rows.append(r)
    assert len(rows) == (7 + 7 * 5) * 40 * 56
    samples = defaultdict(list)
    for r in rows:
        samples[r['stage'], r['sample']].append(r)
    sample_summary = []
    for (stage, sample), tests in samples.items():
        assert len(tests) == 2240 and len({(r['block'], r['triplet']) for r in tests}) == 2240
        for mark in MARKS:
            relevant = [r for r in tests if mark in r['triplet'].split('/')]
            assert len(relevant) == 40 * 21
            n = sum(float(r['observed_hoi']) < 0 and r['q'] <= .05 for r in relevant)
            sample_summary.append(dict(stage=stage, sample=sample, mark=mark,
                                       n_significant=n, n_possible=840, mean_fraction=n / 840))
    stages = EARLY + LATE
    summary = []
    for stage in stages:
        for mark in MARKS:
            z = [r for r in sample_summary if r['stage'] == stage and r['mark'] == mark]
            assert len(z) == (1 if stage in EARLY else 5)
            if stage in LATE:
                assert {r['sample'] for r in z} == {stage + '_' + t for t in tissues}
            summary.append(dict(stage=stage, mark=mark, mean_fraction=mean(r['mean_fraction'] for r in z),
                                n_samples=len(z)))
    original_reference = [r for r in read(BATCH / 'same_stage_mark_summary.tsv') if int(r['mask']) == 0]
    for r in original_reference:
        new = next(x for x in summary if x['stage'] == r['stage'] and x['mark'] == r['mark'])
        assert abs(new['mean_fraction'] - float(r['mean_fraction'])) < 1e-12
    # Check that averaging the seven late rows reproduces the prior Late reference.
    old_late = read(BATCH / 'data/late_common_tissue_means.tsv')
    for mark in MARKS:
        new = mean(r['mean_fraction'] for r in summary if r['stage'] in LATE and r['mark'] == mark)
        old = mean(float(r['mean_fraction']) for r in old_late if r['mark'] == mark)
        assert abs(new - old) < 1e-12
    write(ROOT / 'data/original_development_mark_participation.tsv', summary)
    write(ROOT / 'data/original_development_mark_participation_by_sample.tsv', sample_summary)
    matrix = np.array([[next(r['mean_fraction'] for r in summary if r['stage'] == s and r['mark'] == m)
                        for m in MARKS] for s in stages])
    plt.rcParams.update({'font.size': 5, 'axes.labelsize': 5, 'axes.titlesize': 5,
                         'xtick.labelsize': 5, 'ytick.labelsize': 5, 'pdf.fonttype': 42, 'ps.fonttype': 42})
    fig, ax = plt.subplots(figsize=(18 / 2.54, 12 / 2.54))
    im = ax.imshow(matrix, aspect='auto', cmap='YlOrRd', vmin=0, vmax=max(.5, float(matrix.max())))
    ax.set_xticks(range(8), ['ATAC', 'DNA methylation'] + MARKS[2:], rotation=40, ha='right')
    ax.set_yticks(range(14), ['Zygote', '2-cell', '4-cell', '8-cell', 'Morula', 'ICM', 'TE'] + LATE)
    ax.axhline(6.5, color='black', linewidth=.8)
    ax.set_ylabel('Developmental stage')
    ax.set_title('Synergistic triplet participation across development (original signals)')
    for i in range(14):
        for j in range(8):
            ax.text(j, i, f'{matrix[i, j]:.1%}', ha='center', va='center',
                    color='white' if matrix[i, j] / im.norm.vmax > .58 else 'black', fontsize=5)
    cb = fig.colorbar(im, ax=ax, fraction=.035, pad=.035)
    cb.ax.yaxis.set_major_formatter(PercentFormatter(1))
    cb.set_label('Significant synergistic triplets / 21 possible per signal')
    fig.subplots_adjust(left=.12, right=.88, bottom=.22, top=.91)
    for ext in ['pdf', 'png']:
        fig.savefig(ROOT / 'plots' / f'original_development_mark_participation_heatmap.{ext}', dpi=300)
    plt.close(fig)
    validation = dict(status='complete', criterion='observed O-information < 0 and BH56 q <= 0.05',
                      n_permutations=10000, n_signals=8, n_stages=14, n_blocks=40,
                      n_triplets_per_signal_per_block=21, late_tissues=tissues,
                      aggregation='Mean across 40 blocks per sample, then equal mean of five common tissues per late stage',
                      early_sources='Original mask 0 for 2cell/8cell/ICM; current mouse_preimp BH56 for other four stages',
                      late_source='Current mouse_after_E11.5 BH56 tests for five common tissues',
                      original_bar_values_match=True, mean_late_matches_original_Late_reference=True)
    (ROOT / 'data/original_development_heatmap_validation.json').write_text(json.dumps(validation, indent=2) + '\n')
    print(json.dumps(validation, indent=2))


if __name__ == '__main__':
    main()
