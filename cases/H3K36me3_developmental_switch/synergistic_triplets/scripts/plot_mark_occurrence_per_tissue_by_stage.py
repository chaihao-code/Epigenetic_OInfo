#!/usr/bin/env python3
"""Original-only mark occurrence in BH56-significant triplets, early and late panels.

Count each mark once for each significant region × triplet record across 40 blocks;
average sample totals over the five common tissues at each late stage. Eight shared
signals, observed O-information < 0 and FDR <= 0.05, 10,000 permutations.
Input is the validated per-sample table from plot_original_development_heatmap.py.
"""

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
import json
import os
from pathlib import Path
from statistics import mean

os.environ.setdefault('MPLCONFIGDIR', '/tmp/bulk16-support-mpl')
os.environ.setdefault('XDG_CACHE_HOME', '/tmp/bulk16-support-cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MARKS = ['DNA_methylation', 'H3K9me3', 'H3K27me3', 'ATAC', 'H3K4me1', 'H3K4me3', 'H3K27ac', 'H3K36me3']
LABELS = ['DNA methylation'] + MARKS[1:]
# Bright, moderately saturated cool colors for DNA methylation/H3K9me3/H3K27me3;
# warm colors for ATAC/H3K4me1/H3K4me3/H3K27ac; neutral gray for H3K36me3.
COLORS = ['#69A9D3', '#9B8ACB', '#53B5B3', '#F2BD61', '#EBA16C', '#E58078', '#D888AA', '#AAAAAA']
EARLY = ['zygote', '2cell', '4cell', '8cell', 'morula', 'ICM', 'TE']
LATE = ['E11.5', 'E12.5', 'E13.5', 'E14.5', 'E15.5', 'E16.5', 'P0']


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write(path, rows):
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def main():
    samples = read(ROOT / 'data/original_development_mark_participation_by_sample.tsv')
    fractions = read(ROOT / 'data/original_development_mark_participation.tsv')
    tissues = _release_manifest(json.loads((ROOT / 'inputs/source_manifest.json').read_text()))['late_reference_tissues']
    assert len(samples) == 42 * 8
    assert len({(r['sample'], r['mark']) for r in samples}) == len(samples)
    summary = []
    for stage in EARLY + LATE:
        row = dict(stage=stage, n_samples=1 if stage in EARLY else 5)
        for mark in MARKS:
            z = [r for r in samples if r['stage'] == stage and r['mark'] == mark]
            assert len(z) == row['n_samples']
            if stage in LATE:
                assert {r['sample'] for r in z} == {stage + '_' + t for t in tissues}
            assert all(int(r['n_possible']) == 840 for r in z)
            row[mark] = mean(int(r['n_significant']) for r in z)
            reference = next(float(r['mean_fraction']) for r in fractions if r['stage'] == stage and r['mark'] == mark)
            assert abs(row[mark] / 840 - reference) < 1e-12
        summary.append(row)
    arial_path = font_manager.findfont('Arial', fallback_to_default=False)
    plt.rcParams.update({'font.family': 'Arial', 'font.size': 8,
                         'axes.labelsize': 8, 'axes.titlesize': 10,
                         'xtick.labelsize': 8, 'ytick.labelsize': 8,
                         'legend.fontsize': 8, 'pdf.fonttype': 42, 'ps.fonttype': 42})
    for stages, period, title in [(EARLY, 'preimplantation', 'Preimplantation'),
                                   (LATE, 'post_E11p5', 'E11.5–P0')]:
        subset = [r for r in summary if r['stage'] in stages]
        write(ROOT / 'data' / f'mark_occurrence_per_tissue_by_stage_{period}.tsv', subset)
        fig, ax = plt.subplots(figsize=(10 / 2.54, 6 / 2.54))
        x = np.arange(len(stages))
        width = .095
        for i, (mark, label, color) in enumerate(zip(MARKS, LABELS, COLORS)):
            values = [next(r[mark] for r in subset if r['stage'] == s) for s in stages]
            ax.bar(x + (i - 3.5) * width, values, width=width, label=label,
                   color=color, edgecolor='#606060', linewidth=.3)
        ax.set_xticks(x, ['Zygote', '2-cell', '4-cell', '8-cell', 'Morula', 'ICM', 'TE'] if stages == EARLY else stages)
        ax.set_ylabel('Occurrence per sample' if stages == EARLY else 'Mean occurrence per tissue')
        ax.set_title(f'{title}: mark occurrence\nin synergistic triplets', pad=5)
        ax.set_ylim(0, 350)
        ax.set_yticks([0, 100, 200, 300])
        ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(width=.5, length=2)
        handles, labels = ax.get_legend_handles_labels()
        fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(.5, .005),
                   frameon=False, ncol=3, columnspacing=.7, handletextpad=.35,
                   handlelength=1.1, labelspacing=.35, borderaxespad=0)
        fig.subplots_adjust(left=.16, right=.99, top=.79, bottom=.38)
        for ext in ['pdf', 'png']:
            fig.savefig(ROOT / 'plots' / f'mark_occurrence_per_tissue_by_stage_{period}.{ext}', dpi=300)
        plt.close(fig)
    # Combined publication figure: two panels with one shared legend and y scale.
    fig, axes = plt.subplots(1, 2, figsize=(18 / 2.54, 6 / 2.54), sharey=True)
    for ax, stages, title in zip(axes, [EARLY, LATE], ['Preimplantation', 'E11.5–P0']):
        x = np.arange(len(stages))
        for i, (mark, label, color) in enumerate(zip(MARKS, LABELS, COLORS)):
            values = [next(r[mark] for r in summary if r['stage'] == s) for s in stages]
            ax.bar(x + (i - 3.5) * .095, values, width=.095, label=label,
                   color=color, edgecolor='#606060', linewidth=.3)
        ax.set_xticks(x, ['Zygote', '2-cell', '4-cell', '8-cell', 'Morula', 'ICM', 'TE'] if stages == EARLY else stages)
        ax.set_title(title, pad=5)
        ax.set_ylim(0, 350)
        ax.set_yticks([0, 100, 200, 300])
        ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(width=.5, length=2)
    axes[0].set_ylabel('Mean occurrence per sample')
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(.5, .005),
               frameon=False, ncol=4, columnspacing=1., handletextpad=.35,
               handlelength=1.1, labelspacing=.35, borderaxespad=0)
    fig.suptitle('Mark occurrence in synergistic triplets', fontsize=10, y=.98)
    fig.subplots_adjust(left=.085, right=.99, top=.76, bottom=.29, wspace=.12)
    for ext in ['pdf', 'png']:
        fig.savefig(ROOT / 'plots' / f'mark_occurrence_per_tissue_by_stage_combined.{ext}', dpi=300)
    plt.close(fig)
    write(ROOT / 'data/mark_occurrence_per_tissue_by_stage.tsv', summary)
    (ROOT / 'data/mark_occurrence_validation.json').write_text(json.dumps(dict(
        status='complete', n_signals=8, n_early_stages=7, n_late_stages=7,
        n_early_samples=7, n_late_samples=35, original_only=True, late_tissues=tissues,
        criterion='observed O-information < 0 and BH56 FDR <= 0.05',
        unit='Number of significant genomic annotation block × triplet records involving each mark',
        aggregation='Sum across 40 annotation blocks within sample; equal mean of five tissues per late stage',
        checked_against_heatmap=True, conversion_from_heatmap='occurrence = fraction × 40 × 21',
        figure_size_cm=[10, 6], combined_figure_size_cm=[18, 6], font_family='Arial', font_file=arial_path,
        font_size_pt=8, title_font_size_pt=10,
        palette=dict(zip(MARKS, COLORS))), indent=2) + '\n')
    print('Saved separate and combined occurrence figures (PDF/PNG), tables and validation JSON.')


if __name__ == '__main__':
    main()
