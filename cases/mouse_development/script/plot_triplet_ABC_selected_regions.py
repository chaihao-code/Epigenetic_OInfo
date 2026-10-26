#!/usr/bin/env python3
"""Overlay seven-signal A/B/C trajectories in six representative regions."""
import csv
import gzip
import hashlib
import itertools
import json
import math
import statistics
from collections import Counter, defaultdict

from plot_triplet_A_B_C_counts_trajectories import OUTPUT, STAGES, COLORS

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np

PANELS = [
    ('Exon / SINE', 'exon_purified_bins_AND_SINE'),
    ('Distal enhancer CTCF− / SINE', 'distal_not_ctcf_bound_enhancer_purified_bins_AND_SINE'),
    ('Proximal enhancer CTCF+ / CpG island', 'proximal_ctcf_bound_enhancer_purified_bins_AND_CpG_island'),
    ('Intron / SINE', 'intron_purified_bins_AND_SINE'),
    ('Distal enhancer CTCF− / LINE', 'distal_not_ctcf_bound_enhancer_purified_bins_AND_LINE'),
    ('None / None', 'None_AND_None'),
]
MARKS = {'ATAC', 'CpG_methyl', 'H3K9me3', 'H3K27me3', 'H3K4me1', 'H3K4me3', 'H3K27ac'}
REPRESSIVE = {'H3K9me3', 'H3K27me3'}
MAXIMA = {'A': 9, 'B': 16, 'C': 10}
PREFIX = 'triplet_ABC_selected_region_trajectories_7signals'


def classify(marks):
    repressive = bool(set(marks) & REPRESSIVE)
    return 'A' if 'CpG_methyl' in marks and repressive else 'B' if repressive else 'C'


def write_tsv(path, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def main():
    source = OUTPUT / 'all_triplets_fdr_10000.tsv.gz'
    universe = {frozenset(t) for t in itertools.combinations(MARKS, 3)}
    assert dict(Counter(classify(t) for t in universe)) == MAXIMA
    blocks = {b for _, b in PANELS}
    families = defaultdict(list)
    samples = defaultdict(set)
    with gzip.open(source, 'rt') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            marks = frozenset(row['triplet'].split('/'))
            assert marks in universe and int(row['order']) == 3 and int(row['n_perm']) == 10000
            assert row['triplet_class'] == classify(marks)
            assert row['stage'] == row['sample'].split('_')[0]
            samples[row['stage']].add(row['sample'])
            if row['block'] in blocks:
                families[row['sample'], row['block']].append(row)
    assert set(samples) == set(STAGES) and sum(map(len, samples.values())) == 82
    assert len(families) == 82 * len(PANELS)
    per_sample = []
    for (sample, block), family in sorted(families.items()):
        assert len(family) == 35 and {frozenset(r['triplet'].split('/')) for r in family} == universe
        ordered = sorted(family, key=lambda r: float(r['p_down']))
        q = np.minimum(1, np.minimum.accumulate((np.array([float(r['p_down']) for r in ordered])
                                                * 35 / np.arange(1, 36))[::-1])[::-1])
        assert np.allclose(q, [float(r['q_down_sample_block_bh']) for r in ordered], rtol=0, atol=1e-12)
        for cls in 'ABC':
            tests = [r for r in family if r['triplet_class'] == cls]
            assert len(tests) == MAXIMA[cls]
            count = sum(float(r['observed_hoi']) < 0 and float(r['q_down_sample_block_bh']) <= .05 for r in tests)
            per_sample.append(dict(sample=sample, stage=sample.split('_')[0], block=block,
                                   triplet_class=cls, count=count, theoretical_max=MAXIMA[cls],
                                   relative=count / MAXIMA[cls]))
    positions = {s: i for i, s in enumerate(STAGES[:10])}
    # Keep day-based late spacing, with one day occupying one early category step.
    positions.update({s: 10 + ((19.5 if s == 'P0' else float(s[1:])) - 11.5) for s in STAGES[10:]})
    grouped = defaultdict(list)
    for row in per_sample:
        grouped[row['block'], row['triplet_class'], row['stage']].append(row)
    summary, lookup = [], {}
    for _, block in PANELS:
        for cls in 'ABC':
            for stage in STAGES:
                records = grouped[block, cls, stage]
                assert {r['sample'] for r in records} == samples[stage]
                values = [r['relative'] for r in records]
                n = len(values)
                std = statistics.pstdev(values)
                record = dict(block=block, triplet_class=cls, stage=stage, x_pos=positions[stage],
                              mean_count=statistics.mean(r['count'] for r in records), theoretical_max=MAXIMA[cls],
                              mean_relative=statistics.mean(values), std_relative=std,
                              sem_relative=std / math.sqrt(n), n=n)
                summary.append(record)
                lookup[block, cls, stage] = record
    # Check every displayed mean against the preserved all-region count tables.
    for cls in 'ABC':
        with (OUTPUT / f'triplet_{cls}_synergistic_count_by_stage.tsv').open() as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                if row['region'] in blocks:
                    new = lookup[row['region'], cls, row['stage']]
                    assert abs(new['mean_count'] - float(row['mean_fdr_significant_negative_count'])) < 1e-12
                    assert new['n'] == int(row['n_samples'])
    write_tsv(OUTPUT / (PREFIX + '_by_sample.tsv'), per_sample)
    write_tsv(OUTPUT / (PREFIX + '_summary.tsv'), summary)
    arial = font_manager.findfont('Arial', fallback_to_default=False)
    plt.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'axes.titlesize': 8,
                         'axes.labelsize': 8, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
                         'legend.fontsize': 8, 'figure.titlesize': 8,
                         'pdf.fonttype': 42, 'ps.fonttype': 42})
    # Compress both the overall figure and each panel to 3/4 of the previous height.
    old_panel_width_cm = 28 * (.99 - .065) / (4 + 3 * .13)
    original_panel_aspect = (9 * (.74 - .25)) / old_panel_width_cm
    height_scale = .75
    panel_aspect = original_panel_aspect * height_scale
    width_cm, left_cm, right_cm, top_cm = 18, 1.25, .20, 1.75
    row_gap_cm = .55
    panel_width_cm = (width_cm - left_cm - right_cm) / (3 + 2 * .13)
    panel_height_cm = panel_width_cm * panel_aspect
    previous_height_cm = 2.05 + 2 * panel_width_cm * original_panel_aspect + 1.0 + 2.25
    height_cm = previous_height_cm * height_scale
    bottom_cm = height_cm - 2 * panel_height_cm - row_gap_cm - top_cm
    fig, axes = plt.subplots(2, 3, figsize=(width_cm / 2.54, height_cm / 2.54), sharex=True, sharey=True)
    axes = axes.ravel()
    fig.subplots_adjust(left=left_cm / width_cm, right=1 - right_cm / width_cm,
                        bottom=bottom_cm / height_cm, top=1 - top_cm / height_cm,
                        wspace=.13, hspace=row_gap_cm / panel_height_cm)
    display_stages = ['Zygote', '2-cell', '4-cell', '8-cell', 'Morula', 'ICM', 'TE'] + STAGES[7:]
    handles = []
    for panel, (ax, (title, block)) in enumerate(zip(axes, PANELS)):
        for cls in 'ABC':
            values = [lookup[block, cls, s] for s in STAGES]
            line = ax.errorbar([v['x_pos'] for v in values], [v['mean_relative'] for v in values],
                               yerr=[v['sem_relative'] for v in values], color=COLORS[cls],
                               marker='o', markersize=2.5, markeredgewidth=.4, linewidth=.9,
                               elinewidth=.5, capsize=1.2, capthick=.5,
                               label=f'Triplet {cls} (/{MAXIMA[cls]})')
            if panel == 0:
                handles.append(line)
        ax.set_title(title, fontsize=8, pad=5)
        ax.set_ylim(-.04, 1.06)
        ax.set_yticks([0, .25, .5, .75, 1])
        ax.set_xticks([positions[s] for s in STAGES], display_stages, rotation=90, ha='center')
        ax.set_xlim(-.4, positions['P0'] + .4)
        ax.tick_params(length=2, pad=1.5, labelsize=8)
        ax.spines[['top', 'right']].set_visible(False)
        ax.axhline(0, color='#999999', linewidth=.4, alpha=.4)
        if all(lookup[block, cls, s]['mean_count'] == 0 for cls in 'ABC' for s in STAGES):
            ax.text(.5, .5, 'A = B = C = 0', transform=ax.transAxes,
                    ha='center', va='center', fontsize=8, color='#555555')
        for split in (6.5, 9.5):
            ax.axvline(split, color='#888888', linestyle='--', linewidth=.6, alpha=.65)
        if panel >= 3:
            ax.set_xlabel('Stage', fontsize=8, labelpad=3)
            ax.xaxis.set_label_coords(.5, -1.0 / panel_height_cm)
    fig.supylabel('Significant triplets / class total', fontsize=8, x=.15 / width_cm,
                  y=(bottom_cm + panel_height_cm + row_gap_cm / 2) / height_cm)
    fig.suptitle('Synergistic triplets across development: seven-signal panel', fontsize=8, y=1 - .12 / height_cm)
    fig.text(.5, 1 - .62 / height_cm, '10,000 permutations; O-information < 0; BH35 q ≤ 0.05; mean ± SEM across available tissues',
             ha='center', va='center', fontsize=8)
    legend = fig.legend(handles=handles, labels=[h.get_label() for h in handles],
                        loc='upper center', bbox_to_anchor=(.5, 1 - .90 / height_cm), ncol=3, frameon=False,
                        fontsize=8, handlelength=2, columnspacing=2.5, borderaxespad=0)
    fig.text(.5, .14 / height_cm, 'n = 1 (preimplantation), 2–4 (gastrulation), 8–12 (E11.5–P0); dashed lines separate periods.',
             ha='center', va='bottom', fontsize=8)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for ax in axes:
        for label in [ax.title, ax.xaxis.label, ax.yaxis.label, *ax.get_xticklabels(), *ax.get_yticklabels()]:
            if label.get_visible() and label.get_text():
                bbox = label.get_window_extent(renderer)
                assert bbox.x0 >= 0 and bbox.y0 >= 0 and bbox.x1 <= fig.bbox.width and bbox.y1 <= fig.bbox.height, label.get_text()
        ticks = [t for t in ax.get_xticklabels() if t.get_visible()]
        for i, first in enumerate(ticks):
            for second in ticks[i + 1:]:
                assert not first.get_window_extent(renderer).overlaps(second.get_window_extent(renderer)), (first.get_text(), second.get_text())
        bbox = ax.get_window_extent(renderer)
        assert abs(bbox.height / bbox.width - panel_aspect) < 1e-10
    from matplotlib.text import Text
    for label in fig.findobj(Text):
        if label.get_visible() and label.get_text():
            assert label.get_fontsize() == 8, label.get_text()
            bbox = label.get_window_extent(renderer)
            assert bbox.x0 >= 0 and bbox.y0 >= 0 and bbox.x1 <= fig.bbox.width and bbox.y1 <= fig.bbox.height, label.get_text()
    legend_bbox = legend.get_window_extent(renderer)
    assert legend_bbox.x0 >= 0 and legend_bbox.x1 <= fig.bbox.width
    for ext in ('png', 'pdf'):
        path = OUTPUT / (PREFIX + '.' + ext)
        fig.savefig(path, dpi=300)
        print('Saved', path, flush=True)
    plt.close(fig)
    metadata = dict(status='complete', source=str(source), source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                    signals=sorted(MARKS), repressive_marks=sorted(REPRESSIVE), n_signals=7,
                    H3K36me3_included=False, theoretical_max=MAXIMA,
                    class_definitions={'A': 'DNA methylation and at least one H3K9me3/H3K27me3',
                                       'B': 'At least one H3K9me3/H3K27me3, without DNA methylation',
                                       'C': 'No H3K9me3 or H3K27me3'},
                    panels=[dict(title=t, block=b) for t, b in PANELS], stages=STAGES,
                    region_selection='Exon/intron SINE, distal CTCF-negative enhancer SINE/LINE, proximal CTCF-positive enhancer CpG island, and None_AND_None',
                    all_zero_panels=[b for _, b in PANELS if all(lookup[b,c,s]['mean_count'] == 0 for c in 'ABC' for s in STAGES)],
                    n_samples=82, samples_per_stage={s:len(samples[s]) for s in STAGES},
                    sample_selection='All available samples within each stage, matching the original seven-signal figure',
                    criterion='observed_hoi < 0 and BH35 q <= 0.05', n_permutations=10000,
                    aggregation='Arithmetic mean of count / class maximum; population SD / sqrt(n) SEM, matching original figure',
                    x_spacing='One unit per early category and per late day; P0 at day 19.5',
                    layout='Two rows of three panels, A/B/C overlaid; shared x and y axes; bottom-row stage labels',
                    figure_size_cm=[width_cm, height_cm], panel_size_cm=[panel_width_cm, panel_height_cm],
                    panel_height_to_width=panel_aspect, original_panel_aspect_preserved=False,
                    panel_height_relative_to_previous=height_scale, figure_height_relative_to_previous=height_scale,
                    panel_titles='Single line', shared_y_axis_title=True,
                    all_font_sizes_pt=8, early_tick_layout='Single vertical lane; all 17 stage labels retained',
                    font_family='Arial', font_file=arial,
                    n_summary_rows=len(summary), n_sample_rows=len(per_sample),
                    complete_35_triplet_families_checked=True, BH35_recomputed_and_matched=True,
                    all_means_match_preserved_count_tables=True, text_bounds_checked=True, x_tick_overlap_checked=True)
    (OUTPUT / (PREFIX + '_metadata.json')).write_text(json.dumps(metadata, indent=2) + '\n')


if __name__ == '__main__':
    main()
