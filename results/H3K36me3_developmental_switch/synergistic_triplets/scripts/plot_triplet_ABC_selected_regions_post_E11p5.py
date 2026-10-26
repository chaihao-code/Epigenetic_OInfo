#!/usr/bin/env python3
"""Eight-signal late companion to the selected-region seven-signal figure."""
import csv
import gzip
import hashlib
import itertools
import json
import math
from collections import Counter, defaultdict
from statistics import mean, pstdev

from plot_triplet_A_B_C_post_E11p5 import (
    ROOT, SOURCE, STAGES, MARKS, REPRESSIVE, COLORS, classify, write_tsv, plt, np,
)
from matplotlib import font_manager

PANELS = [
    ('Exon / SINE', 'exon_purified_bins_AND_SINE'),
    ('Distal enhancer CTCF− / SINE', 'distal_not_ctcf_bound_enhancer_purified_bins_AND_SINE'),
    ('Proximal enhancer CTCF+ / CpG island', 'proximal_ctcf_bound_enhancer_purified_bins_AND_CpG_island'),
    ('Intron / SINE', 'intron_purified_bins_AND_SINE'),
    ('Distal enhancer CTCF− / LINE', 'distal_not_ctcf_bound_enhancer_purified_bins_AND_LINE'),
    ('None / None', 'None_AND_None'),
]
MAXIMA = {'A': 15, 'B': 31, 'C': 10}
PREFIX = 'triplet_ABC_selected_region_trajectories_8signals_post_E11p5_n5'


def main():
    tissues = json.loads((ROOT / 'inputs/source_manifest.json').read_text())['late_reference_tissues']
    universe = {frozenset(t) for t in itertools.combinations(MARKS, 3)}
    assert dict(Counter(classify(t) for t in universe)) == MAXIMA
    assert classify(['DNA_methylation', 'H3K36me3', 'ATAC']) == 'A'
    assert classify(['H3K36me3', 'H3K4me1', 'ATAC']) == 'B'
    blocks = {b for _, b in PANELS}
    families = defaultdict(list)
    with gzip.open(SOURCE, 'rt') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            assert row['stage'] in STAGES and row['tissue'] in tissues
            assert row['sample'] == row['stage'] + '_' + row['tissue']
            assert int(row['n_perm']) == 10000
            assert frozenset(row['triplet'].split('/')) in universe
            if row['block'] in blocks:
                families[row['sample'], row['block']].append(row)
    assert len(families) == 35 * len(PANELS)
    per_sample = []
    for (sample, block), family in sorted(families.items()):
        assert len(family) == 56 and {frozenset(r['triplet'].split('/')) for r in family} == universe
        ordered = sorted(family, key=lambda r: float(r['p_down']))
        q = np.minimum(1, np.minimum.accumulate((np.array([float(r['p_down']) for r in ordered])
                                                * 56 / np.arange(1, 57))[::-1])[::-1])
        assert np.allclose(q, [float(r['q_down_block_bh56']) for r in ordered], rtol=0, atol=1e-12)
        stage, tissue = sample.split('_', 1)
        for cls in 'ABC':
            selected = [r for r in family if classify(r['triplet'].split('/')) == cls]
            assert len(selected) == MAXIMA[cls]
            count = sum(float(r['observed_hoi']) < 0 and float(r['q_down_block_bh56']) <= .05 for r in selected)
            per_sample.append(dict(sample=sample, stage=stage, tissue=tissue, block=block, triplet_class=cls,
                                   count=count, theoretical_max=MAXIMA[cls], relative=count / MAXIMA[cls]))
    grouped = defaultdict(list)
    for r in per_sample:
        grouped[r['block'], r['triplet_class'], r['stage']].append(r)
    positions = {s: 19.5 if s == 'P0' else float(s[1:]) for s in STAGES}
    summary, lookup = [], {}
    for _, block in PANELS:
        for cls in 'ABC':
            for stage in STAGES:
                records = grouped[block, cls, stage]
                assert len(records) == 5 and {r['tissue'] for r in records} == set(tissues)
                values = [r['relative'] for r in records]
                std = pstdev(values)
                record = dict(block=block, triplet_class=cls, stage=stage, x_pos=positions[stage],
                              mean_count=mean(r['count'] for r in records), theoretical_max=MAXIMA[cls],
                              mean_relative=mean(values), std_relative=std, sem_relative=std / math.sqrt(5), n=5)
                summary.append(record)
                lookup[block, cls, stage] = record
    checked = 0
    with (ROOT / 'data/triplet_ABC_stacked_post_E11p5_n5_summary.tsv').open() as handle:
        for r in csv.DictReader(handle, delimiter='\t'):
            if r['block'] in blocks:
                new = lookup[r['block'], r['triplet_class'], r['stage']]
                assert abs(new['mean_count'] - float(r['mean_count'])) < 1e-12
                checked += 1
    assert checked == len(summary) == 7 * 3 * len(PANELS)
    write_tsv(ROOT / 'data' / (PREFIX + '_by_sample.tsv'), per_sample)
    write_tsv(ROOT / 'data' / (PREFIX + '_summary.tsv'), summary)
    arial = font_manager.findfont('Arial', fallback_to_default=False)
    plt.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'axes.titlesize': 8,
                         'axes.labelsize': 8, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
                         'legend.fontsize': 8, 'figure.titlesize': 8,
                         'pdf.fonttype': 42, 'ps.fonttype': 42})
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
    handles = []
    for panel, (ax, (title, block)) in enumerate(zip(axes, PANELS)):
        for cls in 'ABC':
            values = [lookup[block, cls, s] for s in STAGES]
            line = ax.errorbar([v['x_pos'] for v in values], [v['mean_relative'] for v in values],
                               yerr=[v['sem_relative'] for v in values], color=COLORS[cls],
                               marker='o', markersize=2.5, markeredgewidth=.4, linewidth=.9,
                               elinewidth=.5, capsize=1.2, capthick=.5, label=f'Triplet {cls} (/{MAXIMA[cls]})')
            if panel == 0:
                handles.append(line)
        ax.set_title(title, fontsize=8, pad=5)
        ax.set_ylim(-.04, 1.06)
        ax.set_yticks([0, .25, .5, .75, 1])
        ax.set_xticks([positions[s] for s in STAGES], STAGES, rotation=90, ha='center')
        ax.set_xlim(11.25, 19.75)
        ax.tick_params(length=2, pad=1.5, labelsize=8)
        ax.spines[['top', 'right']].set_visible(False)
        ax.axhline(0, color='#999999', linewidth=.4, alpha=.4)
        if all(lookup[block, cls, s]['mean_count'] == 0 for cls in 'ABC' for s in STAGES):
            ax.text(.5, .5, 'A = B = C = 0', transform=ax.transAxes,
                    ha='center', va='center', fontsize=8, color='#555555')
        if panel >= 3:
            ax.set_xlabel('Stage', fontsize=8, labelpad=3)
            ax.xaxis.set_label_coords(.5, -1.0 / panel_height_cm)
    fig.supylabel('Significant triplets / class total', fontsize=8, x=.15 / width_cm,
                  y=(bottom_cm + panel_height_cm + row_gap_cm / 2) / height_cm)
    fig.suptitle('Synergistic triplets at E11.5–P0: eight-signal panel including H3K36me3', fontsize=8, y=1 - .12 / height_cm)
    fig.text(.5, 1 - .62 / height_cm, '10,000 permutations; O-information < 0; BH56 q ≤ 0.05; mean ± SEM across five common tissues',
             ha='center', va='center', fontsize=8)
    legend = fig.legend(handles=handles, labels=[h.get_label() for h in handles],
                        loc='upper center', bbox_to_anchor=(.5, 1 - .90 / height_cm), ncol=3, frameon=False,
                        fontsize=8, handlelength=2, columnspacing=2.5, borderaxespad=0)
    fig.text(.5, .14 / height_cm, 'Five common tissues per stage (n = 5); repressive group: H3K9me3, H3K27me3 and H3K36me3.',
             ha='center', va='bottom', fontsize=8)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for ax in axes:
        for text in [ax.title, ax.xaxis.label, ax.yaxis.label, *ax.get_xticklabels(), *ax.get_yticklabels()]:
            if text.get_visible() and text.get_text():
                bbox = text.get_window_extent(renderer)
                assert bbox.x0 >= 0 and bbox.y0 >= 0 and bbox.x1 <= fig.bbox.width and bbox.y1 <= fig.bbox.height, text.get_text()
        ticks = [t for t in ax.get_xticklabels() if t.get_visible()]
        for i, first in enumerate(ticks):
            for second in ticks[i + 1:]:
                assert not first.get_window_extent(renderer).overlaps(second.get_window_extent(renderer))
        bbox = ax.get_window_extent(renderer)
        assert abs(bbox.height / bbox.width - panel_aspect) < 1e-10
    from matplotlib.text import Text
    for text in fig.findobj(Text):
        if text.get_visible() and text.get_text():
            assert text.get_fontsize() == 8, text.get_text()
            bbox = text.get_window_extent(renderer)
            assert bbox.x0 >= 0 and bbox.y0 >= 0 and bbox.x1 <= fig.bbox.width and bbox.y1 <= fig.bbox.height, text.get_text()
    bbox = legend.get_window_extent(renderer)
    assert bbox.x0 >= 0 and bbox.x1 <= fig.bbox.width
    for ext in ('png', 'pdf'):
        path = ROOT / 'plots' / (PREFIX + '.' + ext)
        fig.savefig(path, dpi=300)
        print('Saved', path, flush=True)
    plt.close(fig)
    metadata = dict(status='complete', source=str(SOURCE), source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
                    signals=MARKS, n_signals=8, H3K36me3_included=True, repressive_marks=sorted(REPRESSIVE),
                    theoretical_max=MAXIMA, stages=STAGES, tissues=tissues, n_samples=35,
                    panels=[dict(title=t, block=b) for t, b in PANELS],
                    region_selection='Same regions as seven-signal companion; distal CTCF-negative enhancer SINE/LINE contrast and None_AND_None',
                    all_zero_panels=[b for _, b in PANELS if all(lookup[b,c,s]['mean_count'] == 0 for c in 'ABC' for s in STAGES)],
                    criterion='observed_hoi < 0 and BH56 q <= 0.05', n_permutations=10000,
                    aggregation='Equal mean across five common tissues; population SD / sqrt(5) SEM, matching reference',
                    P0_plot_day=19.5, figure_size_cm=[width_cm, height_cm],
                    panel_size_cm=[panel_width_cm, panel_height_cm], panel_height_to_width=panel_aspect,
                    original_panel_aspect_preserved=False,
                    panel_height_relative_to_previous=height_scale, figure_height_relative_to_previous=height_scale,
                    panel_titles='Single line', shared_y_axis_title=True,
                    all_font_sizes_pt=8, font_family='Arial', font_file=arial,
                    layout='Two rows of three panels, A/B/C overlaid; shared x and y axes; same regions and style as seven-signal companion',
                    n_summary_rows=len(summary), n_sample_rows=len(per_sample),
                    complete_56_triplet_families_checked=True, BH56_recomputed_and_matched=True,
                    H3K36me3_classification_checked=True, all_means_match_stacked_count_table=True,
                    text_bounds_checked=True, x_tick_overlap_checked=True)
    (ROOT / 'data' / (PREFIX + '_metadata.json')).write_text(json.dumps(metadata, indent=2) + '\n')


if __name__ == '__main__':
    main()
