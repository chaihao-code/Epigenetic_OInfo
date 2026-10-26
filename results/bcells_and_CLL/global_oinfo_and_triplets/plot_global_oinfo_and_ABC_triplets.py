#!/usr/bin/env python3
"""Align the order-7 O-info lollipops and the existing FDR-filtered ABC means."""
import hashlib
import json
from pathlib import Path
import os
import sys

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT / 'scripts'))
from plot_abc_pooled_combined import BLOCKS, CAPACITIES
from plot_abc_triplet_classes import REGION_LABELS

import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.text import Text
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
PLOTS = HERE / 'plots'
DATA = HERE / 'data'
PREFIX = 'global_oinfo_and_ABC_triplets_18x8p5cm'
GROUPS = {
    'R': ['CEMT0004', 'CEMT0005', 'CEMT0025', 'CEMT0026',
          'CEMT0027', 'CEMT0028', 'CEMT0029', 'CEMT0030'],
    'S': ['CEMT0006', 'CEMT0094', 'CEMT0095', 'CEMT0096', 'CEMT0097'],
}
COLORS = {'Higher in R': '#E74C3C', 'Higher in S': '#3498DB',
          'Not significant': '#95A5A6'}
LEGEND_LABELS = {'Higher in R': 'Higher in CLL-R', 'Higher in S': 'Higher in CLL-S',
                 'Not significant': 'Not significant'}
ROW_LABELS = ['CLL-R', 'CLL-S', 'pB']


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_block(block):
    region, feature = block.split('_AND_')
    region = region.replace('_purified_bins', '')
    if region == 'None_func':
        region = 'None'
    if feature == 'None_rep':
        feature = 'None'
    return region + '_AND_' + feature


def global_oinfo_statistics():
    reference = json.loads((HERE / 'reference_config.json').read_text())
    assert reference['GROUP_DEFS'] == {'Group2': GROUPS['R'], 'Group3': GROUPS['S']}
    axis_reference = reference
    axis_blocks = [canonical_block(r + '_AND_' + f)
                   for r in axis_reference['REGION_ORDER'] for f in axis_reference['FEATURE_ORDER']]
    assert axis_blocks == [canonical_block(b) for b in BLOCKS]
    mapping = pd.read_csv(PROJECT / 'ABC_triplet_classes/sample_labels.tsv', sep='\t')
    for group, members in GROUPS.items():
        assert set(mapping.loc[mapping.pooled_group.eq(group), 'sample_id']) == set(members)
    frames, sources = [], {}
    for group, members in GROUPS.items():
        for member in members:
            sample = 'CEEHRC.' + member
            root = PROJECT / 'results_oinfo_n_bins_5' / sample
            source = root / (sample + '_oinfo.csv')
            df = pd.read_csv(source, usecols=['block', 'order', 'hoi', 'remaining_names'])
            df = df[df.order.eq(7)].copy()
            assert len(df) == 40 and set(df.block) == set(BLOCKS)
            assert not df.block.duplicated().any() and np.isfinite(df.hoi).all()
            for names in df.remaining_names:
                signals = names.split('/')
                assert len(signals) == 7
                assert all(sum(mark in name for name in signals) == 1 for mark in
                           ['H3K4me1', 'H3K4me3', 'H3K9me3', 'H3K27me3', 'H3K36me3', 'H3K27ac', 'WGB-Seq'])
            check_source = root / (sample + '_permutation_summary.csv')
            reference_values = pd.read_csv(check_source)
            check = df.merge(reference_values[['block', 'observed_hoi']], on='block', validate='one_to_one')
            assert len(check) == 40
            assert np.allclose(check.hoi, check.observed_hoi, rtol=0, atol=1e-12)
            df['sample'], df['group'] = sample, group
            frames.append(df.drop(columns='remaining_names'))
            sources[str(source.relative_to(PROJECT))] = sha256(source)
            sources[str(check_source.relative_to(PROJECT))] = sha256(check_source)
    samples = pd.concat(frames, ignore_index=True)
    matrix = samples.pivot(index='sample', columns='block', values='hoi').reindex(columns=BLOCKS)
    r = matrix.reindex(['CEEHRC.' + s for s in GROUPS['R']]).to_numpy()
    s = matrix.reindex(['CEEHRC.' + s for s in GROUPS['S']]).to_numpy()
    assert r.shape == (8, 40) and s.shape == (5, 40)
    test = stats.ttest_ind(r, s, axis=0, equal_var=False)
    mean_diff = r.mean(axis=0) - s.mean(axis=0)
    # Preserve exactly the directional transformation in the supplied lollipop script.
    p = np.where(mean_diff > 0, test.pvalue / 2, 1 - test.pvalue / 2)
    independent = stats.ttest_ind(r, s, axis=0, equal_var=False, alternative='greater')
    assert np.allclose(p, independent.pvalue, atol=1e-14, rtol=1e-12)
    assert np.isfinite(p).all() and ((p > 0) & (p <= 1)).all()
    summary = pd.DataFrame({'block': BLOCKS, 'order': 7, 'n_R': 8, 'n_S': 5,
                            'mean_R': r.mean(axis=0), 'mean_S': s.mean(axis=0),
                            'mean_diff_R_minus_S': mean_diff, 'welch_t': test.statistic,
                            'p_two_sided': test.pvalue, 'p_value': p,
                            'log10_p_value': np.log10(p), 'category': 'Not significant'})
    summary.loc[(p < .05) & (mean_diff > 0), 'category'] = 'Higher in R'
    summary.loc[(p < .05) & (mean_diff < 0), 'category'] = 'Higher in S'
    return samples, summary, sources


def render(summary, counts):
    font = Path(os.environ.get('BCELLS_ARIAL_FONT', str(Path.home() / '.local/share/fonts/arial.ttf')))
    if font.exists():
        font_manager.fontManager.addfont(str(font))
    plt.rcParams.update({'font.family': 'Arial' if font.exists() else 'DejaVu Sans',
                         'font.size': 8, 'axes.labelsize': 8, 'axes.titlesize': 8,
                         'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 8,
                         'pdf.fonttype': 42, 'ps.fonttype': 42, 'axes.linewidth': .6})
    width, height = 18., 8.5
    left, plot_width = 2.40, 14.05
    heat_bottom, heat_height, class_height = 2.65, 2.4, .8
    top_bottom, top_height = 5.50, 2.45
    fig = plt.figure(figsize=(width/2.54, height/2.54), dpi=300)
    top = fig.add_axes([left/width, top_bottom/height, plot_width/width, top_height/height])
    for i, row in enumerate(summary.itertuples()):
        color = COLORS[row.category]
        top.vlines(i, 0, row.log10_p_value, color=color, lw=.5, alpha=.5, zorder=1)
        top.scatter(i, row.log10_p_value, c=color, s=20, edgecolor='black', lw=.3, alpha=.9, zorder=3)
    top.axhline(np.log10(.05), color='gray', ls='--', lw=.5, alpha=.7)
    top.text(39.1, np.log10(.05)-.10, 'p = 0.05', ha='right', va='bottom', fontsize=8, color='gray')
    for k in range(1, 8):
        top.axvline(k*5-.5, color='gray', lw=.5, alpha=.6)
    top.set_xlim(-.5, 39.5)
    # True log10(p) values: negative ticks, with more significant points higher.
    assert summary.log10_p_value.min() >= -8.5
    top.set_ylim(0, -8.5)
    top.set_yticks([0, -2, -4, -6, -8])
    top.set_ylabel('log$_{10}$(p value)', fontsize=8, labelpad=4)
    top.set_xticks([])
    top.tick_params(axis='y', length=2.5, width=.6, pad=2)
    top.spines[['top', 'right']].set_visible(False)
    handles = [Line2D([], [], linestyle='', marker='o', markersize=4,
                      markerfacecolor=color, markeredgecolor='black', markeredgewidth=.3, label=LEGEND_LABELS[name])
               for name, color in COLORS.items()]
    fig.legend(handles=handles, loc='center', bbox_to_anchor=((left+plot_width/2)/width, 8.25/height),
               ncol=3, frameon=False, handletextpad=.45, columnspacing=1.5, borderpad=0)
    axes, cells, colorbars, matrices = [], [], [], {}
    for panel, cls in enumerate('ABC'):
        bottom = heat_bottom + (2-panel)*class_height
        sub = counts[counts['class'].eq(cls)]
        matrix = sub.pivot(index='pooled_group', columns='block', values='mean_count').reindex(
            index=['R', 'S', 'pB'], columns=BLOCKS).to_numpy()
        assert matrix.shape == (3, 40) and np.isfinite(matrix).all()
        assert ((matrix >= 0) & (matrix <= CAPACITIES[cls])).all()
        matrices[cls] = matrix.tolist()
        ax = fig.add_axes([left/width, bottom/height, plot_width/width, class_height/height], sharex=top)
        axes.append(ax)
        im = ax.imshow(matrix, cmap='YlOrRd', vmin=0, vmax=CAPACITIES[cls],
                       aspect='auto', interpolation='nearest')
        ax.set_xlim(-.5, 39.5)
        ax.set_yticks(range(3), ROW_LABELS)
        ax.tick_params(axis='y', length=0, pad=2, labelsize=8)
        ax.set_xticks([])
        for k in range(1, 8):
            ax.axvline(k*5-.5, color='gray', lw=.5)
        for i in range(3):
            for j in range(40):
                cells.append(ax.text(j, i, f'{matrix[i,j]:.1f}', ha='center', va='center', fontsize=3.5,
                                     color='white' if matrix[i,j] >= CAPACITIES[cls]*.62 else 'black'))
        fig.text(.67/width, (bottom+class_height/2)/height, 'Triplet '+cls,
                 ha='center', va='center', fontsize=8)
        cax = fig.add_axes([16.65/width, (bottom+.15)/height, .13/width, .50/height])
        cb = fig.colorbar(im, cax=cax, ticks=[0, CAPACITIES[cls]])
        cb.ax.tick_params(labelsize=8, length=1.5, pad=1.5)
        colorbars.append(cb)
    fig.text(17.70/width, (heat_bottom+heat_height/2)/height, 'mean occurence',
             rotation=90, ha='center', va='center', fontsize=8)
    bottom_ax = axes[-1]
    region_texts, feature_texts = [], []
    for k, label in enumerate(REGION_LABELS):
        arrow_y = -.15/class_height
        bottom_ax.annotate('', xy=(k*5+4.4, arrow_y), xytext=(k*5-.4, arrow_y),
                           xycoords=bottom_ax.get_xaxis_transform(),
                           arrowprops=dict(arrowstyle='<->', color='gray', lw=.5), annotation_clip=False)
        region_texts.append(bottom_ax.text(k*5+2, -.28/class_height, label,
                                          transform=bottom_ax.get_xaxis_transform(),
                                          ha='center', va='top', fontsize=8, clip_on=False))
    for j in range(40):
        feature_texts.append(bottom_ax.text(j, -1.08/class_height,
                            ['CpG island', 'LINE', 'SINE', 'LTR', 'None'][j%5],
                            transform=bottom_ax.get_xaxis_transform(), ha='center', va='top',
                            rotation=90, fontsize=8, clip_on=False))
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    cell_set = set(cells)
    for text in fig.findobj(Text):
        if text.get_visible() and text.get_text():
            assert text.get_fontsize() == (3.5 if text in cell_set else 8), text.get_text()
            box = text.get_window_extent(renderer)
            assert box.x0 >= 0 and box.y0 >= 0 and box.x1 <= fig.bbox.width and box.y1 <= fig.bbox.height, text.get_text()
            if text in cell_set:
                panel_box = text.axes.get_window_extent(renderer)
                assert box.width < panel_box.width/40 and box.height < panel_box.height/3
    row_labels = [t for ax in axes for t in ax.get_yticklabels()]
    for first, second in zip(row_labels, row_labels[1:]):
        assert not first.get_window_extent(renderer).overlaps(second.get_window_extent(renderer))
    for first, second in zip(region_texts, region_texts[1:]):
        assert not first.get_window_extent(renderer).overlaps(second.get_window_extent(renderer))
    assert min(t.get_window_extent(renderer).y0 for t in region_texts) > max(
        t.get_window_extent(renderer).y1 for t in feature_texts)
    reference_x = top.transData.transform(np.column_stack([np.arange(40), np.zeros(40)]))[:, 0]
    for ax in axes:
        mapped = ax.transData.transform(np.column_stack([np.arange(40), np.zeros(40)]))[:, 0]
        assert np.allclose(mapped, reference_x, rtol=0, atol=1e-8)
    for ext in ['pdf', 'png']:
        path = PLOTS / (PREFIX + '.' + ext)
        fig.savefig(path, dpi=300)
        print(path)
    plt.close(fig)
    return dict(figure_size_cm=[width, height], label_font_size_pt=8, cell_font_size_pt=3.5,
                shared_block_order=BLOCKS, shared_x_positions_verified=True,
                complete_x_labels_shown_once=True, text_bounds_checked=True,
                row_labels_do_not_overlap=True, cell_numbers_fit_cells=True,
                heatmap_height_cm=heat_height, heatmap_matrix_values=matrices,
                heatmap_color_limits={c:[0,CAPACITIES[c]] for c in 'ABC'},
                heatmap_colorbar_label='mean occurence', heatmap_display_row_labels=ROW_LABELS,
                lollipop_legend_labels=list(LEGEND_LABELS.values()),
                heatmap_numeric_transform='None: original mean counts and color scales retained',
                lollipop_y_values='log10(p_value)', lollipop_y_label='log10(p value)',
                lollipop_y_axis_limits=[0,-8.5], lollipop_more_significant_points='higher')


def main():
    for path in [PLOTS, DATA]:
        path.mkdir(exist_ok=True)
    abc_plots = PROJECT / 'ABC_triplet_classes/plots'
    source = abc_plots / 'ABC_pooled_average_counts.tsv'
    retained = [source, abc_plots / 'ABC_pooled_average_heatmap.pdf', abc_plots / 'ABC_pooled_average_heatmap.png']
    before = {str(p): sha256(p) for p in retained}
    counts = pd.read_csv(source, sep='\t')
    assert len(counts) == 360 and not counts.duplicated(['class','pooled_group','block']).any()
    assert counts.n_samples.eq(counts.pooled_group.map({'R':8,'S':5,'pB':4})).all()
    samples, summary, sources = global_oinfo_statistics()
    layout = render(summary, counts)
    samples.to_csv(DATA / 'global_oinfo_order7_by_sample.tsv', sep='\t', index=False)
    summary.to_csv(DATA / 'global_oinfo_order7_R_vs_S_by_region.tsv', sep='\t', index=False)
    for path, digest in before.items():
        assert sha256(Path(path)) == digest
    validation = dict(status='complete', statistics_source='Existing raw order-7 O-information; no genomic recomputation',
                      group_mapping={'R':'Group2','S':'Group3'}, group_members=GROUPS,
                      global_test='Welch t-test, same as supplied plot_volcano_g2_vs_g3.py',
                      p_value='P(R > S); p_two_sided/2 if mean_R>mean_S, else 1-p_two_sided/2',
                      p_value_matches_supplied_lollipop=True, independently_verified_with_scipy_alternative_greater=True,
                      global_significance='Nominal p < 0.05, unchanged from supplied lollipop; no BH change',
                      global_Oinfo_matches_observed_system_permutation_summaries=True,
                      global_source_sha256=sources, ABC_source=str(source), ABC_source_sha256=before[str(source)],
                      ABC_criterion='observed_hoi < 0 AND q_down_global_BH < 0.05; 23,800-test global BH family',
                      ABC_original_outputs_and_counts_unchanged=True,
                      n_regions=40, n_global_donors=13,
                      categories=summary.category.value_counts().to_dict(), **layout)
    (DATA / (PREFIX+'_validation.json')).write_text(json.dumps(validation, indent=2)+'\n')
    print('Verified 40 aligned regions, original ABC counts/colors, fonts, text bounds and exact figure dimensions.')


if __name__ == '__main__':
    main()
