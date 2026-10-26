#!/usr/bin/env python3
"""Place triplet and core trajectories side by side in an 18 cm wide figure."""
import json
from trajectory_utils import ROOT, MARKS, config, read

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D


def main():
    manifest = config()
    font_path = font_manager.findfont('Arial', fallback_to_default=False)
    plt.rcParams.update({'font.family': 'Arial', 'font.size': 8,
                         'axes.labelsize': 8, 'axes.titlesize': 8,
                         'xtick.labelsize': 8, 'ytick.labelsize': 8,
                         'legend.fontsize': 8, 'pdf.fonttype': 42, 'ps.fonttype': 42})
    stages = ['zygote', '2cell', '4cell', '8cell', 'morula', 'E3.5'] + manifest['late_stages']
    # Categorical spacing keeps all 13 stage labels legible at the requested width.
    positions = {s: i + (.5 if i >= 6 else 0) for i, s in enumerate(stages)}
    fig, axes = plt.subplots(1, 2, figsize=(18 / 2.54, 6.5 / 2.54), sharex=True)
    fig.subplots_adjust(left=.08, right=.98, bottom=.47, top=.91, wspace=.18)
    for ax, kind, title, ylabel in zip(axes,
            ['synergistic_triplet', 'redundancy_core'],
            ['Synergistic triplets', 'Redundancy cores'], ['Triplets', 'Cores']):
        rows = read(ROOT / 'results' / (kind + '_mark_trajectories.tsv'))
        assert len(rows) == 13 * 8
        for mark in MARKS:
            z = [next(r for r in rows if r['stage'] == s and r['mark'] == mark) for s in stages]
            highlight = mark == 'H3K36me3'
            color = '#000000' if highlight else manifest['palette'][mark]
            ax.errorbar([positions[s] for s in stages], [float(r['mean_count']) for r in z],
                        yerr=[float(r['sem_count']) for r in z], color=color, marker='o',
                        markersize=2.8,
                        markeredgecolor='#000000' if highlight else '#606060', markeredgewidth=.3,
                        linewidth=1.1 if highlight else .9,
                        elinewidth=.5, capsize=1.,
                        zorder=10 if highlight else 3)
        ax.axvline(5.75, color='#999999', linestyle='--', linewidth=.5)
        ax.set_ylim(bottom=0)
        ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(width=.5, length=2, pad=2)
        ax.set_title(title, fontsize=8, pad=2)
        ax.set_ylabel(ylabel, fontsize=8, labelpad=1)
        ax.set_xticks([positions[s] for s in stages],
                      ['Zygote', '2-cell', '4-cell', '8-cell', 'Morula', 'E3.5'] + manifest['late_stages'],
                      rotation=90, ha='center')
        ax.set_xlabel('Stage', fontsize=8, labelpad=2)
    axes[0].set_yticks([0, 150, 300])
    axes[1].set_yticks([0, 20, 40])
    axes[1].set_ylim(top=40)
    handles = [Line2D([], [], color='#000000' if m == 'H3K36me3' else manifest['palette'][m],
                      marker='o', markersize=2.8,
                      linewidth=1.1 if m == 'H3K36me3' else .9,
                      markeredgecolor='#606060', markeredgewidth=.3) for m in MARKS]
    legend = fig.legend(handles, [m.replace('_', ' ') for m in MARKS], loc='lower center',
                        bbox_to_anchor=(.5, .02), ncol=4, frameon=False,
                        columnspacing=2.2, handlelength=1.6, handletextpad=.5,
                        labelspacing=.35, borderaxespad=0, borderpad=0)
    fig.canvas.draw()
    # Verify every visible label and the shared legend at the exported size.
    renderer = fig.canvas.get_renderer()
    from matplotlib.text import Text
    for text in fig.findobj(Text):
        if text.get_visible() and text.get_text():
            assert text.get_fontsize() == 8, text.get_text()
            bbox = text.get_window_extent(renderer)
            assert bbox.x0 >= 0 and bbox.y0 >= 0 and bbox.x1 <= fig.bbox.width and bbox.y1 <= fig.bbox.height, text.get_text()
    legend_bbox = legend.get_window_extent(renderer)
    for ax in axes:
        ticks = ax.get_xticklabels()
        for i, first in enumerate(ticks):
            for second in ticks[i + 1:]:
                assert not first.get_window_extent(renderer).overlaps(second.get_window_extent(renderer))
        assert legend_bbox.y1 < ax.xaxis.label.get_window_extent(renderer).y0
    for ext in ['pdf', 'png']:
        fig.savefig(ROOT / 'plots' / ('combined_triplet_core_trajectories.' + ext), dpi=300)
    plt.close(fig)
    (ROOT / 'results/combined_trajectories_validation.json').write_text(json.dumps(dict(
        status='complete', figure_size_cm=[18, 6.5], font_family='Arial', font_file=font_path,
        all_font_sizes_pt=8, shared_legend=True, legend_columns=4, legend_rows=2,
        orientation='two horizontally arranged panels', stage_spacing='categorical with half-step inter-period gap',
        data_changed=False, H3K36me3_color='#000000', H3K36me3_linewidth_pt=1.1,
        H3K36me3_markersize_pt=2.8, text_bounds_checked=True,
        tick_overlap_checked=True, legend_axis_label_separation_checked=True), indent=2) + '\n')
    print('Saved 18 × 6.5 cm combined figure (PDF/PNG), two-row, four-column shared legend and all text 8 pt.')


if __name__ == '__main__':
    main()
