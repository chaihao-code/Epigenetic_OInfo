#!/usr/bin/env python3
"""
Plot occurrence frequency of top redundancy cores (bar + cumulative line + dot matrix).
Input: all_supreme_cores_filtered.csv
"""

import os
import csv
from collections import Counter

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.gridspec import GridSpec

ARIAL_PATH = '/usr/share/fonts/truetype/msttcorefonts/Arial.ttf'
try:
    arial_prop = FontProperties(fname=ARIAL_PATH, size=5)
    from matplotlib import font_manager
    font_manager.fontManager.addfont(ARIAL_PATH)
    plt.rcParams['font.family'] = 'Arial'
except Exception:
    arial_prop = FontProperties(size=5)

plt.rcParams['font.size'] = 5
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

IN_CSV = '/media/dell/data5/sequence_feature/sequence_entropy/result/hoi_analysis/projects/mouse_preimp/supreme_cores_batch_v4_replacement_bg_mad_z2.5_bic_k1k2_eta05/all_supreme_cores_filtered.csv'
OUT_DIR = '/media/dell/data5/sequence_feature/sequence_entropy/result/hoi_analysis/projects/mouse_preimp/supreme_cores_batch_v4_replacement_bg_mad_z2.5_bic_k1k2_eta05'
OUT_PREFIX = os.path.join(OUT_DIR, 'core_occurrence_frequency_top25')
os.makedirs(OUT_DIR, exist_ok=True)

TOP_N = 25


def normalize_signals(sig_str):
    parts = sig_str.split(' + ')
    parts = ['DNA methylation' if p == 'WGBS' else p for p in parts]
    return sorted(parts)


def main():
    df = pd.read_csv(IN_CSV)
    total_cores = len(df)
    combo_counter = Counter()
    for sig in df['signals']:
        combo_counter['/'.join(normalize_signals(sig))] += 1

    top_combos = combo_counter.most_common(TOP_N)
    n_plot = len(top_combos)

    counts = np.array([c for _, c in top_combos])
    proportions = counts / total_cores * 100
    cumulative = np.cumsum(proportions)

    # Fixed vertical order: active marks by frequency, then repressive/DNA/H3K9me3 at bottom
    MARK_ORDER = ['H3K27ac', 'ATAC', 'H3K4me3', 'H3K4me1', 'H2A.Z',
                  'H3K36me3', 'DNA methylation', 'H3K27me3', 'H3K9me3']

    combo_names = [c for c, _ in top_combos]
    dot_matrix = np.zeros((len(MARK_ORDER), n_plot))
    for j, combo in enumerate(combo_names):
        marks = set(combo.split('/'))
        for i, mark in enumerate(MARK_ORDER):
            if mark in marks:
                dot_matrix[i, j] = 1

    # Figure: fixed 9 cm width
    width_cm = 9.0
    fig = plt.figure(figsize=(width_cm / 2.54, 6.0 / 2.54))
    gs = GridSpec(2, 1, figure=fig, height_ratios=[1.4, 1.0], hspace=0.08,
                  left=0.10, right=0.88, top=0.86, bottom=0.22)

    ax_bar = fig.add_subplot(gs[0, 0])
    ax_cum = ax_bar.twinx()
    ax_dot = fig.add_subplot(gs[1, 0], sharex=ax_bar)

    x = np.arange(n_plot)

    # Bar plot
    bars = ax_bar.bar(x, proportions, color='steelblue', edgecolor='black', linewidth=0.4, alpha=0.7)
    ax_bar.set_ylabel('Proportion (%)', color='steelblue', fontsize=5, fontproperties=arial_prop)
    ax_bar.tick_params(axis='y', labelcolor='steelblue', labelsize=4.5)
    ax_bar.set_ylim(0, max(proportions) * 1.2)
    ax_bar.set_title('Occurrence frequency of the top redundancy cores', fontsize=6, fontproperties=arial_prop, pad=3)
    ax_bar.tick_params(axis='x', bottom=False, labelbottom=False)

    # Cumulative line
    ax_cum.plot(x, cumulative, color='crimson', marker='o', markersize=3, linewidth=0.8)
    ax_cum.set_ylabel('Cumulative (%)', color='crimson', fontsize=5, fontproperties=arial_prop)
    ax_cum.tick_params(axis='y', labelcolor='crimson', labelsize=4.5)
    ax_cum.set_ylim(0, 105)
    ax_cum.axhline(y=100, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)

    # Dot matrix
    dot_spacing = 1.1
    for i, mark in enumerate(MARK_ORDER):
        y = np.full(n_plot, (len(MARK_ORDER) - 1 - i) * dot_spacing)
        present = dot_matrix[i, :].astype(bool)
        ax_dot.axhline(y[0], color='lightgray', linewidth=0.4, zorder=0)
        ax_dot.scatter(x[present], y[present], color='black', s=3.0, zorder=2)

    # vertical lines connecting dots in same combo
    for j in range(n_plot):
        present_ys = []
        for i in range(len(MARK_ORDER)):
            y = (len(MARK_ORDER) - 1 - i) * dot_spacing
            if dot_matrix[i, j]:
                present_ys.append(y)
        present_ys.sort()
        for k in range(len(present_ys) - 1):
            ax_dot.plot([j, j], [present_ys[k], present_ys[k + 1]], color='black', linewidth=0.5, zorder=1)

    ax_dot.set_xlim(-0.5, n_plot - 0.5)
    ax_dot.set_ylim(-0.5, (len(MARK_ORDER) - 1) * dot_spacing + 0.5)
    ax_dot.set_yticks(np.arange(len(MARK_ORDER)) * dot_spacing)
    ax_dot.set_yticklabels(list(reversed(MARK_ORDER)), fontsize=4.5, fontproperties=arial_prop)
    ax_dot.set_xticks(x)
    ax_dot.set_xticklabels([])
    ax_dot.tick_params(axis='x', bottom=False)

    for spine in ax_dot.spines.values():
        spine.set_visible(False)

    ax_bar.grid(False)
    ax_dot.grid(False)

    plt.savefig(f'{OUT_PREFIX}.pdf', dpi=300, bbox_inches='tight', pad_inches=0.02)
    plt.savefig(f'{OUT_PREFIX}.png', dpi=300, bbox_inches='tight', pad_inches=0.02)
    print(f'Saved: {OUT_PREFIX}.pdf')
    print(f'Saved: {OUT_PREFIX}.png')

    # Save TSV
    tsv_path = os.path.join(OUT_DIR, 'core_occurrence_frequency_top25.tsv')
    with open(tsv_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        writer.writerow(['rank', 'combination', 'count', 'proportion_pct', 'cumulative_pct'])
        for rank, ((combo, cnt), prop, cum) in enumerate(zip(top_combos, proportions, cumulative), start=1):
            writer.writerow([rank, combo, cnt, f'{prop:.2f}', f'{cum:.2f}'])
    print(f'Saved: {tsv_path}')

    plt.close(fig)


if __name__ == '__main__':
    main()
