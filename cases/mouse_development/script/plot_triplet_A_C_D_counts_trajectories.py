#!/usr/bin/env python3
"""
绘制 Triplet A / B / C 的 synergistic triplet relative count trajectory 拼图 (E11.5+ with H3K36me3)
（纵轴为 count / theoretical maximum）
格式参考：combined_hoi_trajectories_2x3_n5.pdf

数据来源：results_f1_all_stages/triplet_{A,B,C}_all_stages.tsv
"""

# Portable repository paths; no analysis is executed by this helper.
from pathlib import Path as _ReleasePath
_RELEASE_ROOT = next(p for p in _ReleasePath(__file__).resolve().parents
                     if (p / ".epigenetic_oinfo_root").is_file())
def _release_path(value):
    return value.replace("@ROOT@", str(_RELEASE_ROOT))
def _release_config(value):
    if isinstance(value, dict):
        return {k: _release_config(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_release_config(v) for v in value]
    return _release_path(value) if isinstance(value, str) else value

import os
import sys
import re
import csv

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec

# ---------- 全局字体 ----------
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['font.size'] = 5
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['axes.titlesize'] = 5
plt.rcParams['axes.labelsize'] = 5
plt.rcParams['xtick.labelsize'] = 4.5
plt.rcParams['ytick.labelsize'] = 4.5
plt.rcParams['legend.fontsize'] = 5

# ---------- 路径 ----------
F1_DIR = _release_path('@ROOT@/cases/mouse_after_E11.5/2026-07-14_synergistic_triplet_counts_H3K36me3/results_f1_all_stages')
OUT_DIR = _release_path('@ROOT@/cases/mouse_after_E11.5/2026-07-14_synergistic_triplet_counts_H3K36me3')
OUT_PREFIX = os.path.join(OUT_DIR, 'triplet_A_B_C_counts_trajectories_n5')
os.makedirs(OUT_DIR, exist_ok=True)

REGIONS = [
    ('distal enh. CTCF+ / CpG island', 'distal_ctcf_bound_enhancer_purified_bins_AND_CpG_island'),
    ('proximal enh. CTCF- / LINE', 'proximal_not_ctcf_bound_enhancer_purified_bins_AND_LINE'),
    ('promoter / LTR', 'promoter_AND_LTR'),
    ('exon / LINE', 'exon_purified_bins_AND_LINE'),
]

TRIPLET_CLASSES = ['A', 'B', 'C']
TRIPLET_ROW_COLORS = {
    'A': '#ca2a28',
    'B': '#349939',
    'C': '#ef7c21',
}
# Theoretical maximum number of synergistic order-3 combinations per triplet class
# (10 marks total: 1 DNA methylation, 3 repressive, 6 active)
TRIPLET_THEORETICAL_MAX = {
    'A': 1 * 3 * 6 + 1 * 3,   # DNA + (1 repressive + 1 active) + (2 repressive)
    'B': 3 * 15 + 3 * 6 + 1,  # (1 repressive + 2 active) + (2 repressive + 1 active) + (3 repressive)
    'C': 20 + 15,             # no repressive: 3 active (20) + DNA + 2 active (15)
}

# ---------- 加载数据 ----------
class_count_dfs = {}
for cls in TRIPLET_CLASSES:
    tsv_path = os.path.join(F1_DIR, f'triplet_{cls}_all_stages.tsv')
    df = pd.read_csv(tsv_path, sep='\t')
    class_count_dfs[cls] = df

preimpl_color = '#2E86AB'

gastrulation_nodes = ['E6.5_Epi', 'E6.5_ExE', 'E7.0_A', 'E7.0_ExE', 'E7.0_P',
                      'E7.5_A', 'E7.5_AM', 'E7.5_ExE', 'E7.5_P']
organ_tissues = ['Embryonic_Facial_Prominence', 'Forebrain', 'Heart', 'Hindbrain',
                 'Intestine', 'Kidney', 'Limb', 'Liver', 'Lung', 'Midbrain',
                 'Neural_Tube', 'Stomach']

gastrulation_colors = plt.cm.hsv(np.linspace(0, 1, len(gastrulation_nodes) + 1))[:-1]
organ_colors = plt.cm.hsv(np.linspace(0, 1, len(organ_tissues) + 2))[1:-1]
gastrulation_color_map = {node: gastrulation_colors[i] for i, node in enumerate(gastrulation_nodes)}
organ_color_map = {tissue: organ_colors[i] for i, tissue in enumerate(organ_tissues)}


# ---------- 工具函数 ----------
def classify_stage(cell_type):
    if any(x in cell_type for x in ['zygote', '2cell', '4cell', '8cell', 'morula', 'ICM', 'TE']):
        return 'preimplementation'
    elif cell_type.startswith('E6.5') or cell_type.startswith('E7.0') or cell_type.startswith('E7.5'):
        return 'gastrulation'
    else:
        return 'organogenesis'


def get_x_label(cell_type):
    m = re.search(r'(E\d+\.\d+|P0)', cell_type)
    return m.group(1) if m else cell_type


def get_embryonic_day(cell_type):
    """Map stage name to approximate embryonic day for x-axis positioning."""
    if cell_type.startswith('E'):
        try:
            return float(cell_type.split('_')[0][1:])
        except (ValueError, IndexError):
            return 0.0
    if cell_type.startswith('P0_'):
        return 19.5
    return 0.0


def get_count(cls, cell_type, region):
    """从 triplet_{cls}_all_stages.tsv 中读取 count。"""
    df = class_count_dfs[cls]
    row = df[(df['stage'] == cell_type) & (df['region'] == region)]
    if row.empty:
        return None
    return int(row['count'].iloc[0])


STAGE_GROUP_COLORS = {
    'preimplementation': '#2E86AB',
    'gastrulation': '#A23B72',
    'organogenesis': '#F18F01',
}


def collect_stage_averages(block_label, triplet_cls):
    """按 stage 分组，计算每个 stage 的平均 relative count (count / theoretical max)。"""
    df = class_count_dfs[triplet_cls]
    df_region = df[df['region'] == block_label]
    if df_region.empty:
        return []

    theoretical_max = TRIPLET_THEORETICAL_MAX[triplet_cls]
    records = []
    for _, row in df_region.iterrows():
        cell_type = row['stage']
        count = int(row['count'])
        relative = count / theoretical_max
        stage = classify_stage(cell_type)
        e_day = get_embryonic_day(cell_type)
        x_label = get_x_label(cell_type)
        records.append({
            'cell_type': cell_type,
            'stage': stage,
            'e_day': e_day,
            'x_label': x_label,
            'relative': relative,
        })
    if not records:
        return []

    # 按 x_label 分组并计算平均
    from collections import defaultdict
    groups = defaultdict(list)
    for r in records:
        groups[r['x_label']].append(r)

    all_days = {r['e_day'] for r in records}
    preimpl_days = sorted([d for d in all_days if d <= 3.5])
    gastr_days = sorted([d for d in all_days if 6 <= d <= 8])
    organ_days = sorted([d for d in all_days if d > 8])
    e_day_to_x = {}
    for i, d in enumerate(preimpl_days):
        e_day_to_x[d] = i
    for i, d in enumerate(gastr_days):
        e_day_to_x[d] = len(preimpl_days) + i
    if organ_days:
        base_x = len(preimpl_days) + len(gastr_days)
        min_od = min(organ_days)
        for d in organ_days:
            e_day_to_x[d] = base_x + (d - min_od) * 2

    stage_avgs = []
    for x_label, grp in groups.items():
        relatives = [r['relative'] for r in grp]
        mean_e_day = np.mean([r['e_day'] for r in grp])
        stage = grp[0]['stage']
        stage_avgs.append({
            'x_label': x_label,
            'stage': stage,
            'mean_e_day': mean_e_day,
            'x_pos': e_day_to_x[mean_e_day],
            'mean_relative': np.mean(relatives),
            'std_relative': np.std(relatives, ddof=0),
            'sem_relative': np.std(relatives, ddof=0) / np.sqrt(len(relatives)) if len(relatives) > 1 else 0,
            'n': len(relatives),
        })

    stage_avgs.sort(key=lambda x: x['x_pos'])
    return stage_avgs


def plot_stage_averages(ax, stage_avgs, color='#333333', linewidth=0.8, markersize=3.5, row_ymax=None):
    if not stage_avgs:
        return

    xs = [d['x_pos'] for d in stage_avgs]
    ys = [d['mean_relative'] for d in stage_avgs]
    yerrs = [d['sem_relative'] for d in stage_avgs]

    # 误差线
    ax.errorbar(xs, ys, yerr=yerrs, fmt='none', ecolor=color, capsize=1.5,
                capthick=0.5, elinewidth=0.5, alpha=0.4)

    # 连接线
    ax.plot(xs, ys, color=color, linewidth=linewidth, alpha=0.7, zorder=1)

    # 点
    for x, y in zip(xs, ys):
        ax.plot(x, y, 'o', color=color,
                markersize=markersize, markeredgecolor='black', markeredgewidth=0.3,
                alpha=0.9, zorder=2)

    # x 轴
    x_info = {d['x_pos']: d['x_label'] for d in stage_avgs}
    sorted_x = sorted(x_info.keys())
    ax.set_xticks(sorted_x)
    ax.set_xticklabels([x_info[x] for x in sorted_x], rotation=45, ha='right', fontsize=4)

    ax.axhline(y=0, color='gray', linestyle='-', alpha=0.3, linewidth=0.4)
    if row_ymax is not None:
        ax.set_ylim(-row_ymax * 0.05, row_ymax * 1.05)
    else:
        ax.set_ylim(bottom=0)

    preimpl_max = max((d['x_pos'] for d in stage_avgs if d['stage'] == 'preimplementation'), default=None)
    gastr_min = min((d['x_pos'] for d in stage_avgs if d['stage'] == 'gastrulation'), default=None)
    gastr_max = max((d['x_pos'] for d in stage_avgs if d['stage'] == 'gastrulation'), default=None)
    organ_min = min((d['x_pos'] for d in stage_avgs if d['stage'] == 'organogenesis'), default=None)
    if preimpl_max is not None and gastr_min is not None:
        ax.axvline(x=(preimpl_max + gastr_min) / 2, color='gray', linestyle='--', alpha=0.5, linewidth=0.6)
    if gastr_max is not None and organ_min is not None:
        ax.axvline(x=(gastr_max + organ_min) / 2, color='gray', linestyle='--', alpha=0.5, linewidth=0.6)


# ---------- 创建画布 ----------
fig = plt.figure(figsize=(10.5/2.54, 10.5/2.54))
gs = gridspec.GridSpec(
    3, 4, figure=fig,
    width_ratios=[1, 1, 1, 1],
    wspace=0.08, hspace=0.22,
    left=0.06, right=0.98, top=0.91, bottom=0.14
)

# ---------- 统一 y 轴范围（相对值 0-1） ----------
row_ymax = {cls: 1.0 for cls in TRIPLET_CLASSES}

# ---------- 绘制 3×3 子图 ----------
all_summary_rows = []

for row_idx, triplet_cls in enumerate(TRIPLET_CLASSES):
    row_color = TRIPLET_ROW_COLORS[triplet_cls]
    for col_idx, (title, block_label) in enumerate(REGIONS):
        ax = fig.add_subplot(gs[row_idx, col_idx])

        stage_avgs = collect_stage_averages(block_label, triplet_cls)
        if stage_avgs:
            plot_stage_averages(ax, stage_avgs, color=row_color, row_ymax=row_ymax[triplet_cls])
            for d in stage_avgs:
                all_summary_rows.append({
                    'triplet_class': triplet_cls,
                    'block': block_label,
                    'stage': d['stage'],
                    'x_label': d['x_label'],
                    'x_pos': d['x_pos'],
                    'mean_relative': d['mean_relative'],
                    'std_relative': d['std_relative'],
                    'sem_relative': d['sem_relative'],
                    'n': d['n'],
                })

        if row_idx == 0:
            ax.set_title(title, fontsize=5, pad=2)
        else:
            ax.set_title('')

        ax.tick_params(axis='both', which='major', pad=1, length=1.5)

        if col_idx == 0:
            ax.set_ylabel(f'Triplet {triplet_cls}\noccurrence proportion', fontsize=5, labelpad=1, color=row_color)
            ax.tick_params(axis='y', colors=row_color)
        else:
            ax.set_yticklabels([])

        if row_idx == 2:
            ax.set_xlabel('Stage', fontsize=5, labelpad=1)
        else:
            ax.set_xticklabels([])

# ---------- 保存图像 ----------
plt.savefig(f'{OUT_PREFIX}.pdf', dpi=300, bbox_inches='tight', pad_inches=0.02)
plt.savefig(f'{OUT_PREFIX}.png', dpi=300, bbox_inches='tight', pad_inches=0.02)
print(f'Saved: {OUT_PREFIX}.pdf')
print(f'Saved: {OUT_PREFIX}.png')

# ---------- 保存汇总 TSV ----------
if all_summary_rows:
    summary_tsv = os.path.join(OUT_DIR, 'triplet_A_B_C_relative_counts_summary.tsv')
    with open(summary_tsv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['triplet_class', 'block', 'stage', 'x_label', 'x_pos', 'mean_relative', 'std_relative', 'sem_relative', 'n'], delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(all_summary_rows)
    print(f'Saved summary TSV: {summary_tsv}')

plt.close(fig)
