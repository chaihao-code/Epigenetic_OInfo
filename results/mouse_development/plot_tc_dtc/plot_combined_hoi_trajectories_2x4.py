#!/usr/bin/env python3
"""
分别绘制 TC 和 DTC trajectory 拼图（3 行 × 2 列）
去掉 exon 行，图例统一放在右侧。
"""

import argparse
import csv
import re
from pathlib import Path

import networkx as nx

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches

# ---------- 全局字体 ----------
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['font.size'] = 5
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['axes.titlesize'] = 5
plt.rcParams['axes.labelsize'] = 5
plt.rcParams['xtick.labelsize'] = 4.5
plt.rcParams['ytick.labelsize'] = 4.5
plt.rcParams['legend.fontsize'] = 5

# ---------- 路径 ----------
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
RESULTS_DIR = PROJECT_DIR / 'results_oinfo_n_bins_5'
TREE_PATH = PROJECT_DIR.parent / 'mouse_development_fine_te/script/plot_O_info/mouse_differentiation_tree.txt'

BLOCKS = [
    ('promoter', 'CpG island', 'promoter_AND_CpG_island'),
    ('promoter', 'LINE', 'promoter_AND_LINE'),
    ('proximal enhancer\nCTCF not bound', 'CpG island', 'proximal_not_ctcf_bound_enhancer_purified_bins_AND_CpG_island'),
    ('proximal enhancer\nCTCF not bound', 'LINE', 'proximal_not_ctcf_bound_enhancer_purified_bins_AND_LINE'),
    ('distal enhancer\nCTCF not bound', 'CpG island', 'distal_not_ctcf_bound_enhancer_purified_bins_AND_CpG_island'),
    ('distal enhancer\nCTCF not bound', 'LINE', 'distal_not_ctcf_bound_enhancer_purified_bins_AND_LINE'),
]

preimpl_color = '#2E86AB'

# 全局颜色映射
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
    if cell_type in ['zygote', '2cell', '4cell', '8cell', 'morula']:
        return cell_type
    elif cell_type in ['ICM', 'TE']:
        return 'E3.5'
    elif 'E6.5' in cell_type:
        return 'E6.5'
    elif 'E7.0' in cell_type:
        return 'E7.0'
    elif 'E7.5' in cell_type:
        return 'E7.5'
    else:
        m = re.search(r'(E\d+\.\d+|P0)', cell_type)
        return m.group(1) if m else cell_type

def get_embryonic_day(cell_type):
    preimpl_days = {'zygote': 0, '2cell': 1, '4cell': 2, '8cell': 2.5,
                    'morula': 3, 'ICM': 3.5, 'TE': 3.5}
    if cell_type in preimpl_days:
        return preimpl_days[cell_type]
    if cell_type.startswith('P0_'):
        return 19.5
    return float(cell_type.split('_')[0][1:])


def load_inputs(results_dir, tree_path, output_dir):
    """Read only overall TC/DTC values and record every missing input."""
    if not results_dir.is_dir():
        raise FileNotFoundError(f'Results directory not found: {results_dir}')
    tree = nx.DiGraph()
    with tree_path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parent, children = line.split('\t')
            for child in children.split(','):
                tree.add_edge(parent, child.strip())
    if 'zygote' not in tree or not nx.is_arborescence(tree):
        raise ValueError('Expected a directed differentiation tree rooted at zygote')
    if nx.descendants(tree, 'zygote') | {'zygote'} != set(tree):
        raise ValueError('All differentiation nodes must be reachable from zygote')
    values = {}
    issues = []
    for cell_type in sorted(tree):
        for _, _, block_label in BLOCKS:
            path = results_dir / cell_type / block_label / f'{cell_type}_{block_label}_oinfo_overall.csv'
            if not path.is_file():
                reason = 'sample directory missing' if not (results_dir / cell_type).is_dir() else 'overall CSV missing'
                issues.append((cell_type, block_label, 'tc,dtc', reason, str(path)))
                continue
            with path.open(newline='') as handle:
                rows = list(csv.DictReader(handle))
            if len(rows) != 1:
                issues.append((cell_type, block_label, 'tc,dtc', f'expected 1 overall row, found {len(rows)}', str(path)))
                continue
            for metric in ('tc', 'dtc'):
                try:
                    value = float(rows[0][metric])
                    if not np.isfinite(value):
                        raise ValueError('nonfinite value')
                except (KeyError, TypeError, ValueError):
                    issues.append((cell_type, block_label, metric, 'missing or invalid metric', str(path)))
                    continue
                values[cell_type, block_label, metric] = value
    report = output_dir / 'tc_dtc_missing_inputs.csv'
    with report.open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['cell_type', 'block_label', 'metric', 'reason', 'expected_path'])
        writer.writerows(issues)
    missing_samples = sorted({row[0] for row in issues if row[3] == 'sample directory missing'})
    print(f'Loaded {len(values)} metric values; input issues: {len(issues)}; report: {report}')
    if missing_samples:
        print('Missing samples in differentiation tree: ' + ', '.join(missing_samples))
    for metric in ('tc', 'dtc'):
        for _, _, block in BLOCKS:
            if not any(key[1:] == (block, metric) for key in values):
                raise ValueError(f'No valid {metric} data for {block}; see {report}')
    return tree, values


def collect_trajectories(block_label, metric, tree, values, start_cell_type='zygote'):
    leaves = [n for n in tree if tree.out_degree(n) == 0]
    trajectories = {}
    for leaf in leaves:
        path = nx.shortest_path(tree, start_cell_type, leaf)
        tissue = leaf.split('_', 1)[1] if '_' in leaf else leaf
        for cell_type in path:
            value = values.get((cell_type, block_label, metric))
            if value is None:
                continue
            trajectories.setdefault(tissue, []).append({
                'cell_type': cell_type,
                'stage': classify_stage(cell_type),
                'e_day': get_embryonic_day(cell_type),
                'x_label': get_x_label(cell_type),
                'value': value,
            })
    if not trajectories:
        return {}
    all_days = {d['e_day'] for lst in trajectories.values() for d in lst}
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
    for lst in trajectories.values():
        for d in lst:
            d['x_pos'] = e_day_to_x[d['e_day']]
    return trajectories

def plot_on_ax(ax, trajectories, linewidth=0.5, markersize=1.8):
    stage_to_all_days = {}
    for lst in trajectories.values():
        for d in lst:
            stage_to_all_days.setdefault(d['stage'], set()).add(d['e_day'])
    for stage in stage_to_all_days:
        stage_to_all_days[stage] = sorted(stage_to_all_days[stage])

    organ_ts = {tissue for tissue, lst in trajectories.items()
                if any(d['stage'] == 'organogenesis' for d in lst)}

    # organ tissues
    for tissue in sorted(organ_ts):
        data = sorted(trajectories[tissue], key=lambda x: x['x_pos'])
        if len(data) < 2:
            continue
        xs = [d['x_pos'] for d in data]
        ys = [d['value'] for d in data]
        stages = [d['stage'] for d in data]
        cts = [d['cell_type'] for d in data]
        eds = [d['e_day'] for d in data]
        for i in range(len(data) - 1):
            stage = stages[i]
            color = preimpl_color if stage == 'preimplementation' else \
                    gastrulation_color_map.get(cts[i], '#666666') if stage == 'gastrulation' else \
                    organ_color_map.get(tissue, '#666666')
            next_stage = stages[i+1]
            is_stage_transition = (stage != next_stage)
            if is_stage_transition:
                nxt_days = stage_to_all_days.get(next_stage, [])
                is_jump = (eds[i+1] != min(nxt_days)) if nxt_days else False
            else:
                days_in_stage = stage_to_all_days.get(stage, [])
                skipped = [dd for dd in days_in_stage if eds[i] < dd < eds[i+1]]
                is_jump = len(skipped) > 0
            ls = '--' if is_jump else '-'
            ax.plot([xs[i], xs[i+1]], [ys[i], ys[i+1]],
                    color=color, linewidth=linewidth, alpha=0.6, linestyle=ls)
        for i in range(len(data)):
            stage = stages[i]
            color = preimpl_color if stage == 'preimplementation' else \
                    gastrulation_color_map.get(cts[i], '#666666') if stage == 'gastrulation' else \
                    organ_color_map.get(tissue, '#666666')
            ax.plot(xs[i], ys[i], 'o', color=color, markersize=markersize,
                    markeredgecolor='black', markeredgewidth=0.3, alpha=0.9)

    # gastrulation-only
    for tissue, data in trajectories.items():
        if tissue in organ_ts:
            continue
        data = sorted(data, key=lambda x: x['x_pos'])
        if len(data) < 2:
            continue
        if not any(d['stage'] == 'gastrulation' for d in data):
            continue
        if any(d['stage'] == 'organogenesis' for d in data):
            continue
        xs = [d['x_pos'] for d in data]
        ys = [d['value'] for d in data]
        stages = [d['stage'] for d in data]
        cts = [d['cell_type'] for d in data]
        eds = [d['e_day'] for d in data]
        for i in range(len(data) - 1):
            color = preimpl_color if stages[i] == 'preimplementation' else \
                    gastrulation_color_map.get(cts[i], '#666666')
            is_jump = (eds[i+1] - eds[i]) > 3
            ls = '--' if is_jump else '-'
            ax.plot([xs[i], xs[i+1]], [ys[i], ys[i+1]],
                    color=color, linewidth=linewidth, alpha=0.6, linestyle=ls)
        for i in range(len(data)):
            color = preimpl_color if stages[i] == 'preimplementation' else \
                    gastrulation_color_map.get(cts[i], '#666666')
            ax.plot(xs[i], ys[i], 'o', color=color, markersize=markersize,
                    markeredgecolor='black', markeredgewidth=0.3, alpha=0.9)

    # 标注 ICM 和 TE
    icm_pts = []
    te_pts = []
    for lst in trajectories.values():
        for d in lst:
            if d['cell_type'] == 'ICM':
                icm_pts.append((d['x_pos'], d['value']))
            elif d['cell_type'] == 'TE':
                te_pts.append((d['x_pos'], d['value']))
    if icm_pts:
        x_icm = np.mean([p[0] for p in icm_pts])
        y_icm = np.mean([p[1] for p in icm_pts])
        ax.text(x_icm - 0.25, y_icm, 'ICM', fontsize=4.5, ha='right', va='center',
                color='black', fontweight='bold')
    if te_pts:
        x_te = np.mean([p[0] for p in te_pts])
        y_te = np.mean([p[1] for p in te_pts])
        ax.text(x_te + 0.25, y_te, 'TE', fontsize=4.5, ha='left', va='center',
                color='black', fontweight='bold')

    # x 轴
    x_info = {}
    for lst in trajectories.values():
        for d in lst:
            x_info[d['x_pos']] = d['x_label']
    sorted_x = sorted(x_info.keys())
    ax.set_xticks(sorted_x)
    ax.set_xticklabels([x_info[x] for x in sorted_x], rotation=45, ha='right', fontsize=4)

    ax.axhline(y=0, color='gray', linestyle='-', alpha=0.3, linewidth=0.4)

    preimpl_max = max((d['x_pos'] for lst in trajectories.values() for d in lst if d['stage'] == 'preimplementation'), default=None)
    gastr_min = min((d['x_pos'] for lst in trajectories.values() for d in lst if d['stage'] == 'gastrulation'), default=None)
    gastr_max = max((d['x_pos'] for lst in trajectories.values() for d in lst if d['stage'] == 'gastrulation'), default=None)
    organ_min = min((d['x_pos'] for lst in trajectories.values() for d in lst if d['stage'] == 'organogenesis'), default=None)
    if preimpl_max is not None and gastr_min is not None:
        ax.axvline(x=(preimpl_max + gastr_min) / 2, color='gray', linestyle='--', alpha=0.5, linewidth=0.6)
    if gastr_max is not None and organ_min is not None:
        ax.axvline(x=(gastr_max + organ_min) / 2, color='gray', linestyle='--', alpha=0.5, linewidth=0.6)

def plot_metric(metric, tree, values, output_dir):
    # ---------- 创建画布 ----------
    fig = plt.figure(figsize=(17/2.54, 10/2.54))
    gs = gridspec.GridSpec(
        3, 3, figure=fig,
        width_ratios=[1, 1, 0.55],
        wspace=0.1, hspace=0.2,
        left=0.07, right=0.98, top=0.93, bottom=0.08
    )

    # ---------- 绘制 6 个子图 ----------
    for idx, (feat1, feat2, block_label) in enumerate(BLOCKS):
        row = idx // 2
        col = idx % 2
        ax = fig.add_subplot(gs[row, col])

        trajectories = collect_trajectories(block_label, metric, tree, values)
        if trajectories:
            plot_on_ax(ax, trajectories)

        ax.set_title(f'{feat1} + {feat2}', fontsize=5, pad=2)
        ax.tick_params(axis='both', which='major', pad=1, length=1.5)

        if col == 0:
            ax.set_ylabel(metric.upper(), fontsize=5, labelpad=1)

        if row == 2:
            ax.set_xlabel('Stage', fontsize=5, labelpad=1)
        else:
            ax.set_xticklabels([])

    # ---------- 右侧统一图例 ----------
    ax_leg = fig.add_subplot(gs[:, 2])
    ax_leg.axis('off')

    legend_elements = []
    legend_elements.append(mpatches.Patch(facecolor='white', edgecolor='none', label='Preimplantation'))
    legend_elements.append(Line2D([0], [0], color=preimpl_color, linewidth=1.0, label='  Preimpl.', alpha=0.6))
    legend_elements.append(mpatches.Patch(facecolor='white', edgecolor='none', label=''))

    legend_elements.append(mpatches.Patch(facecolor='white', edgecolor='none', label='Gastrulation'))
    for node in gastrulation_nodes:
        legend_elements.append(Line2D([0], [0], color=gastrulation_color_map[node], linewidth=1.0, label=f'  {node}', alpha=0.6))
    legend_elements.append(mpatches.Patch(facecolor='white', edgecolor='none', label=''))

    legend_elements.append(mpatches.Patch(facecolor='white', edgecolor='none', label='Organogenesis'))
    for tissue in organ_tissues:
        legend_elements.append(Line2D([0], [0], color=organ_color_map[tissue], linewidth=1.0, label=f'  {tissue}', alpha=0.6))

    ax_leg.legend(
        handles=legend_elements,
        loc='upper left',
        frameon=False,
        borderaxespad=0,
        labelspacing=0.6,
        handletextpad=0.4,
        handlelength=1.2
    )

    fig.suptitle(metric.upper(), fontsize=7, y=0.99)

    # ---------- 保存 ----------
    output_prefix = output_dir / f'combined_{metric}_trajectories_3x2_n5'
    plt.savefig(f'{output_prefix}.pdf', dpi=300, bbox_inches='tight', pad_inches=0.02)
    plt.savefig(f'{output_prefix}.png', dpi=300, bbox_inches='tight', pad_inches=0.02)
    print(f'Saved: {output_prefix}.pdf')
    print(f'Saved: {output_prefix}.png')

    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results-dir', type=Path, default=RESULTS_DIR)
    parser.add_argument('--tree-path', type=Path, default=TREE_PATH)
    parser.add_argument('--output-dir', type=Path, default=SCRIPT_DIR)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    tree, values = load_inputs(args.results_dir, args.tree_path, args.output_dir)
    for metric in ('tc', 'dtc'):
        plot_metric(metric, tree, values, args.output_dir)


if __name__ == '__main__':
    main()
