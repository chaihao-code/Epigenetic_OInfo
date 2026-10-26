#!/usr/bin/env python3
"""
绘制单个 Triplet 组合出现次数的 bar + dot matrix 图
参考：Occurrence counts of individual Triplet B combinations

统计所有 stage/region 中，synergistic 的 order-3 triplet 各具体组合出现的总次数。
此版本在原始信号基础上额外包含 H3K9ac 和 H3K4me2。
"""

import os
import sys
import csv
import glob
import re
from collections import Counter

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

# ---------- 字体 ----------
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

# ---------- 路径 ----------
RESULTS_DIR = '/media/dell/data4/ch_sequence_entropy/mouse_development_after_E11.5/results_oinfo_n_bins_5'
OUT_DIR = '/media/dell/data5/sequence_feature/sequence_entropy/result/hoi_analysis/projects/mouse_development_E11.5_plus/2026-07-27_synergistic_triplet_counts_all_marks'
os.makedirs(OUT_DIR, exist_ok=True)

# ---------- 配置 ----------
# 核心信号（包含 H3K36me3，以及新增的 H3K9ac、H3K4me2）
KNOWN_MARKS = ['ATAC', 'H3K27ac', 'H3K4me3', 'H3K4me2', 'H3K4me1', 'H3K27me3', 'H3K9ac', 'H3K9me3', 'H3K36me3']
DNA_METH_MARKS = {'CpG_methyl', 'WGBS'}
ALL_MARKS = KNOWN_MARKS + ['DNA_methylation']

# dot matrix 里从上到下的顺序（显示用标签）
MARK_ORDER = ['ATAC', 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K27ac', 'H3K9ac', 'H3K9me3', 'H3K27me3', 'H3K36me3', 'DNA methylation']


def normalize_marks(marks):
    """把 CpG_methyl / WGBS 统一为 DNA_methylation。"""
    normalized = set()
    for m in marks:
        if m in DNA_METH_MARKS:
            normalized.add('DNA_methylation')
        else:
            normalized.add(m)
    return normalized


def extract_marks(remaining_names):
    """从 remaining_names 字符串中提取 mark 名称。"""
    marks = set()
    # 先替换路径分隔符
    text = remaining_names.replace('/', '_')
    for mark in KNOWN_MARKS:
        if mark in text:
            marks.add(mark)
    for mark in DNA_METH_MARKS:
        if mark in text:
            marks.add(mark)
    return normalize_marks(marks)


def classify_triplet(marks):
    """分类 triplet，返回 'A', 'B', 'C' 或 None（H3K36me3 也视为 repressive-like）。"""
    has_dna = 'DNA_methylation' in marks
    has_repressive = any(m in marks for m in ['H3K9me3', 'H3K27me3', 'H3K36me3'])

    if has_dna and has_repressive:
        return 'A'
    if not has_dna and has_repressive:
        return 'B'
    if not has_repressive:
        return 'C'
    return None


def count_combinations(triplet_cls):
    """统计指定 triplet class 的各组合 synergistic 出现次数。"""
    counter = Counter()

    for stage in sorted(os.listdir(RESULTS_DIR)):
        stage_dir = os.path.join(RESULTS_DIR, stage)
        if not os.path.isdir(stage_dir):
            continue
        for region in sorted(os.listdir(stage_dir)):
            region_dir = os.path.join(stage_dir, region)
            if not os.path.isdir(region_dir):
                continue
            pattern = os.path.join(region_dir, '*_oinfo_missing_*.csv')
            for csv_path in glob.glob(pattern):
                with open(csv_path, newline='') as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        if row.get('order', '') != '3':
                            continue
                        if row.get('interpretation', '') != 'Synergy':
                            continue
                        marks = extract_marks(row.get('remaining_names', ''))
                        if len(marks) != 3:
                            continue
                        if classify_triplet(marks) != triplet_cls:
                            continue
                        combo = '/'.join(sorted(marks))
                        counter[combo] += 1

    return counter


def parse_combo(combo_str):
    """把 'ATAC/H3K4me1/H3K9me3' 解析成 mark 集合。"""
    return set(combo_str.split('/'))


def plot_occurrence(triplet_cls, min_count=0, top_n=25):
    counter = count_combinations(triplet_cls)
    # 过滤并排序
    combos = [(combo, count) for combo, count in counter.items() if count >= min_count]
    combos.sort(key=lambda x: (-x[1], x[0]))
    # 只保留 Top N，避免组合太多画不下
    combos = combos[:top_n]

    if not combos:
        print(f'No Triplet {triplet_cls} combinations found.')
        return

    n_combos = len(combos)
    combo_names = [c[0] for c in combos]
    counts = np.array([c[1] for c in combos])

    # 准备 dot matrix 数据
    dot_matrix = np.zeros((len(MARK_ORDER), n_combos))
    for j, combo in enumerate(combo_names):
        marks = parse_combo(combo)
        for i, mark in enumerate(MARK_ORDER):
            if mark in marks:
                dot_matrix[i, j] = 1

    # ---------- 画图 ----------
    # Figure size: 9 cm x 6 cm
    fig = plt.figure(figsize=(9/2.54, 6/2.54))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.15, 1.15], hspace=0.06,
                          left=0.16, right=0.98, top=0.86, bottom=0.18)

    ax_bar = fig.add_subplot(gs[0, 0])
    ax_dot = fig.add_subplot(gs[1, 0], sharex=ax_bar)

    # bar plot
    x = np.arange(n_combos)
    ax_bar.bar(x, counts, color='skyblue', edgecolor='black', linewidth=0.4)
    ax_bar.set_title(f'Triplet {triplet_cls} combination occurrence', fontsize=5, fontproperties=arial_prop, pad=2)
    ax_bar.tick_params(axis='x', bottom=False, labelbottom=False)
    ax_bar.tick_params(axis='y', labelsize=4.5, pad=1)
    ax_bar.set_ylabel('Count', fontsize=5, fontproperties=arial_prop, labelpad=1)

    # dot matrix
    dot_spacing = 1.1
    for i, mark in enumerate(MARK_ORDER):
        y = np.full(n_combos, (len(MARK_ORDER) - 1 - i) * dot_spacing)
        present = dot_matrix[i, :].astype(bool)
        # horizontal gray reference line for each mark row
        ax_dot.axhline(y[0], color='lightgray', linewidth=0.4, zorder=0)
        # dots, half size of previous
        ax_dot.scatter(x[present], y[present], color='black', s=2.5, zorder=2)

    # vertical lines connecting adjacent dots within each combination/column
    for j in range(n_combos):
        present_ys = []
        for i in range(len(MARK_ORDER)):
            y = (len(MARK_ORDER) - 1 - i) * dot_spacing
            if dot_matrix[i, j]:
                present_ys.append(y)
        present_ys.sort()
        for k in range(len(present_ys) - 1):
            ax_dot.plot([j, j], [present_ys[k], present_ys[k + 1]],
                        color='black', linewidth=0.5, zorder=1)

    ax_dot.set_xlim(-0.5, n_combos - 0.5)
    ax_dot.set_ylim(-0.5, (len(MARK_ORDER) - 1) * dot_spacing + 0.5)
    ax_dot.set_yticks(np.arange(len(MARK_ORDER)) * dot_spacing)
    ax_dot.set_yticklabels(list(reversed(MARK_ORDER)), fontsize=4.5, fontproperties=arial_prop)
    ax_dot.set_xticks(x)
    ax_dot.set_xticklabels([])
    ax_dot.tick_params(axis='x', bottom=False)

    # remove dot matrix border
    for spine in ax_dot.spines.values():
        spine.set_visible(False)

    # hide grid
    ax_bar.grid(False)
    ax_dot.grid(False)

    out_path = os.path.join(OUT_DIR, f'triplet_{triplet_cls}_combination_occurrence')
    plt.savefig(f'{out_path}.pdf', dpi=300, bbox_inches='tight', pad_inches=0.02)
    plt.savefig(f'{out_path}.png', dpi=300, bbox_inches='tight', pad_inches=0.02)
    print(f'Saved: {out_path}.pdf')
    print(f'Saved: {out_path}.png')

    # 保存 TSV
    tsv_path = os.path.join(OUT_DIR, f'triplet_{triplet_cls}_combination_occurrence.tsv')
    with open(tsv_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        writer.writerow(['combination', 'count'])
        writer.writerows(combos)
    print(f'Saved: {tsv_path}')

    plt.close(fig)


if __name__ == '__main__':
    cls = sys.argv[1] if len(sys.argv) > 1 else 'D'
    plot_occurrence(cls)
