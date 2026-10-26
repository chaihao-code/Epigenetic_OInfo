#!/usr/bin/env python3
"""
Plot heatmap of joint MI / sum of pairwise MI ratios across blocks and anchors.

X-axis: double-layer labels (region groups + vertical feature labels)
Font: Arial, 5 pt
Figure width: 14 cm
PDF fonttype: 42 (embedded)
"""

import os
import glob
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable

# ============================================================================
# Font setup
# ============================================================================
# Try to find Arial, fallback to DejaVu Sans if not available
arial_found = False
for font_path in font_manager.findSystemFonts(fontpaths=None, fontext='ttf'):
    try:
        prop = font_manager.FontProperties(fname=font_path)
        if prop.get_name() == 'Arial' or 'arial' in os.path.basename(font_path).lower():
            font_manager.fontManager.addfont(font_path)
            arial_found = True
            break
    except Exception:
        continue

plt.rcParams['font.family'] = 'Arial' if arial_found else 'DejaVu Sans'
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# ============================================================================
# Configuration
# ============================================================================
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_DIR = os.path.join(BASE_DIR, "plots")
os.makedirs(OUT_DIR, exist_ok=True)

N_BINS_LIST = [3, 5, 10]
ANCHORS = ['ATAC', 'H2A.Z', 'CpG_methyl', 'H3K27ac', 'H3K27me3', 'H3K4me1', 'H3K4me3', 'H3K36me3', 'H3K9me3']

# Layout order for x-axis labels (must match actual block name prefixes)
REGION_ORDER = [
    'promoter', 'proximal_ctcf_bound_enhancer_purified_bins',
    'proximal_not_ctcf_bound_enhancer_purified_bins',
    'distal_ctcf_bound_enhancer_purified_bins',
    'distal_not_ctcf_bound_enhancer_purified_bins',
    'exon_purified_bins', 'intron_purified_bins', 'None'
]
FEATURE_ORDER = ['CpG_island', 'LINE', 'SINE', 'LTR', 'None']

REGION_LABELS = ['promoter', 'proximal\nCTCF+', 'proximal\nCTCF-', 'distal\nCTCF+',
                 'distal\nCTCF-', 'exon', 'intron', 'None']
FEATURE_LABELS_DISP = ['CpG island', 'LINE', 'SINE', 'LTR', 'None']

# Font sizes (5 pt for everything)
FS = 5
FS_TITLE = 5
FS_LABEL = 5
FS_TICK = 5
FS_CBAR = 5


EARLY_EMBRYO_STAGES = ['zygote', '2cell', '4cell', '8cell', 'morula', 'ICM', 'TE']

def get_stages(n_bins):
    mi_dir = os.path.join(BASE_DIR, f"pairwise_mi_n{n_bins}")
    stages = [d for d in os.listdir(mi_dir) if os.path.isdir(os.path.join(mi_dir, d)) and d in EARLY_EMBRYO_STAGES]
    return sorted(stages)


def compute_ratio_matrix(n_bins):
    stages = get_stages(n_bins)
    mi_dir = os.path.join(BASE_DIR, f"pairwise_mi_n{n_bins}")

    # Build ordered block list: region first, then feature
    ORDERED_BLOCKS = [f'{r}_AND_{f}' for r in REGION_ORDER for f in FEATURE_ORDER]
    # Exclude None_AND_None if present
    ORDERED_BLOCKS = [b for b in ORDERED_BLOCKS if b != 'None_AND_None']

    # Collect all blocks present in data
    all_blocks = set()
    for stage in stages:
        joint_file = os.path.join(mi_dir, stage, f"{stage}_mi_joint.csv")
        if os.path.exists(joint_file):
            df = pd.read_csv(joint_file)
            all_blocks.update(df['block'].unique())

    # Filter to only blocks that exist in our ordered list
    blocks = [b for b in ORDERED_BLOCKS if b in all_blocks]

    # Compute ratios per stage
    ratios = []
    for stage in stages:
        jf = os.path.join(mi_dir, stage, f"{stage}_mi_joint.csv")
        if not os.path.exists(jf):
            continue
        block_names = pd.read_csv(jf)['block'].unique()

        for block_name in block_names:
            if block_name not in blocks:
                continue

            pf = os.path.join(mi_dir, stage, f"{stage}_mi_pairwise.csv")
            if not os.path.exists(pf):
                continue

            df_joint = pd.read_csv(jf)
            df_pair = pd.read_csv(pf)
            
            joint_block = df_joint[df_joint['block'] == block_name]
            pair_block = df_pair[df_pair['block'] == block_name]
            
            for _, row in joint_block.iterrows():
                feature = row['feature']
                anchor = feature.split('.', 1)[-1]
                H_X = row['H_X']
                I_X_joint = row['I_X_joint']
                
                if H_X == 0:
                    continue
                
                joint_ratio = I_X_joint / H_X
                
                pair_rows = pair_block[
                    (pair_block['feature_1'] == feature) | (pair_block['feature_2'] == feature)
                ]
                pairwise_sum = pair_rows['I_XY'].sum()
                
                if pairwise_sum == 0:
                    continue
                
                ratio = joint_ratio / pairwise_sum
                ratios.append({
                    'stage': stage,
                    'block': block_name,
                    'anchor': anchor,
                    'ratio': ratio
                })
    
    df_ratios = pd.DataFrame(ratios)
    df_avg = df_ratios.groupby(['block', 'anchor'])['ratio'].mean().reset_index()
    
    # Build matrix: rows=anchors, cols=blocks
    n_rows = len(ANCHORS)
    n_cols = len(blocks)
    matrix = np.full((n_rows, n_cols), np.nan)
    for i, anchor in enumerate(ANCHORS):
        for j, block in enumerate(blocks):
            val = df_avg[(df_avg['anchor'] == anchor) & (df_avg['block'] == block)]['ratio']
            if len(val) > 0:
                matrix[i, j] = val.values[0]
    
    return matrix, blocks


def plot_heatmap(matrix, blocks, n_bins):
    n_rows, n_cols = matrix.shape
    
    # Figure size: width = 14 cm = 5.51 inches
    # Height: proportional to rows, with extra space for labels
    fig_w = 14 / 2.54  # 5.51 inches
    fig_h = 5.5 / 2.54  # 5.5 cm fixed height
    
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.subplots_adjust(bottom=0.42, left=0.14, right=0.82, top=0.92)
    
    # Replace AND with / for display
    x_labels = [b.replace('_AND_', '/') for b in blocks]
    
    # Color: RdBu_r centered at log2(1)=0
    log_mat = np.log2(matrix)
    vmin_log, vmax_log = -3.5, 3.5
    
    im = ax.imshow(log_mat, aspect='auto', cmap='RdBu_r',
                   vmin=vmin_log, vmax=vmax_log,
                   interpolation='nearest',
                   extent=(-0.5, n_cols - 0.5, n_rows - 0.5, -0.5))
    
    # Colorbar (linear log2 scale)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="4%", pad=0.08)
    cbar = fig.colorbar(im, cax=cax)
    cbar_ticks = [-3, -2, -1, 0, 1, 2, 3]
    cbar.set_ticks(cbar_ticks)
    cbar.set_ticklabels([str(t) for t in cbar_ticks])
    cbar.ax.tick_params(labelsize=FS_CBAR)
    cbar.set_label(r'$\log_2 [ I(Y; X_{\rm rest}) / \sum_i I(Y; X_i) ]$',
                   fontsize=FS_CBAR, rotation=270, labelpad=14)
    
    # Y-axis ticks
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(ANCHORS, fontsize=FS_TICK)
    ax.set_ylabel('Feature (anchor)', fontsize=FS_LABEL)
    
    # Hide default x-axis
    ax.set_xticks([])
    ax.set_xticklabels([])
    
    # Group separators: between different regions
    for i in range(1, n_cols):
        prev_region = blocks[i - 1].split('_AND_')[0]
        curr_region = blocks[i].split('_AND_')[0]
        if prev_region != curr_region:
            ax.axvline(i - 0.5, color='black', linewidth=0.8, alpha=0.3)
    
    # Region group labels with double arrows (upper layer)
    for r_idx, r in enumerate(REGION_ORDER):
        indices = [i for i, b in enumerate(blocks) if b.startswith(f'{r}_AND_')]
        if not indices:
            continue
        left = indices[0] - 0.4
        right = indices[-1] + 0.4
        ax.annotate('', xy=(right, -0.08), xytext=(left, -0.08),
                    xycoords=ax.get_xaxis_transform(),
                    textcoords=ax.get_xaxis_transform(),
                    arrowprops=dict(arrowstyle='<->', color='gray', lw=1.5))
        ax.text((left + right) / 2, -0.12, REGION_LABELS[r_idx],
                ha='center', va='top', fontsize=FS_TICK, color='black',
                transform=ax.get_xaxis_transform())
    
    # Feature labels: vertical text below each column (lower layer)
    feature_y = -0.22
    for i, block in enumerate(blocks):
        f = block.split('_AND_', 1)[1]
        label = f.replace('_', ' ')
        ax.text(i, feature_y, label, ha='center', va='top', fontsize=FS_TICK,
                color='black', transform=ax.get_xaxis_transform(), rotation=90)
    
    stage_label = 'early_embryo' if len(EARLY_EMBRYO_STAGES) < 50 else 'all_stages'
    ax.set_title(f'Joint MI / Sum of Pairwise MI  (n_bins={n_bins}, averaged across {stage_label})',
                 fontsize=FS_TITLE, pad=10)
    ax.set_xlim(-0.5, n_cols - 0.5)
    
    out_base = os.path.join(OUT_DIR, f"joint_vs_pairwise_ratio_heatmap_n{n_bins}_{stage_label}_v3")
    plt.savefig(f"{out_base}.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{out_base}.pdf", format='pdf', bbox_inches='tight')
    plt.close()
    print(f"Saved: {out_base}.png/.pdf")


if __name__ == '__main__':
    for n_bins in N_BINS_LIST:
        print(f"Processing n_bins={n_bins}...")
        matrix, blocks = compute_ratio_matrix(n_bins)
        plot_heatmap(matrix, blocks, n_bins)
