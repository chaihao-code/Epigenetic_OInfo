#!/usr/bin/env python3
"""
Plot normalized MI trends for EACH BLOCK and EACH ANCHOR separately.

One figure per (block, anchor) pair, showing all 6 developmental stages.
Each figure contains curves for the anchor vs every other signal + joint total.

Usage:
    python plot_by_block_all_anchors.py --n-bins 5 --block distal_ctcf_bound_enhancer_AND_CpG_island --anchor H3K4me3
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

import argparse
from functools import lru_cache
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.font_manager import FontProperties

# ============================================================================
# Font setup
# ============================================================================
from matplotlib import rcParams
from matplotlib.font_manager import FontProperties, fontManager
import matplotlib.pyplot as plt

font_path = _release_path('@ROOT@/data/fonts/arial.ttf')

# 检查文件是否存在
import os
if os.path.exists(font_path):
    # 注册字体
    fontManager.addfont(font_path)
    # 获取字体的名称（可能是 'Arial' 或 'ArialMT'）
    prop = FontProperties(fname=font_path)
    font_name = prop.get_name()
    rcParams['font.family'] = font_name
    print(f"使用字体: {font_name}")
else:
    rcParams['font.family'] = 'DejaVu Sans'
    print("Arial 未找到，使用 DejaVu Sans")

# 其他设置
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42

# ============================================================================
# Configuration
# ============================================================================

BASE_DIR = Path(__file__).resolve().parents[1]
STAGES = ['zygote', '2cell', '4cell', '8cell', 'morula', 'ICM']
STAGE_LABELS = ['zygote', '2cell', '4cell', '8cell', 'morula', 'ICM']

FS = 5  # global font size in pt

TARGET_COLORS = {
    'ATAC': '#999999',
    'H3K27ac': '#E41A1C',
    'H3K27me3': '#377EB8',
    'H3K4me1': '#4DAF4A',
    'H3K4me3': '#984EA3',
    'H3K9me3': '#FF7F00',
    'DNA methylation': '#A65628',
    'total': '#000000',
}
TARGET_LINESTYLES = {
    'ATAC': '-',
    'H3K27ac': '-',
    'H3K27me3': '-',
    'H3K4me1': '-',
    'H3K4me3': '-',
    'H3K9me3': '-',
    'DNA methylation': '-',
    'total': '--',
}


def get_result_dir(n_bins: int) -> Path:
    return BASE_DIR / f"pairwise_mi_n{n_bins}"


@lru_cache(maxsize=None)
def get_feature_list(n_bins: int) -> list:
    """Get list of base feature names from the first stage."""
    res_dir = get_result_dir(n_bins)
    stage_dir = res_dir / STAGES[0]
    self_path = stage_dir / f"{STAGES[0]}_mi_self.csv"
    df = pd.read_csv(self_path)
    bases = set()
    for name in df['feature'].unique():
        if '.' in name:
            bases.add(name.split('.', 1)[1])
        else:
            bases.add(name)
    return sorted(bases)


def get_block_list(n_bins: int) -> list:
    res_dir = get_result_dir(n_bins)
    stage_dir = res_dir / STAGES[0]
    pairwise_path = stage_dir / f"{STAGES[0]}_mi_pairwise.csv"
    df = pd.read_csv(pairwise_path)
    return sorted(df['block'].unique())


@lru_cache(maxsize=None)
def load_stage_block_data(stage: str, block: str, n_bins: int):
    res_dir = get_result_dir(n_bins)
    stage_dir = res_dir / stage
    df_p = pd.read_csv(stage_dir / f"{stage}_mi_pairwise.csv")
    df_j = pd.read_csv(stage_dir / f"{stage}_mi_joint.csv")
    sub_p = df_p[df_p['block'] == block]
    sub_j = df_j[df_j['block'] == block]
    return sub_p, sub_j


def extract_pairwise_norm(sub_p: pd.DataFrame, anchor_kw: str, target_kw: str):
    if sub_p.empty:
        return np.nan
    for _, row in sub_p.iterrows():
        f1, f2 = row['feature_1'], row['feature_2']
        if anchor_kw in f1 and target_kw in f2:
            H_a, H_t = row['H_X'], row['H_Y']
            I_xy = row['I_XY']
            break
        elif anchor_kw in f2 and target_kw in f1:
            H_a, H_t = row['H_Y'], row['H_X']
            I_xy = row['I_XY']
            break
    else:
        return np.nan
    denom = H_a
    return I_xy / denom if denom > 0 else 0.0


def extract_joint_norm(sub_j: pd.DataFrame, anchor_kw: str):
    if sub_j.empty:
        return np.nan
    mask = sub_j['feature'].astype(str).str.contains(anchor_kw)
    sub = sub_j[mask]
    if sub.empty:
        return np.nan
    row = sub.iloc[0]
    H_x = row['H_X']
    H_minus = row['H_X_minus_i']
    I_joint = row['I_X_joint']
    denom = H_x
    return I_joint / denom if denom > 0 else 0.0


def plot_block_anchor(block: str, anchor: str, n_bins: int, output_dir: Path):
    all_features = get_feature_list(n_bins)
    targets = [f for f in all_features if f != anchor] + ['total']
    
    x = np.arange(len(STAGE_LABELS))
    fig_w = 10 / 2.54   # 10 cm -> inches
    fig_h = 7 / 2.54    # 7 cm -> inches
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    csv_rows = []
    
    for target in targets:
        y_vals = []
        for stage, label in zip(STAGES, STAGE_LABELS):
            sub_p, sub_j = load_stage_block_data(stage, block, n_bins)
            if target == 'total':
                val = extract_joint_norm(sub_j, anchor)
            else:
                val = extract_pairwise_norm(sub_p, anchor, target)
            y_vals.append(val)
            csv_rows.append({
                'anchor': anchor,
                'target': target,
                'n_bins': n_bins,
                'stage': label,
                'block': block,
                'normalized_mi': val,
            })
        
        if all(np.isnan(y_vals)):
            continue
        
        ax.plot(
            x, y_vals,
            label=target,
            color=TARGET_COLORS.get(target, None),
            linestyle=TARGET_LINESTYLES.get(target, '-'),
            linewidth=1.0 if target == 'total' else 0.8,
            marker='o',
            markersize=3,
            markerfacecolor='white',
            markeredgewidth=0.8,
            markeredgecolor=TARGET_COLORS.get(target, None),
        )
    
    ax.set_xticks(x)
    ax.set_xticklabels(STAGE_LABELS, rotation=30, ha='right', fontsize=FS)
    ax.set_ylabel(r'$I(X;Y) / H(X)$', fontsize=FS)
    ax.set_title(f'{anchor} | {block} | n={n_bins}', fontsize=FS)
    ax.tick_params(axis='both', which='major', labelsize=FS)
    ax.legend(
        loc='upper center',
        bbox_to_anchor=(0.5, -0.28),
        ncol=4,
        frameon=True,
        fontsize=FS,
    )
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.30)
    
    safe_block = block.replace('/', '_').replace('&', '_AND_')
    png_path = output_dir / f"{safe_block}_{anchor}_mi_trend.png"
    pdf_path = output_dir / f"{safe_block}_{anchor}_mi_trend.pdf"
    plt.savefig(png_path, dpi=300)
    plt.savefig(pdf_path, format='pdf', bbox_inches='tight')
    plt.close()
    
    csv_path = output_dir / f"{safe_block}_{anchor}_mi_trend.csv"
    pd.DataFrame(csv_rows).to_csv(csv_path, index=False)
    
    return png_path.name


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-bins", type=int, required=True, choices=[3, 5, 10])
    parser.add_argument("--block", type=str, default=None,
                        help="Specific block to plot (e.g., distal_ctcf_bound_enhancer_AND_CpG_island). If omitted, plot all.")
    parser.add_argument("--anchor", type=str, default=None,
                        help="Specific anchor feature to plot (e.g., H3K4me3). If omitted, plot all.")
    args = parser.parse_args()
    
    n_bins = args.n_bins
    output_dir = BASE_DIR / "plots" / f"n{n_bins}"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    blocks = get_block_list(n_bins)
    anchors = get_feature_list(n_bins)
    
    if args.block:
        if args.block not in blocks:
            raise ValueError(f"Block '{args.block}' not found. Available: {blocks}")
        blocks = [args.block]
    
    if args.anchor:
        if args.anchor not in anchors:
            raise ValueError(f"Anchor '{args.anchor}' not found. Available: {anchors}")
        anchors = [args.anchor]
    
    total = len(blocks) * len(anchors)
    print(f"n_bins={n_bins}: {len(blocks)} blocks × {len(anchors)} anchors = {total} figures")
    
    cnt = 0
    for block in blocks:
        for anchor in anchors:
            fname = plot_block_anchor(block, anchor, n_bins, output_dir)
            cnt += 1
            if cnt % 50 == 0 or cnt == total:
                print(f"  [{cnt}/{total}] {fname}")
    
    print(f"All done! Output: {output_dir}")


if __name__ == "__main__":
    main()
