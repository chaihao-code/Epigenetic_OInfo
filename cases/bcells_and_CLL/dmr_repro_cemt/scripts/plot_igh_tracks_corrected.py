#!/usr/bin/env python3
"""
Corrected version of the IGH methylation track plot
(scripts/plot_methylation_diff_chr14_with_gene_counts.py).

Differences from the original:
  - uses the NA-corrected 1kb matrix (dmr_analysis_v2) instead of .bp files:
    zero-coverage bins are NaN -> tracks show gaps instead of dropping to 0
  - group means are nanmeans over covered samples only
  - gold marks |delta_beta| > 0.2 on corrected means (not the zero-filled diff)
  - adds a coverage strip for the three dropout samples (CEMT0006/0095/0097)
"""

from pathlib import Path
import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch

mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams.update({'font.size': 5, 'axes.titlesize': 5, 'axes.labelsize': 5,
                     'xtick.labelsize': 5, 'ytick.labelsize': 5, 'legend.fontsize': 5,
                     'pdf.fonttype': 42, 'ps.fonttype': 42})

ROOT = Path(__file__).resolve().parents[1]
MATRIX = ROOT / 'input/all_bins_1kb_cemt_g2_g3.csv'
GTF = ROOT / 'input/gencode.v47.annotation.gtf'
OUT_PDF = str(ROOT / 'output/chr14_105800000_106500000_methylation_tracks_corrected.pdf')
OUT_PNG = OUT_PDF.replace('.pdf', '.png')

GROUP2 = ['CEMT0004', 'CEMT0005', 'CEMT0025', 'CEMT0026',
          'CEMT0027', 'CEMT0028', 'CEMT0029', 'CEMT0030']
GROUP3 = ['CEMT0006', 'CEMT0094', 'CEMT0095', 'CEMT0096', 'CEMT0097']
ALL_SAMPLES = GROUP2 + GROUP3

# IGHV mutation status (M = mutated, UM = unmutated)
IGHV = {'CEMT0004': 'UM', 'CEMT0005': 'M', 'CEMT0006': 'M',
        'CEMT0025': 'M', 'CEMT0026': 'M', 'CEMT0027': 'UM',
        'CEMT0028': 'M', 'CEMT0029': 'UM', 'CEMT0030': 'UM',
        'CEMT0094': 'M', 'CEMT0095': 'UM', 'CEMT0096': 'M',
        'CEMT0097': 'M'}

CHROM, START, END, BIN = 'chr14', 105_800_000, 106_500_000, 1000


def parse_args():
    """Optional: START END [suffix] to zoom into a sub-region."""
    import sys
    global START, END, OUT_PDF, OUT_PNG
    if len(sys.argv) >= 3:
        START = int(sys.argv[1])
        END = int(sys.argv[2])
        suffix = sys.argv[3] if len(sys.argv) > 3 else f'{START}_{END}'
        base = str(ROOT / f'output/chr14_{suffix}_methylation_tracks_corrected')
        OUT_PDF = base + '.pdf'
        OUT_PNG = base + '.png'


def load_data():
    df = pd.read_csv(MATRIX)
    mask = (df['chrom'] == CHROM) & (df['start'] >= START) & (df['end'] <= END)
    sub = df.loc[mask].copy()
    sub['g2_mean'] = sub[GROUP2].mean(axis=1, skipna=True) * 100
    sub['g3_mean'] = sub[GROUP3].mean(axis=1, skipna=True) * 100
    # number of covered G3/G2 samples per bin
    sub['g3_n'] = sub[GROUP3].notna().sum(axis=1)
    sub['g2_n'] = sub[GROUP2].notna().sum(axis=1)
    # group-level mean only where at least 3 of 5 (G3) / 5 of 8 (G2) covered;
    # bins with fewer samples are exploratory -> separate faint series
    sub['g3_mean_grp'] = sub['g3_mean'].where(sub['g3_n'] >= 3)
    sub['g2_mean_grp'] = sub['g2_mean'].where(sub['g2_n'] >= 5)
    sub['diff'] = sub['g2_mean_grp'] - sub['g3_mean_grp']
    x = (sub['start'].values + 500) / 1e6
    return sub, x


def count_genes_per_bin(chrom, start, end, bin_size=1000):
    n_bins = (end - start) // bin_size
    bin_starts = start + np.arange(n_bins) * bin_size
    cats = {
        'IG_C': ['IG_C_gene', 'IG_C_pseudogene'],
        'IG_J': ['IG_J_gene', 'IG_J_pseudogene'],
        'IG_D': ['IG_D_gene', 'IG_D_pseudogene'],
        'IG_V': ['IG_V_gene', 'IG_V_pseudogene'],
    }
    counts = {k: np.zeros(n_bins, dtype=int) for k in cats}
    with open(GTF) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            cols = line.strip().split('\t')
            if cols[0] != chrom or cols[2] != 'gene':
                continue
            gstart, gend = int(cols[3]), int(cols[4])
            if gend < start or gstart > end:
                continue
            m = re.search(r'gene_type "([^"]+)"', cols[8])
            gtype = m.group(1) if m else 'unknown'
            bs, be = max(gstart, start), min(gend, end)
            i0 = (bs - start) // bin_size
            i1 = min((be - start) // bin_size, n_bins - 1)
            for cat, types in cats.items():
                if gtype in types:
                    counts[cat][i0:i1 + 1] += 1
                    break
    return bin_starts, counts


def main():
    ROOT.joinpath('output').mkdir(parents=True, exist_ok=True)
    parse_args()
    sub, x = load_data()
    bin_starts, counts = count_genes_per_bin(CHROM, START, END, BIN)
    x_bar = (bin_starts + BIN // 2) / 1e6

    fig = plt.figure(figsize=(18.0 / 2.54, 8.0 / 2.54))
    gs = GridSpec(3, 1, height_ratios=[2.5, 2.8, 0.4],
                  hspace=0.08, left=0.15, right=0.96,
                  top=0.94, bottom=0.30)

    ax = fig.add_subplot(gs[0])
    ax_cov = fig.add_subplot(gs[1], sharex=ax)
    ax_bar = fig.add_subplot(gs[2], sharex=ax)

    # --- methylation tracks: NaN -> gaps ---
    for s in GROUP2:
        ax.plot(x, sub[s].values * 100, color='red', alpha=0.25, linewidth=0.6)
    for s in GROUP3:
        vals = sub[s].values * 100
        covered = sub[s].notna().values
        ax.plot(x, vals, color='blue', alpha=0.25, linewidth=0.6)
        # bright green at the baseline where this sample has no coverage
        miss = np.where(covered, np.nan, np.zeros(len(vals)))
        ax.plot(x, miss, color='lime', alpha=0.6, linewidth=1.2)
    ax.plot([], [], color='lime', linewidth=1.2, label='G3 no-coverage bin')
    # per-bin mean over covered samples only (NaN excluded); group-level series
    # masked to bins with enough covered samples
    ax.plot(x, sub['g2_mean'].values, color='darkred', linewidth=0.7,
            alpha=0.6, label='G2 covered-only mean')
    ax.plot(x, sub['g3_mean'].values, color='darkblue', linewidth=0.7,
            alpha=0.6, label='G3 covered-only mean')
    ax.plot(x, sub['g3_mean_grp'].values, color='darkblue', linewidth=0.7,
            alpha=0.9, label='G3 mean (>=3/5 samples)')
    ax.plot(x, sub['g2_mean_grp'].values, color='darkred', linewidth=0.7,
            alpha=0.9, label='G2 mean (>=5/8 samples)')
    # exploratory bins with only 1-2 covered G3 samples: faint dashed
    g3_low_n = sub['g3_mean'].where(sub['g3_n'] < 3)
    ax.plot(x, g3_low_n.values, color='darkblue', linewidth=0.7, alpha=0.35,
            linestyle=(0, (2, 2)), label='G3 mean, n<3 (exploratory)')
    # 10-bin rolling mean to suppress sparse-CpG noise
    g2_smooth = sub['g2_mean_grp'].rolling(window=10, min_periods=5, center=True).mean()
    g3_smooth = sub['g3_mean_grp'].rolling(window=10, min_periods=5, center=True).mean()
    g3_smooth_low = g3_low_n.rolling(window=10, min_periods=5, center=True).mean()
    ax.plot(x, g2_smooth.values, color='darkred', linewidth=1.4,
            label='G2 mean, 10kb smooth')
    ax.plot(x, g3_smooth.values, color='darkblue', linewidth=1.4,
            label='G3 mean, 10kb smooth')
    ax.plot(x, g3_smooth_low.values, color='darkblue', linewidth=1.0,
            alpha=0.4, linestyle=(0, (2, 2)))

    # mark dropout zone of the 3 low-coverage samples
    ax.axvspan(105.85, 106.25, color='lightgrey', alpha=0.25, zorder=0)
    ax.text(105.87, 96, 'CEMT0006/0095/0097 no data\n(alignment dropout); G3 here = '
            'n=2 (CEMT0094/0096),\ndashed line = exploratory',
            fontsize=5, fontname='Arial', va='top',
            color='dimgray')

    # gold: |delta| > 20% on corrected means
    for sb in sub.loc[sub['diff'].abs() > 20, 'start'].values:
        ax.axvspan(sb / 1e6, (sb + 1000) / 1e6, color='gold', alpha=0.18, zorder=0)

    ax.set_ylabel('Methylation (%)', fontname='Arial', fontsize=5)
    ax.set_title(f'{CHROM}:{START/1e6:.2f}-{END/1e6:.2f} Mb (corrected: missing = gap, not 0)',
                 fontname='Arial', fontsize=5)
    ax.set_xlim(START / 1e6, END / 1e6)
    ax.set_ylim(0, 100)
    ax.tick_params(axis='y', which='major', labelsize=5)
    ax.tick_params(axis='x', which='both', top=False, labeltop=False,
                   labelbottom=False)
    ax.set_yticks([0, 20, 40, 60, 80, 100])
    ax.grid(True, alpha=0.3, linewidth=0.3)

    # --- coverage strip for all 13 samples (G2 red labels, G3 blue) ---
    n = len(ALL_SAMPLES)
    ax_cov.set_ylim(-0.5, n - 0.5)
    ax_cov.set_yticks(range(n))
    ax_cov.set_yticklabels([f'{IGHV[s]} {s}' for s in ALL_SAMPLES],
                           fontsize=5, fontname='Arial')
    ax_cov.set_xlim(START / 1e6, END / 1e6)
    ax_cov.tick_params(axis='x', which='both', bottom=False, labelbottom=False,
                       top=False, labeltop=False)
    ax_cov.tick_params(axis='y', which='both', left=False, labelleft=True)
    for sp in ax_cov.spines.values():
        sp.set_visible(False)
    for i, s in enumerate(ALL_SAMPLES):
        cov = sub[s].notna().values.astype(float)
        color = 'red' if s in GROUP2 else 'blue'
        ax_cov.fill_between(x, i - 0.4, i - 0.4 + cov * 0.8,
                            color=color, alpha=0.6, linewidth=0)
    ax_cov.text(START / 1e6 + 0.003, -0.45, 'coverage (all samples)',
                fontsize=5, fontname='Arial', va='bottom', color='dimgray')

    # --- gene annotation strip ---
    colors = {'IG_C': '#e41a1c', 'IG_J': '#377eb8',
              'IG_D': '#4daf4a', 'IG_V': '#984ea3'}
    ax_bar.set_ylim(-0.5, 0.5)
    ax_bar.set_yticks([])
    ax_bar.set_xlim(START / 1e6, END / 1e6)
    ax_bar.text(START / 1e6 + 0.003, 0, 'genes',
                fontsize=5, fontname='Arial', ha='left', va='center')
    ax_bar.tick_params(axis='x', which='both', top=False, labeltop=False,
                       bottom=False, labelbottom=False)
    ax_bar.tick_params(axis='y', which='both', left=False, labelleft=False)
    for sp in ax_bar.spines.values():
        sp.set_visible(False)
    for cat in ['IG_C', 'IG_J', 'IG_D', 'IG_V']:
        for bi, v in enumerate(counts[cat]):
            if v > 0:
                ax_bar.add_patch(plt.Rectangle(
                    (x_bar[bi] - (BIN / 2) / 1e6, -0.35),
                    BIN / 1e6, 0.7,
                    facecolor=colors[cat], edgecolor='none', alpha=0.85))

    # x labels only on bottom
    ax_bar.set_xlabel('Position (Mb)', fontname='Arial', fontsize=5)
    ax_bar.tick_params(axis='x', which='both', bottom=True, labelbottom=True)
    xt = np.arange(np.ceil(START / 1e5) * 1e5, END + 1, 100000)
    ax_bar.set_xticks(xt / 1e6)
    ax_bar.set_xticklabels([f'{v/1e6:g}' for v in xt], fontsize=5)

    # Both legends stay within the fixed 18 x 8 cm canvas.
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.55, 0.075),
               ncol=4, fontsize=5, frameon=False, columnspacing=1.4,
               handlelength=2.0, labelspacing=0.4)
    gene_handles = [Patch(facecolor='none', edgecolor='none', label='Genes:')]
    gene_handles += [Patch(facecolor=color, edgecolor='none', alpha=0.85,
                           label=category) for category, color in colors.items()]
    fig.legend(handles=gene_handles, loc='lower center', bbox_to_anchor=(0.55, 0.005),
               ncol=5, fontsize=5, frameon=False, columnspacing=1.8)

    fig.savefig(OUT_PDF)
    fig.savefig(OUT_PNG, dpi=300)
    print(f'Saved: {OUT_PDF}')
    print(f'Saved: {OUT_PNG}')


if __name__ == '__main__':
    main()
