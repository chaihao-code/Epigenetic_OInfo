#!/usr/bin/env python3
"""Signed delta-beta versus nominal P, highlighting BH-FDR significance."""
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/bcells-dmr-mpl')
os.environ.setdefault('XDG_CACHE_HOME', '/tmp/bcells-dmr-cache')
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / 'output'


def main():
    source = OUTPUT / 'all_bins_1kb_cemt_g2_g3_limma_results.csv'
    if not source.exists():
        source = source.with_suffix('.csv.gz')
    data = pd.read_csv(source, usecols=['chrom', 'start', 'end', 'delta_beta', 'pvalue', 'qvalue'])
    valid = np.isfinite(data.delta_beta) & np.isfinite(data.pvalue) & data.pvalue.between(0, 1)
    plotted = data.loc[valid].copy()
    plotted['minus_log10_p'] = -np.log10(plotted.pvalue.clip(lower=np.finfo(float).tiny))
    fdr = plotted.qvalue.lt(.05)
    hits = fdr & plotted.delta_beta.abs().gt(.2)
    only_fdr = fdr & ~hits
    high_r = hits & plotted.delta_beta.gt(0)
    high_s = hits & plotted.delta_beta.lt(0)
    assert len(data) == 2697259 and int(fdr.sum()) == 29
    assert int(high_r.sum()) == 11 and int(high_s.sum()) == 4
    # One global BH family gives a data-specific raw-P boundary; this is not P=.05.
    p_boundary = float(plotted.loc[fdr, 'pvalue'].max())
    boundary_y = -np.log10(p_boundary)
    font = Path.home() / '.local/share/fonts/arial.ttf'
    if font.exists():
        font_manager.fontManager.addfont(str(font))
        plt.rcParams['font.family'] = 'Arial'
    plt.rcParams.update({'font.size': 5, 'axes.titlesize': 5, 'axes.labelsize': 5,
                         'xtick.labelsize': 5, 'ytick.labelsize': 5,
                         'legend.fontsize': 5, 'pdf.fonttype': 42, 'ps.fonttype': 42})
    fig, ax = plt.subplots(figsize=(8 / 2.54, 8 / 2.54))
    background = plotted.loc[~fdr]
    ax.scatter(background.delta_beta, background.minus_log10_p, s=.35,
               c='#858585', alpha=.12, linewidths=0, rasterized=True, zorder=1)
    for mask, color, marker in [(only_fdr, '#e5a335', 'o'),
                                (high_s, '#3b74af', 'o'), (high_r, '#c83b3b', 'o')]:
        subset = plotted.loc[mask]
        ax.scatter(subset.delta_beta, subset.minus_log10_p, s=10, color=color,
                   edgecolors='white', linewidths=.45, marker=marker, zorder=4)
    ax.axvline(0, color='#999999', linewidth=.65, zorder=2)
    for value in (-.2, .2):
        ax.axvline(value, color='#666666', linestyle='--', linewidth=.85, zorder=2)
    ax.axhline(boundary_y, color='#6a3d9a', linestyle='--', linewidth=.9, zorder=2)
    limit = max(.3, np.ceil(float(plotted.delta_beta.abs().max())*10)/10)
    ax.set_xlim(-limit-.025, limit+.025)
    ax.set_ylim(0, max(10, float(plotted.minus_log10_p.max())+.6))
    ax.set_xlabel('Δβ = mean methylation in CLL-R − CLL-S')
    ax.set_ylabel(r'$-\log_{10}$(nominal P-value)')
    ax.set_title('Differential methylation in 1 kb windows: CLL-R versus CLL-S\n'
                 'limma on M-values, adjusting age and sex', fontsize=5, pad=5)
    legend = [
        Line2D([], [], marker='o', linestyle='none', color='#999999', markersize=3,
               label=f'FDR ≥ 0.05 / unavailable (n={len(plotted)-int(fdr.sum()):,})'),
        Line2D([], [], marker='o', linestyle='none', color='#e5a335', markersize=3,
               label=f'FDR < 0.05, |Δβ| ≤ 0.2 (n={int(only_fdr.sum())})'),
        Line2D([], [], marker='o', linestyle='none', color='#c83b3b', markersize=3,
               label=f'R higher: FDR < 0.05, Δβ > 0.2 (n={int(high_r.sum())})'),
        Line2D([], [], marker='o', linestyle='none', color='#3b74af', markersize=3,
               label=f'S higher: FDR < 0.05, Δβ < −0.2 (n={int(high_s.sum())})'),
        Line2D([], [], linestyle='--', color='#6a3d9a', linewidth=.9,
               label=f'BH FDR boundary: P ≤ {p_boundary:.2g}'),
    ]
    fig.legend(handles=legend, loc='lower left', bbox_to_anchor=(.14, .07), fontsize=5, frameon=False,
              facecolor='white', framealpha=.95, edgecolor='#dddddd')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.text(.16, .025,
             f'{len(plotted):,} windows plotted; {len(data)-len(plotted)} missing P or Δβ.',
             fontsize=5, ha='left', va='bottom')
    fig.subplots_adjust(left=.16, right=.96, top=.89, bottom=.40)
    for extension in ('pdf', 'png'):
        target = OUTPUT / ('cemt_R_vs_S_dmr_volcano.'+extension)
        fig.savefig(target, dpi=300)
        print('Saved', target, flush=True)
    plt.close(fig)
    plotted.loc[fdr].sort_values('qvalue').to_csv(OUTPUT/'volcano_FDR005_windows.tsv', sep='\t', index=False)
    metadata = dict(figure_size_cm=[8, 8], font_size_pt=5, source=str(source), x='Signed delta_beta = Group2/R - Group3/S',
                    y='-log10 nominal P-value from age/sex-adjusted limma M-value model',
                    n_input=len(data), n_plotted=len(plotted), n_excluded=len(data)-len(plotted),
                    n_FDR005=int(fdr.sum()), n_FDR005_effect_large=int(hits.sum()),
                    n_R_higher=int(high_r.sum()), n_S_higher=int(high_s.sum()),
                    BH_family='All windows retained by the original analysis, not only effect-filtered windows',
                    raw_P_boundary_for_FDR005=p_boundary,
                    horizontal_line='Largest observed raw P among windows with BH q<0.05; data-specific, not universal',
                    background='All non-FDR windows plotted, no sampling; background rasterized in PDF')
    (OUTPUT/'volcano_metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
    print(json.dumps(metadata, indent=2), flush=True)


if __name__ == '__main__':
    main()
