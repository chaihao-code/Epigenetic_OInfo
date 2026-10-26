#!/usr/bin/env python3
"""
Composite figure (transposed layout): 3 columns (n_bins = 3, 5, 10) x 3 rows.
Row 1 (top): density histogram of pairwise / joint normalized MI ratios.
Row 2: MI trend, H3K4me3 anchor in distal_ctcf_bound_enhancer_purified_bins_AND_CpG_island.
Row 3: MI trend, H3K9me3 anchor in distal_not_ctcf_bound_enhancer_purified_bins_AND_SINE.
Each column 5.5 cm wide. Arial 5 pt, PDF fonttype 42.
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

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.font_manager import FontProperties, fontManager
from matplotlib import rcParams

# ----------------------------------------------------------------------------
# Font
# ----------------------------------------------------------------------------
font_path = _release_path('@ROOT@/data/fonts/arial.ttf')
if Path(font_path).exists():
    fontManager.addfont(font_path)
    rcParams['font.family'] = FontProperties(fname=font_path).get_name()
else:
    rcParams['font.family'] = 'DejaVu Sans'
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42

FS = 5

BASE = Path(__file__).resolve().parents[1]
PLOTS = BASE / "plots"
STAGES = ['zygote', '2cell', '4cell', '8cell', 'morula', 'ICM']

N_LIST = [3, 5, 10]

TREND_PANELS = [
    ("distal_ctcf_bound_enhancer_purified_bins_AND_CpG_island", "H3K4me3",
     'H3K4me3 | distal CTCF+ enhancer / CpG island'),
    ("distal_not_ctcf_bound_enhancer_purified_bins_AND_SINE", "H3K9me3",
     'H3K9me3 | distal CTCF- enhancer / SINE'),
]

TARGET_COLORS = {
    'ATAC': '#999999',
    'CpG_methyl': '#A65628',
    'H2A.Z': '#F2C14E',
    'H3K27ac': '#E41A1C',
    'H3K27me3': '#377EB8',
    'H3K4me1': '#4DAF4A',
    'H3K4me3': '#984EA3',
    'H3K36me3': '#72B7B2',
    'H3K9me3': '#FF7F00',
    'total': '#000000',
}


def load_trend_csv(n, block, anchor):
    f = PLOTS / f"n{n}" / f"{block}_{anchor}_mi_trend.csv"
    return pd.read_csv(f)


def load_ratios(n):
    pair, total = [], []
    for f in (PLOTS / f"n{n}").glob("*_mi_trend.csv"):
        df = pd.read_csv(f, usecols=['target', 'normalized_mi']).dropna()
        p = df.loc[df['target'] != 'total', 'normalized_mi'].values
        t = df.loc[df['target'] == 'total', 'normalized_mi'].values
        pair.extend(p[p > 0])
        total.extend(t[t > 0])
    return np.asarray(pair), np.asarray(total)


# ----------------------------------------------------------------------------
# Layout: exact cm geometry, 3 columns x 3 rows
# ----------------------------------------------------------------------------
CM = 1 / 2.54
COL_W = 5.5
LEFT_MARGIN = 1.7
COL_GAP = 0.5
RIGHT_MARGIN = 0.3
FIG_W_CM = LEFT_MARGIN + 3 * COL_W + 2 * COL_GAP + RIGHT_MARGIN  # 18.8

TOP = 0.7            # room for "n = x" labels
ROW_H = 4.6
ROW_GAP = 0.8
LEGEND_H = 1.4
N_ROWS = 3
FIG_H_CM = TOP + N_ROWS * ROW_H + (N_ROWS - 1) * ROW_GAP + LEGEND_H + 0.2

fig = plt.figure(figsize=(FIG_W_CM * CM, FIG_H_CM * CM))

# axes[row][col], row 0 = histogram, rows 1..2 = trend panels
axes = []
for i in range(N_ROWS):
    y0 = 1 - (TOP + i * (ROW_H + ROW_GAP) + ROW_H) / FIG_H_CM
    row_axes = []
    for j in range(3):
        x0 = (LEFT_MARGIN + j * (COL_W + COL_GAP)) / FIG_W_CM
        row_axes.append(fig.add_axes([x0, y0, COL_W / FIG_W_CM, ROW_H / FIG_H_CM]))
    axes.append(row_axes)

x = np.arange(len(STAGES))

for col, n in enumerate(N_LIST):
    # ---------------- Row 0: histogram ----------------
    ax_h = axes[0][col]
    pair, total = load_ratios(n)
    all_vals = np.concatenate([pair, total])
    lower = np.percentile(all_vals, 0.1)
    upper = all_vals.max()
    bins = np.logspace(np.log10(lower), np.log10(upper), 80)
    log_bins = np.log10(bins)
    d_log_x = np.diff(log_bins)[0]

    def dens(v):
        c, _ = np.histogram(v, bins=bins)
        return c / (len(v) * d_log_x)

    centers = 10 ** ((log_bins[:-1] + log_bins[1:]) / 2)
    widths = np.diff(bins)
    ax_h.bar(centers, dens(pair), width=widths * 0.95, color='lightskyblue',
             edgecolor='steelblue', alpha=0.7, align='center', linewidth=0.3)
    ax_h.bar(centers, dens(total), width=widths * 0.95, color='lightcoral',
             edgecolor='crimson', alpha=0.7, align='center', linewidth=0.3)
    ax_h.set_xscale('log')
    ax_h.set_title(f'n = {n}', fontsize=6, fontweight='bold', pad=3)
    ax_h.tick_params(axis='both', labelsize=FS)
    ax_h.grid(True, alpha=0.3, which='both')
    if col == 0:
        ax_h.set_ylabel(r'Density ($\mathrm{d} \log_{10} x$)', fontsize=FS)
    else:
        ax_h.set_yticklabels([])

    med_p, med_t = np.median(pair), np.median(total)
    ax_h.axvline(med_p, color='steelblue', linestyle='--', linewidth=0.9)
    ax_h.axvline(med_t, color='crimson', linestyle='--', linewidth=0.9)
    y_max = ax_h.get_ylim()[1]
    ax_h.text(med_p / 1.25, y_max * 0.90, f'pairwise\n{med_p:.3f}',
              color='steelblue', fontsize=FS, ha='right', va='top')
    ax_h.text(med_t / 1.1, y_max * 0.90, f'joint\n{med_t:.3f}',
              color='crimson', fontsize=FS, ha='right', va='top')

    # ---------------- Rows 1..2: MI trends ----------------
    for row, (block, anchor, row_label) in enumerate(TREND_PANELS, start=1):
        ax_t = axes[row][col]
        df = load_trend_csv(n, block, anchor)
        targets = [t for t in df['target'].unique() if t != 'total'] + ['total']
        for target in targets:
            sub = df[df['target'] == target].set_index('stage').reindex(STAGES)
            y = sub['normalized_mi'].values
            is_total = target == 'total'
            ax_t.plot(
                x, y,
                label=target,
                color=TARGET_COLORS.get(target, None),
                linestyle='--' if is_total else '-',
                linewidth=1.2 if is_total else 0.9,
                marker='o', markersize=2.5,
                markerfacecolor='white', markeredgewidth=0.7,
                markeredgecolor=TARGET_COLORS.get(target, None),
                zorder=5 if is_total else 3,
            )
        ax_t.set_xticks(x)
        ax_t.set_xticklabels(STAGES, rotation=30, ha='right', fontsize=FS)
        ax_t.set_ylim(0, 1.05)
        ax_t.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
        ax_t.tick_params(axis='both', labelsize=FS)
        ax_t.grid(True, alpha=0.3)
        if col == 0:
            ax_t.set_ylabel(r'$I(X;Y) / H(X)$', fontsize=FS)
        else:
            ax_t.set_yticklabels([])

# Row labels (rotated, far left)
row_labels = ['all 40 blocks'] + [lbl for (_, _, lbl) in TREND_PANELS]
for i, lbl in enumerate(row_labels):
    y0 = 1 - (TOP + i * (ROW_H + ROW_GAP) + ROW_H / 2) / FIG_H_CM
    fig.text(0.45 / FIG_W_CM, y0, lbl, rotation=90, ha='center', va='center',
             fontsize=FS)

# Histogram xlabel on bottom of row 0 (in the gap below the row)
fig.text((LEFT_MARGIN + 1.5 * (COL_W + COL_GAP)) / FIG_W_CM,
         1 - (TOP + ROW_H + 0.55) / FIG_H_CM,
         r'$I(X;Y) / H(X)$', ha='center', va='top', fontsize=FS)

# ----------------------------------------------------------------------------
# Shared legend (targets are identical across panels)
# ----------------------------------------------------------------------------
handles, labels = axes[1][0].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', ncol=5, frameon=False,
           fontsize=FS, bbox_to_anchor=(0.5, 0.005), columnspacing=1.2,
           handlelength=1.8, handletextpad=0.6)

out_base = PLOTS / "composite_trend_hist_n3_n5_n10_transposed"
fig.savefig(f"{out_base}.png", dpi=300)
fig.savefig(f"{out_base}.pdf", format='pdf')
print(f"Saved: {out_base}.png/.pdf")
