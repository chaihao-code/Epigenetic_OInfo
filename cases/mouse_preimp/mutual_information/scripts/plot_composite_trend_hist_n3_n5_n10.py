#!/usr/bin/env python3
"""
Composite figure: 3 rows (n_bins = 3, 5, 10) x 2 columns.
Left: MI trend for H3K4me3 anchor in distal_ctcf_bound_enhancer_purified_bins_AND_CpG_island.
Right: density histogram of pairwise / joint normalized MI ratios.
Figure 17 cm wide; each panel 6 cm high. Arial 5 pt, PDF fonttype 42.
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
BLOCK = "distal_ctcf_bound_enhancer_purified_bins_AND_CpG_island"
ANCHOR = "H3K4me3"
STAGES = ['zygote', '2cell', '4cell', '8cell', 'morula', 'ICM']

N_LIST = [3, 5, 10]

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


def load_trend_csv(n):
    f = PLOTS / f"n{n}" / f"{BLOCK}_{ANCHOR}_mi_trend.csv"
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
# Layout: exact cm geometry
# ----------------------------------------------------------------------------
CM = 1 / 2.54
FIG_W_CM = 17.0
LEFT_MARGIN = 1.6
LEFT_W = 5.9
GAP = 0.7
RIGHT_W = 8.3
RIGHT_MARGIN = FIG_W_CM - LEFT_MARGIN - LEFT_W - GAP - RIGHT_W  # 0.5

TOP = 0.5
PANEL_H = 6.0
ROW_GAP = 0.9
LEGEND_H = 1.6
FIG_H_CM = TOP + 3 * PANEL_H + 2 * ROW_GAP + LEGEND_H + 0.4

fig = plt.figure(figsize=(FIG_W_CM * CM, FIG_H_CM * CM))

axes_left, axes_right = [], []
for i in range(3):
    y0 = 1 - (TOP + i * (PANEL_H + ROW_GAP) + PANEL_H) / FIG_H_CM
    h = PANEL_H / FIG_H_CM
    ax_l = fig.add_axes([LEFT_MARGIN / FIG_W_CM, y0, LEFT_W / FIG_W_CM, h])
    ax_r = fig.add_axes([(LEFT_MARGIN + LEFT_W + GAP) / FIG_W_CM, y0,
                         RIGHT_W / FIG_W_CM, h])
    axes_left.append(ax_l)
    axes_right.append(ax_r)

x = np.arange(len(STAGES))

for row, n in enumerate(N_LIST):
    ax_l, ax_r = axes_left[row], axes_right[row]

    # ---------------- Left: MI trend ----------------
    df = load_trend_csv(n)
    targets = [t for t in df['target'].unique() if t != 'total'] + ['total']
    for target in targets:
        sub = df[df['target'] == target].set_index('stage').reindex(STAGES)
        y = sub['normalized_mi'].values
        is_total = target == 'total'
        ax_l.plot(
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
    ax_l.set_xticks(x)
    ax_l.set_xticklabels(STAGES, rotation=30, ha='right', fontsize=FS)
    ax_l.set_ylabel(r'$I(X;Y) / H(X)$', fontsize=FS)
    ax_l.tick_params(axis='both', labelsize=FS)
    ax_l.set_ylim(bottom=0)
    ax_l.grid(True, alpha=0.3)
    ax_l.text(-0.14, 1.02, f'n = {n}', transform=ax_l.transAxes,
              fontsize=6, fontweight='bold', ha='left', va='bottom')

    # ---------------- Right: histogram ----------------
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
    ax_r.bar(centers, dens(pair), width=widths * 0.95, color='lightskyblue',
             edgecolor='steelblue', alpha=0.7, align='center', linewidth=0.3)
    ax_r.bar(centers, dens(total), width=widths * 0.95, color='lightcoral',
             edgecolor='crimson', alpha=0.7, align='center', linewidth=0.3)
    ax_r.set_xscale('log')
    ax_r.set_xlabel(r'$I(X;Y) / H(X)$', fontsize=FS)
    ax_r.set_ylabel(r'Density ($\mathrm{d} \log_{10} x$)', fontsize=FS)
    ax_r.tick_params(axis='both', labelsize=FS)
    ax_r.grid(True, alpha=0.3, which='both')

    med_p, med_t = np.median(pair), np.median(total)
    ax_r.axvline(med_p, color='steelblue', linestyle='--', linewidth=0.9)
    ax_r.axvline(med_t, color='crimson', linestyle='--', linewidth=0.9)
    y_max = ax_r.get_ylim()[1]
    ax_r.text(med_p / 1.2, y_max * 0.90, f'pairwise\n{med_p:.3f}',
              color='steelblue', fontsize=FS, ha='right', va='top')
    ax_r.text(med_t / 1.1, y_max * 0.90, f'joint\n{med_t:.3f}',
              color='crimson', fontsize=FS, ha='right', va='top')

# ----------------------------------------------------------------------------
# Shared legend (targets are identical across rows)
# ----------------------------------------------------------------------------
handles, labels = axes_left[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', ncol=5, frameon=False,
           fontsize=FS, bbox_to_anchor=(0.5, 0.005), columnspacing=1.2,
           handlelength=1.8, handletextpad=0.6)

out_base = PLOTS / "composite_trend_hist_n3_n5_n10"
fig.savefig(f"{out_base}.png", dpi=300)
fig.savefig(f"{out_base}.pdf", format='pdf')
print(f"Saved: {out_base}.png/.pdf")
