#!/usr/bin/env python3
"""
Plot probability density histogram of MI ratios I(X;Y)/H(X) from pre-computed CSVs.
Density is computed with respect to d(log x) so that area under the curve = 1
on a log-scaled x-axis.
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

# ============================================================================
# Font setup
# ============================================================================
arial_found = False
for font_path in font_manager.findSystemFonts(fontpaths=None, fontext='ttf'):
    try:
        prop = font_manager.FontProperties(fname=font_path)
        if prop.get_name() == 'Arial' or 'arial' in Path(font_path).name.lower():
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

BASE_DIR = Path(__file__).resolve().parents[1]
FS = 5


def load_ratios_from_csvs(n_bins: int):
    """Load normalized_mi values for pair and total separately."""
    csv_dir = BASE_DIR / "plots" / f"n{n_bins}"
    files = list(csv_dir.glob("*_mi_trend.csv"))
    print(f"n_bins={n_bins}: found {len(files)} CSV files")
    
    pair_ratios = []
    total_ratios = []
    
    for f in files:
        df = pd.read_csv(f, usecols=['target', 'normalized_mi'])
        df = df.dropna()
        
        pair_vals = df[df['target'] != 'total']['normalized_mi'].values
        total_vals = df[df['target'] == 'total']['normalized_mi'].values
        
        pair_vals = pair_vals[pair_vals > 0]
        total_vals = total_vals[total_vals > 0]
        
        pair_ratios.extend(pair_vals)
        total_ratios.extend(total_vals)
    
    pair_ratios = np.array(pair_ratios, dtype=np.float64)
    total_ratios = np.array(total_ratios, dtype=np.float64)
    
    print(f"  Pair ratios:   {len(pair_ratios):,}  | median={np.median(pair_ratios):.6f}, max={pair_ratios.max():.6f}")
    print(f"  Total ratios:  {len(total_ratios):,}  | median={np.median(total_ratios):.6f}, max={total_ratios.max():.6f}")
    
    return pair_ratios, total_ratios


def plot_log_density_histogram(pair_ratios, total_ratios, n_bins, output_path):
    fig_w = 8 / 2.54
    fig_h = 7 / 2.54
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    
    # Common log-spaced bins
    all_vals = np.concatenate([pair_ratios, total_ratios])
    lower = np.percentile(all_vals, 0.1)
    upper = all_vals.max()
    bins = np.logspace(np.log10(lower), np.log10(upper), 80)
    log_bins = np.log10(bins)
    d_log_x = np.diff(log_bins)[0]
    
    def compute_density(values):
        counts, _ = np.histogram(values, bins=bins)
        density = counts / (len(values) * d_log_x)
        return density
    
    dens_pair = compute_density(pair_ratios)
    dens_total = compute_density(total_ratios)
    
    bin_centers = 10 ** ((log_bins[:-1] + log_bins[1:]) / 2)
    bin_widths = np.diff(bins)
    
    ax.bar(bin_centers, dens_pair, width=bin_widths * 0.95, color='lightskyblue',
           edgecolor='steelblue', alpha=0.7, align='center')
    ax.bar(bin_centers, dens_total, width=bin_widths * 0.95, color='lightcoral',
           edgecolor='crimson', alpha=0.7, align='center')
    
    ax.set_xscale('log')
    ax.set_xlabel(r'$I(X;Y) / H(X)$', fontsize=FS)
    ax.set_ylabel(r'Density ($\mathrm{d} \log_{10} x$)', fontsize=FS)
    ax.grid(True, alpha=0.3, which='both')
    ax.tick_params(axis='both', which='major', labelsize=FS)
    
    # Median vertical lines with colored labels
    med_pair = np.median(pair_ratios)
    med_total = np.median(total_ratios)
    
    ax.axvline(med_pair, color='steelblue', linestyle='--', linewidth=1.0)
    ax.axvline(med_total, color='crimson', linestyle='--', linewidth=1.0)
    
    y_max = ax.get_ylim()[1]
    ax.text(med_pair, y_max * 0.92, f'pairwise\nmedian={med_pair:.3f}',
            color='steelblue', fontsize=FS, ha='center', va='top',
            fontweight='bold')
    ax.text(med_total, y_max * 0.92, f'joint\nmedian={med_total:.3f}',
            color='crimson', fontsize=FS, ha='center', va='top',
            fontweight='bold')
    
    ax.set_title(f'n_bins={n_bins}', fontsize=FS)
    
    plt.tight_layout()
    
    png_path = output_path.with_suffix('.png')
    pdf_path = output_path.with_suffix('.pdf')
    plt.savefig(png_path, dpi=300)
    plt.savefig(pdf_path, format='pdf', bbox_inches='tight')
    plt.close()
    print(f"  Saved: {png_path}")
    print(f"  Saved: {pdf_path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-bins", type=int, required=True, choices=[3, 5, 10])
    args = parser.parse_args()
    
    pair_ratios, total_ratios = load_ratios_from_csvs(args.n_bins)
    
    out_dir = BASE_DIR / "plots"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"mi_ratio_histogram_n{args.n_bins}"
    
    plot_log_density_histogram(pair_ratios, total_ratios, args.n_bins, out_path)


if __name__ == "__main__":
    main()
