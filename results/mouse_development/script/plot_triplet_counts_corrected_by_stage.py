#!/usr/bin/env python3
"""Plot BH35-significant negative triplets by stage and genomic block.

Adapted from project/permutation_test/2026-09-28/scripts/
plot_triplet_counts_corrected_by_stage.py. Input is the consolidated 10,000-permutation FDR table:
35 p_down tests corrected separately within every sample x block,
then require observed_hoi < 0 and q <= 0.05. Stage values average
counts across available tissues. Missing stages are omitted.
Override input run and plot directory with HOI_NULL_DIR and HOI_OUT_DIR.
"""

import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/mouse-development-mpl')
os.environ.setdefault('XDG_CACHE_HOME', '/tmp/mouse-development-cache')
os.makedirs(os.environ['MPLCONFIGDIR'], exist_ok=True)
os.makedirs(os.environ['XDG_CACHE_HOME'], exist_ok=True)
import csv
import glob
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.text import Text

import os as _os
ARIAL_CANDIDATES = [
    '/usr/share/fonts/truetype/msttcorefonts/Arial.ttf',
    _os.path.expanduser('~/.local/share/fonts/arial.ttf'),
    _os.path.expanduser('~/.local/share/fonts/ARIAL.TTF'),
]
ARIAL_PATH = next((p for p in ARIAL_CANDIDATES if _os.path.exists(p)), None)
if ARIAL_PATH is None:
    raise RuntimeError('Arial font not found in any known location')
arial_prop = FontProperties(fname=ARIAL_PATH, size=5)
from matplotlib import font_manager
font_manager.fontManager.addfont(ARIAL_PATH)
plt.rcParams['font.family'] = 'Arial'

plt.rcParams['font.size'] = 5
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

NULL_DIR = os.environ.get(
    'HOI_NULL_DIR',
    '/lustre/home/2101110354/sequence_entropy/project/mouse_development/result_oinfo_n_bins_5_triplet_perm10000')
OUT_DIR = os.environ.get(
    'HOI_OUT_DIR',
    '/lustre/home/2101110354/sequence_entropy/project/mouse_development/synergistic_triplets')
os.makedirs(OUT_DIR, exist_ok=True)

STAGE_ORDER = ['zygote', '2cell', '4cell', '8cell', 'morula', 'ICM', 'TE',
               'E6.5', 'E7.0', 'E7.5',
               'E11.5', 'E12.5', 'E13.5', 'E14.5', 'E15.5', 'E16.5', 'P0']

REGION_ORDER = [
    'promoter',
    'proximal_ctcf_bound_enhancer_purified_bins',
    'proximal_not_ctcf_bound_enhancer_purified_bins',
    'distal_ctcf_bound_enhancer_purified_bins',
    'distal_not_ctcf_bound_enhancer_purified_bins',
    'exon_purified_bins',
    'intron_purified_bins',
    'None',
]
FEATURE_ORDER = ['CpG_island', 'LINE', 'SINE', 'LTR', 'None']
ORDERED_REGIONS = [f'{r}_AND_{f}' for r in REGION_ORDER for f in FEATURE_ORDER]
REGION_LABELS = ['promoter', 'proximal\nCTCF+', 'proximal\nCTCF-',
                 'distal\nCTCF+', 'distal\nCTCF-', 'exon', 'intron', 'None']

REPRESSIVE_MARKS = ['H3K9me3', 'H3K27me3', 'H3K36me3']
P_CUT = 0.05


def has_dna_methylation(names):
    return 'CpG_methyl' in names or 'WGBS' in names


def has_repressive_mark(names):
    return any(m in names for m in REPRESSIVE_MARKS)


def pattern_a(names):
    return has_dna_methylation(names) and has_repressive_mark(names)


def pattern_b(names):
    return not has_dna_methylation(names) and has_repressive_mark(names)


def pattern_c(names):
    return not pattern_a(names) and not pattern_b(names)


def load_counts():
    """Count negative BH35 discoveries, retaining sign-only counts for reference."""
    source = Path(OUT_DIR) / 'all_triplets_fdr_10000.tsv.gz'
    data = pd.read_csv(source, sep='\t')
    assert data['order'].eq(3).all()
    assert data['n_perm'].eq(10000).all()
    assert data.groupby(['sample', 'block']).size().eq(35).all()
    counts = {}
    for (sample, region), gdf in data.groupby(['sample', 'block']):
        for cls, matcher in (('A', pattern_a), ('B', pattern_b), ('C', pattern_c)):
            sub = gdf[gdf['triplet'].map(matcher)]
            negative = sub['observed_hoi'] < 0
            significant = negative & (sub['q_down_sample_block_bh'] <= P_CUT)
            counts[(sample, region, cls)] = (int(significant.sum()), int(negative.sum()))
    assert sum(v[0] for v in counts.values()) == int(((data['observed_hoi'] < 0) & (data['q_down_sample_block_bh'] <= P_CUT)).sum())
    return counts


def plot_heatmap(stages, values, title, cbar_label, out_base):
    n_rows = len(stages)
    n_cols = len(ORDERED_REGIONS)
    mat = np.array(values, dtype=float)

    # fixed inner-matrix size: 11 cm x 3.5 cm
    axes_width_inch = 11 / 2.54
    axes_height_inch = 3.5 / 2.54
    left_margin, right_margin = 0.13, 0.90
    top_margin, bottom_margin = 0.80, 0.28
    width_inch = axes_width_inch / (right_margin - left_margin)
    height_inch = axes_height_inch / (top_margin - bottom_margin)

    fig, ax = plt.subplots(figsize=(width_inch, height_inch))
    fig.subplots_adjust(bottom=bottom_margin, left=left_margin,
                        right=right_margin, top=top_margin)

    cmap = 'Reds'
    vmin, vmax = 0, np.nanmax(mat)
    fmt = '{v:.1f}'
    threshold = vmax * 0.55 if vmax > 0 else 0

    im = ax.imshow(mat, aspect='auto', cmap=cmap, vmin=vmin, vmax=vmax,
                   interpolation='nearest',
                   extent=(-0.5, n_cols - 0.5, n_rows - 0.5, -0.5))

    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(stages, fontsize=5, fontproperties=arial_prop)
    ax.set_xticks([])
    ax.set_xticklabels([])

    for i in range(n_rows):
        for j in range(n_cols):
            val = mat[i, j]
            if np.isnan(val):
                continue
            color = 'white' if val > threshold else 'black'
            txt = f'{val:.0f}' if val >= 10 else f'{val:.1f}'
            ax.text(j, i, txt, ha='center', va='center',
                    fontsize=3, color=color, fontproperties=arial_prop)

    cbar = fig.colorbar(im, ax=ax, shrink=1.0, pad=0.01)
    cbar.set_label(cbar_label, fontsize=5, fontproperties=arial_prop)
    cbar.ax.tick_params(labelsize=5)
    for label in cbar.ax.get_yticklabels():
        label.set_fontproperties(arial_prop)

    for i in range(1, n_cols):
        if ORDERED_REGIONS[i - 1].split('_AND_')[0] != ORDERED_REGIONS[i].split('_AND_')[0]:
            ax.axvline(i - 0.5, color='black', linewidth=0.8, alpha=0.3)

    for r_idx, r in enumerate(REGION_ORDER):
        indices = [i for i, reg in enumerate(ORDERED_REGIONS)
                   if reg.startswith(f'{r}_AND_')]
        if not indices:
            continue
        left = indices[0] - 0.4
        right = indices[-1] + 0.4
        ax.annotate('', xy=(right, -0.08), xytext=(left, -0.08),
                    xycoords=ax.get_xaxis_transform(),
                    textcoords=ax.get_xaxis_transform(),
                    arrowprops=dict(arrowstyle='<->', color='gray', lw=1.5))
        ax.text((left + right) / 2, -0.10, REGION_LABELS[r_idx],
                ha='center', va='top', fontsize=5, color='black',
                transform=ax.get_xaxis_transform())

    feature_y = -0.30
    for i, reg in enumerate(ORDERED_REGIONS):
        f = reg.split('_AND_', 1)[1]
        label = f.replace('_', ' ')
        ax.text(i, feature_y, label, ha='center', va='top', fontsize=5,
                color='black', transform=ax.get_xaxis_transform(), rotation=90)

    ax.set_title(title, fontsize=5, fontproperties=arial_prop, pad=5)
    ax.set_xlim(-0.5, n_cols - 0.5)

    out_pdf = os.path.join(OUT_DIR, f'{out_base}.pdf')
    out_png = os.path.join(OUT_DIR, f'{out_base}.png')
    fig.savefig(out_pdf, format='pdf', dpi=300)
    fig.savefig(out_png, format='png', dpi=300)
    plt.close(fig)
    print(f'[{out_base}] saved {out_pdf}')


def plot_total_by_stage(base_stages, totals, class_totals, out_base):
    """Line plot of total synergistic triplet count vs developmental stage.

    totals: {stage: mean total count per sample}
    class_totals: {cls: {stage: mean class count per sample}}
    """
    width_inch, height_inch = 9 / 2.54, 6 / 2.54
    fig, ax = plt.subplots(figsize=(width_inch, height_inch))
    fig.subplots_adjust(left=0.14, right=0.97, top=0.86, bottom=0.22)

    x = np.arange(len(base_stages))
    class_colors = {'A': '#1f77b4', 'B': '#2ca02c', 'C': '#ff7f0e'}
    for cls in ('A', 'B', 'C'):
        ys = [class_totals[cls].get(s, np.nan) for s in base_stages]
        ax.plot(x, ys, color=class_colors[cls], linewidth=0.8,
                marker='o', markersize=2, label=f'Triplet {cls}')
    y_total = [totals.get(s, np.nan) for s in base_stages]
    ax.plot(x, y_total, color='black', linewidth=1.2,
            marker='o', markersize=2.5, label='Total')

    ax.set_xticks(x)
    ax.set_xticklabels(base_stages, rotation=45, ha='right',
                       fontsize=5, fontproperties=arial_prop)
    for label in ax.get_yticklabels():
        label.set_fontproperties(arial_prop)
    ax.set_ylabel('Mean synergistic triplets per sample',
                  fontsize=5, fontproperties=arial_prop)
    ax.set_title(f'Total synergistic triplets along development '
                 f'(summed over 40 regions, q_down <= {P_CUT}, O-info < 0)',
                 fontsize=5, fontproperties=arial_prop, pad=14)
    ax.legend(fontsize=4.5, prop=arial_prop, frameon=False, ncol=4,
              loc='lower center', bbox_to_anchor=(0.5, 1.01),
              handlelength=1.2, columnspacing=0.8, borderaxespad=0)

    out_pdf = os.path.join(OUT_DIR, f'{out_base}.pdf')
    out_png = os.path.join(OUT_DIR, f'{out_base}.png')
    fig.savefig(out_pdf, format='pdf', dpi=300)
    fig.savefig(out_png, format='png', dpi=300)
    plt.close(fig)
    print(f'[{out_base}] saved {out_pdf}')


def plot_stacked(stages, matrices, base_of, *, write_summary=True):
    """Recreate the three vertically stacked A/B/C heatmaps."""
    titles = {
        'A': 'Triplet A: DNA methylation + (H3K9me3 or H3K27me3) + other',
        'B': 'Triplet B: (H3K9me3 or H3K27me3) + other (no DNA methylation)',
        'C': 'Triplet C: neither A nor B',
    }
    width_cm, height_cm = 18.0, 19.5
    top_margin_cm, bottom_margin_cm, panel_gap_cm = 1.45, 2.65, .65
    panel_height_cm = (height_cm-top_margin_cm-bottom_margin_cm-2*panel_gap_cm)/3
    fig, axes = plt.subplots(3, 1, figsize=(width_cm/2.54, height_cm/2.54), dpi=300)
    fig.subplots_adjust(left=0.13, right=0.94,
                        top=1-top_margin_cm/height_cm, bottom=bottom_margin_cm/height_cm,
                        hspace=panel_gap_cm/panel_height_cm)
    meta = json.loads((Path(OUT_DIR) / 'analysis_metadata.json').read_text())
    fig.suptitle('10,000 permutations; negative O-information: BH within each sample x region (35 triplets), q <= 0.05\n'
                 f"{meta['n_samples']}/{82} samples available; mean count across tissues per stage",
                 fontsize=8, y=1-.15/height_cm)
    labels = [f"{stage} (n={sum(v == stage for v in base_of.values())})" for stage in stages]
    cell_texts = []
    color_limits = {}
    for ax, cls in zip(axes, 'ABC'):
        mat = matrices[cls]
        vmax = max(1.0, float(np.nanmax(mat)))
        color_limits[cls] = [0, vmax]
        im = ax.imshow(mat, aspect='auto', cmap='Reds', vmin=0, vmax=vmax, interpolation='nearest')
        ax.set_yticks(range(len(stages)), labels=labels, fontsize=8)
        ax.set_xticks([])
        ax.set_title(titles[cls], fontsize=8, pad=8)
        for i in range(len(stages)):
            for j in range(len(ORDERED_REGIONS)):
                value = mat[i, j]
                if np.isfinite(value):
                    cell_texts.append(ax.text(j, i, f'{value:.0f}' if value >= 10 else f'{value:.1f}',
                            ha='center', va='center', fontsize=6,
                            color='white' if value > vmax * .55 else 'black'))
        for j in range(5, 40, 5):
            ax.axvline(j - .5, color='black', lw=.8, alpha=.4)
        cbar = fig.colorbar(im, ax=ax, fraction=.025, pad=.01)
        cbar.ax.tick_params(labelsize=8)
        cbar.set_label('Mean count', fontsize=8)
    ax = axes[-1]
    panel_height_cm = ax.get_position().height * height_cm
    region_texts, feature_texts = [], []
    for i, label in enumerate(REGION_LABELS):
        arrow_y = -.15/panel_height_cm
        ax.annotate('', xy=(5*i+4.3, arrow_y), xytext=(5*i-.3, arrow_y),
                    xycoords=ax.get_xaxis_transform(), arrowprops=dict(arrowstyle='<->', color='gray', lw=1.5))
        region_texts.append(ax.text(5*i+2, -.28/panel_height_cm, label,
                            transform=ax.get_xaxis_transform(), ha='center', va='top', fontsize=8))
    for j, region in enumerate(ORDERED_REGIONS):
        feature_texts.append(ax.text(j, -1.05/panel_height_cm, region.split('_AND_')[1].replace('_', ' '),
                transform=ax.get_xaxis_transform(), ha='center', va='top', rotation=90, fontsize=8))
    cells = set(cell_texts)
    for text in fig.findobj(Text):
        if text not in cells:
            text.set_fontsize(8)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for text in fig.findobj(Text):
        if text.get_visible() and text.get_text():
            assert text.get_fontsize() == (6 if text in cells else 8), text.get_text()
            box = text.get_window_extent(renderer)
            assert box.x0 >= 0 and box.y0 >= 0 and box.x1 <= fig.bbox.width and box.y1 <= fig.bbox.height, text.get_text()
            if text in cells:
                panel_box = text.axes.get_window_extent(renderer)
                assert box.width < panel_box.width/len(ORDERED_REGIONS)
                assert box.height < panel_box.height/len(stages)
    for first, second in zip(region_texts, region_texts[1:]):
        assert not first.get_window_extent(renderer).overlaps(second.get_window_extent(renderer))
    for ax in axes:
        labels_y = ax.get_yticklabels()
        for first, second in zip(labels_y, labels_y[1:]):
            assert not first.get_window_extent(renderer).overlaps(second.get_window_extent(renderer)), 'Overlapping stage labels'
    assert not fig._suptitle.get_window_extent(renderer).overlaps(axes[0].title.get_window_extent(renderer))
    assert min(t.get_window_extent(renderer).y0 for t in region_texts) > max(
        t.get_window_extent(renderer).y1 for t in feature_texts)
    for extension in ('png', 'pdf'):
        fig.savefig(Path(OUT_DIR) / f'triplet_ABC_stacked.{extension}', dpi=300)
    plt.close(fig)
    inputs = [Path(OUT_DIR) / f'triplet_{cls}_synergistic_count_by_stage.tsv' for cls in 'ABC']
    (Path(OUT_DIR) / 'triplet_ABC_stacked_layout.json').write_text(json.dumps(dict(
        figure_size_cm=[width_cm, height_cm], original_pdf_height_preserved=False,
        original_pdf_size_cm=[36.01885364583333, 33.86539843055556],
        top_margin_cm=top_margin_cm, bottom_margin_cm=bottom_margin_cm, panel_gap_cm=panel_gap_cm,
        label_font_size_pt=8, cell_number_font_size_pt=6,
        cell_number_format='integer if value >= 10, otherwise one decimal; unchanged',
        source_summary_sha256={path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs},
        stages=stages, region_order=ORDERED_REGIONS, color_limits=color_limits,
        matrix_values={cls: np.asarray(matrices[cls]).tolist() for cls in 'ABC'},
        text_bounds_checked=True, cell_numbers_fit_cells=True,
        region_and_feature_labels_do_not_overlap=True, stage_labels_do_not_overlap=True,
        selection='observed_hoi < 0 and q_down_sample_block_bh <= 0.05',
        aggregation='Unchanged mean per-sample discovery count across available tissues in each stage'),
        indent=2) + '\n')
    print(f'[triplet_ABC_stacked] saved at {width_cm:g} x {height_cm:.6f} cm')
    if not write_summary:
        return
    coverage = [{'stage': stage, 'n_samples': sum(v == stage for v in base_of.values()),
                 'samples': ','.join(sorted(k for k, v in base_of.items() if v == stage))} for stage in stages]
    pd.DataFrame(coverage).to_csv(Path(OUT_DIR) / 'stage_coverage.tsv', sep='\t', index=False)
    (Path(OUT_DIR) / 'triplet_counts_plot_metadata.json').write_text(json.dumps(dict(
        source=str(Path(OUT_DIR) / 'all_triplets_fdr_10000.tsv.gz'),
        analysis_created_at_utc=meta['created_at_utc'], n_samples=meta['n_samples'],
        selection='observed_hoi < 0 and q_down_sample_block_bh <= 0.05',
        family='35 lower-tail p-values per sample x block',
        aggregation='Mean per-sample discovery count across available tissues in each stage'), indent=2) + '\n')


def redraw_saved_stacked():
    """Restyle only the stacked figure, using the existing stage summary tables."""
    root = Path(OUT_DIR)
    coverage = pd.read_csv(root / 'stage_coverage.tsv', sep='\t')
    assert not coverage.stage.duplicated().any()
    stages = [stage for stage in STAGE_ORDER if stage in set(coverage.stage)]
    assert len(stages) == len(coverage)
    base_of = {}
    for row in coverage.itertuples():
        names = row.samples.split(',')
        assert len(names) == row.n_samples
        for sample in names:
            assert sample not in base_of
            base_of[sample] = row.stage
    inputs = [root / f'triplet_{cls}_synergistic_count_by_stage.tsv' for cls in 'ABC']
    inputs.append(root / 'stage_coverage.tsv')
    before = {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs}
    stage_sizes = coverage.set_index('stage').n_samples
    matrices = {}
    for cls, path in zip('ABC', inputs):
        data = pd.read_csv(path, sep='\t')
        assert len(data) == len(stages)*len(ORDERED_REGIONS)
        assert not data.duplicated(['stage', 'region']).any()
        assert data.n_samples.eq(data.stage.map(stage_sizes)).all()
        matrix = data.pivot(index='stage', columns='region',
                            values='mean_fdr_significant_negative_count').reindex(
                                index=stages, columns=ORDERED_REGIONS).to_numpy()
        assert np.isfinite(matrix).all()
        matrices[cls] = matrix
    plot_stacked(stages, matrices, base_of, write_summary=False)
    assert all(hashlib.sha256(path.read_bytes()).hexdigest() == digest for path, digest in before.items())
    print('Saved stage summaries and sample coverage are unchanged.')


def main():
    counts = load_counts()
    samples = sorted(set(k[0] for k in counts))
    print(f'{len(samples)} samples, {len(counts)} (sample, region, class) cells')

    # base stage = first token of the sample name (E12.5_Heart -> E12.5)
    base_of = {s: s.split('_')[0] for s in samples}
    base_stages = [s for s in STAGE_ORDER if s in set(base_of.values())]
    print(f'{len(base_stages)} base stages: {base_stages}')

    classes = [
        ('A', 'DNA methylation + (H3K9me3 or H3K27me3 or H3K36me3) + Others', pattern_a),
        ('B', '(H3K9me3 or H3K27me3 or H3K36me3) + Others (no DNA methylation)', pattern_b),
        ('C', 'Not A and Not B', pattern_c),
    ]

    stage_matrices = {}
    for label, description, _ in classes:
        # per-sample matrix
        mat = np.full((len(samples), len(ORDERED_REGIONS)), np.nan)
        for i, s in enumerate(samples):
            for j, reg in enumerate(ORDERED_REGIONS):
                v = counts.get((s, reg, label))
                if v is not None:
                    mat[i, j] = v[0]

        # stage-averaged matrix
        avg_mat = np.full((len(base_stages), len(ORDERED_REGIONS)), np.nan)
        summary_rows = []
        for bi, bs in enumerate(base_stages):
            stage_samples = [s for s in samples if base_of[s] == bs]
            for j, reg in enumerate(ORDERED_REGIONS):
                vals = [counts[(s, reg, label)][0] for s in stage_samples
                        if (s, reg, label) in counts]
                sign_vals = [counts[(s, reg, label)][1] for s in stage_samples
                             if (s, reg, label) in counts]
                if vals:
                    avg_mat[bi, j] = float(np.mean(vals))
                    summary_rows.append((bs, reg, len(stage_samples),
                                         float(np.mean(vals)),
                                         float(np.mean(sign_vals)) if sign_vals else np.nan))

        stage_matrices[label] = avg_mat
        tsv = os.path.join(OUT_DIR, f'triplet_{label}_synergistic_count_by_stage.tsv')
        with open(tsv, 'w', newline='') as f:
            w = csv.writer(f, delimiter='\t', lineterminator='\n')
            w.writerow(['stage', 'region', 'n_samples',
                        'mean_fdr_significant_negative_count',
                        'mean_signbased_synergistic_count'])
            w.writerows(summary_rows)
        print(f'[Triplet {label}] summary: {tsv}')

        plot_heatmap(
            base_stages, avg_mat,
            f'Triplet {label}: {description}\n'
            f'Mean number of BH35-significant negative triplets per sample '
            f'(q_down <= {P_CUT}, O-info < 0), averaged across available tissues',
            r'Mean number of synergistic triplets',
            f'triplet_{label}_corrected_synergy_count_by_stage')

    plot_stacked(base_stages, stage_matrices, base_of)

    # ---- total heatmap: class-agnostic synergistic triplets ----
    # cell (stage, region) = mean over stage samples of the number of
    # permutation-corrected synergistic triplets of ANY class (A+B+C)
    total_mat = np.full((len(base_stages), len(ORDERED_REGIONS)), np.nan)
    summary_rows = []
    for bi, bs in enumerate(base_stages):
        stage_samples = [s for s in samples if base_of[s] == bs]
        for j, reg in enumerate(ORDERED_REGIONS):
            vals = [sum(counts[(s, reg, cls)][0] for cls in ('A', 'B', 'C')
                        if (s, reg, cls) in counts)
                    for s in stage_samples]
            if vals:
                total_mat[bi, j] = float(np.mean(vals))
                summary_rows.append((bs, reg, len(stage_samples), total_mat[bi, j]))

    tsv = os.path.join(OUT_DIR, 'total_synergistic_count_by_stage_region.tsv')
    with open(tsv, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(['stage', 'region', 'n_samples', 'mean_total_synergistic_count'])
        w.writerows(summary_rows)
    print(f'[total heatmap] summary: {tsv}')

    plot_heatmap(
        base_stages, total_mat,
        f'Total (A+B+C): all BH35-significant negative triplets\n'
        f'Mean number per sample (q_down <= {P_CUT}, O-info < 0), averaged across available tissues',
        r'Mean number of synergistic triplets',
        'total_corrected_synergy_count_by_stage')

    # ---- total synergistic triplet count along development ----
    per_sample_total = {}       # sample -> total over regions
    per_sample_class_total = {} # (sample, cls) -> total over regions
    for s in samples:
        per_sample_total[s] = sum(
            counts[(s, reg, cls)][0]
            for reg in ORDERED_REGIONS
            for cls in ('A', 'B', 'C')
            if (s, reg, cls) in counts)
        for cls in ('A', 'B', 'C'):
            per_sample_class_total[(s, cls)] = sum(
                counts[(s, reg, cls)][0]
                for reg in ORDERED_REGIONS
                if (s, reg, cls) in counts)

    totals = {}
    class_totals = {cls: {} for cls in ('A', 'B', 'C')}
    summary_rows = []
    for bs in base_stages:
        stage_samples = [s for s in samples if base_of[s] == bs]
        totals[bs] = float(np.mean([per_sample_total[s] for s in stage_samples]))
        row = [bs, len(stage_samples), totals[bs]]
        for cls in ('A', 'B', 'C'):
            class_totals[cls][bs] = float(np.mean(
                [per_sample_class_total[(s, cls)] for s in stage_samples]))
            row.append(class_totals[cls][bs])
        summary_rows.append(row)

    tsv = os.path.join(OUT_DIR, 'total_synergistic_triplet_count_by_stage.tsv')
    with open(tsv, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(['stage', 'n_samples', 'mean_total',
                    'mean_triplet_A', 'mean_triplet_B', 'mean_triplet_C'])
        w.writerows(summary_rows)
    print(f'[total] summary: {tsv}')

    plot_total_by_stage(base_stages, totals, class_totals,
                        'total_synergistic_triplet_count_by_stage')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stacked-only', action='store_true',
                        help='Restyle only triplet_ABC_stacked from saved stage tables, without recounting or redrawing other figures')
    args = parser.parse_args()
    if args.stacked_only:
        redraw_saved_stacked()
    else:
        main()
