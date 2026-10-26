#!/usr/bin/env python3
"""Three stacked A/B/C count heatmaps for the eight-signal E11.5–P0 panel."""
import csv
import gzip
import hashlib
import itertools
import json
from collections import Counter, defaultdict
from statistics import mean

from plot_triplet_A_B_C_post_E11p5 import (
    ROOT, SOURCE, STAGES, MARKS, REPRESSIVE, classify, write_tsv, plt, np,
)
from matplotlib.font_manager import findfont
from matplotlib.text import Text

REGIONS = [
    'promoter', 'proximal_ctcf_bound_enhancer_purified_bins',
    'proximal_not_ctcf_bound_enhancer_purified_bins',
    'distal_ctcf_bound_enhancer_purified_bins',
    'distal_not_ctcf_bound_enhancer_purified_bins',
    'exon_purified_bins', 'intron_purified_bins', 'None',
]
REGION_LABELS = ['promoter', 'proximal\nCTCF+', 'proximal\nCTCF−',
                 'distal\nCTCF+', 'distal\nCTCF−', 'exon', 'intron', 'None']
FEATURES = ['CpG_island', 'LINE', 'SINE', 'LTR', 'None']
BLOCKS = [r + '_AND_' + f for r in REGIONS for f in FEATURES]
PREFIX = 'triplet_ABC_stacked_post_E11p5_n5'


def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def draw_heatmaps(lookup, *, stages, stage_labels, prefix, period_title,
                  aggregation_title, colorbar_label='Mean count'):
    """Use the same fixed-size layout for the early and late eight-signal panels."""
    assert len(stages) == len(stage_labels) == 7
    plt.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'axes.titlesize': 8,
                         'axes.labelsize': 8, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
                         'figure.titlesize': 8, 'pdf.fonttype': 42, 'ps.fonttype': 42})
    font_path = findfont('Arial', fallback_to_default=False)
    width_cm, height_cm = 18, 12.3
    left_cm, bottom_cm, panel_width_cm, panel_height_cm, row_gap_cm = 1.85, 2.70, 14.30, 2.20, .65
    fig = plt.figure(figsize=(width_cm / 2.54, height_cm / 2.54))
    axes = [fig.add_axes([left_cm / width_cm,
                          (bottom_cm + (2 - i) * (panel_height_cm + row_gap_cm)) / height_cm,
                          panel_width_cm / width_cm, panel_height_cm / height_cm]) for i in range(3)]
    fig.suptitle(f'Significant synergistic triplets: {period_title}; eight shared signals\n'
                 f'10,000 permutations; O-information < 0; BH56 q ≤ 0.05; {aggregation_title}\n'
                 'Repressive group: H3K9me3, H3K27me3, H3K36me3', fontsize=8, y=1 - .12 / height_cm)
    titles = {'A': 'Triplet A: DNA methylation + at least one repressive mark (15 possible)',
              'B': 'Triplet B: at least one repressive mark; no DNA methylation (31 possible)',
              'C': 'Triplet C: no repressive mark (10 possible)'}
    cell_texts = []
    for ax, cls in zip(axes, 'ABC'):
        matrix = np.array([[lookup[cls, s, b] for b in BLOCKS] for s in stages])
        vmax = max(1., float(matrix.max()))
        im = ax.imshow(matrix, aspect='auto', cmap='Reds', vmin=0, vmax=vmax, interpolation='nearest')
        ax.set_yticks(range(7), stage_labels, fontsize=8)
        ax.set_xticks([])
        ax.set_title(titles[cls], fontsize=8, pad=5)
        for i in range(7):
            for j in range(40):
                value = matrix[i, j]
                cell_texts.append(ax.text(j, i, f'{value:.1f}', ha='center', va='center', fontsize=4,
                                          color='white' if value > .55 * vmax else 'black'))
        for j in range(5, 40, 5):
            ax.axvline(j - .5, color='black', linewidth=.8, alpha=.4)
        box = ax.get_position()
        cax = fig.add_axes([(left_cm + panel_width_cm + .15) / width_cm,
                            box.y0, .22 / width_cm, box.height])
        cb = fig.colorbar(im, cax=cax)
        cb.ax.tick_params(labelsize=8, length=2, pad=2)
        cb.set_label(colorbar_label, fontsize=8)
    ax = axes[-1]
    for i, label in enumerate(REGION_LABELS):
        ax.annotate('', xy=(5*i+4.3, -.28 / panel_height_cm), xytext=(5*i-.3, -.28 / panel_height_cm),
                    xycoords=ax.get_xaxis_transform(), arrowprops=dict(arrowstyle='<->', color='gray', lw=1))
        ax.text(5*i+2, -.40 / panel_height_cm, label,
                transform=ax.get_xaxis_transform(), ha='center', va='top', fontsize=8)
    for j, block in enumerate(BLOCKS):
        ax.text(j, -1.10 / panel_height_cm, block.split('_AND_')[1].replace('_', ' '),
                transform=ax.get_xaxis_transform(), ha='center', va='top', rotation=90, fontsize=8)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    cells = set(cell_texts)
    for text in fig.findobj(Text):
        if text.get_visible() and text.get_text():
            assert text.get_fontsize() == (4 if text in cells else 8), text.get_text()
            bbox = text.get_window_extent(renderer)
            assert bbox.x0 >= 0 and bbox.y0 >= 0 and bbox.x1 <= fig.bbox.width and bbox.y1 <= fig.bbox.height, text.get_text()
            if text in cells:
                panel_box = text.axes.get_window_extent(renderer)
                assert bbox.width < panel_box.width / 40 and bbox.height < panel_box.height / 7, text.get_text()
    for panel in axes:
        ticks = panel.get_yticklabels()
        for first, second in zip(ticks, ticks[1:]):
            assert not first.get_window_extent(renderer).overlaps(second.get_window_extent(renderer))
    for ext in ('png', 'pdf'):
        fig.savefig(ROOT / 'plots' / (prefix + '.' + ext), dpi=300)
    plt.close(fig)
    return titles, dict(figure_size_cm=[width_cm, height_cm], font_family='Arial', font_file=font_path,
                        label_font_size_pt=8, cell_number_font_size_pt=4, text_bounds_checked=True,
                        cell_numbers_fit_cells=True, stage_label_overlap_checked=True)


def main():
    tissues = json.loads((ROOT / 'inputs/source_manifest.json').read_text())['late_reference_tissues']
    universe = {frozenset(t) for t in itertools.combinations(MARKS, 3)}
    maxima = dict(Counter(classify(t) for t in universe))
    assert maxima == {'A': 15, 'B': 31, 'C': 10}
    # Explicitly guard the requested treatment of H3K36me3.
    assert classify(['DNA_methylation', 'H3K36me3', 'ATAC']) == 'A'
    assert classify(['H3K36me3', 'H3K4me1', 'ATAC']) == 'B'
    families = defaultdict(list)
    with gzip.open(SOURCE, 'rt') as handle:
        for r in csv.DictReader(handle, delimiter='\t'):
            assert r['stage'] in STAGES and r['tissue'] in tissues
            assert r['sample'] == r['stage'] + '_' + r['tissue']
            assert r['block'] in BLOCKS and int(r['n_perm']) == 10000
            assert frozenset(r['triplet'].split('/')) in universe
            families[r['sample'], r['block']].append(r)
    assert len(families) == 35 * 40
    per_sample, counts = [], {}
    for stage in STAGES:
        for tissue in tissues:
            sample = stage + '_' + tissue
            for block in BLOCKS:
                rows = families[sample, block]
                assert len(rows) == 56 and {frozenset(r['triplet'].split('/')) for r in rows} == universe
                ordered = sorted(rows, key=lambda r: float(r['p_down']))
                q = np.minimum(1, np.minimum.accumulate((np.array([float(r['p_down']) for r in ordered])
                                                        * 56 / np.arange(1, 57))[::-1])[::-1])
                assert np.allclose(q, [float(r['q_down_block_bh56']) for r in ordered], rtol=0, atol=1e-12)
                significant = [r for r in rows if float(r['observed_hoi']) < 0 and float(r['q_down_block_bh56']) <= .05]
                by_class = Counter(classify(r['triplet'].split('/')) for r in significant)
                assert sum(by_class.values()) == len(significant)
                for cls in 'ABC':
                    count = by_class[cls]
                    assert 0 <= count <= maxima[cls]
                    counts[sample, block, cls] = count
                    per_sample.append(dict(sample=sample, stage=stage, tissue=tissue, block=block,
                                           triplet_class=cls, count=count, theoretical_max=maxima[cls]))
    summary, lookup = [], {}
    for cls in 'ABC':
        for stage in STAGES:
            for block in BLOCKS:
                value = mean(counts[stage + '_' + t, block, cls] for t in tissues)
                record = dict(triplet_class=cls, stage=stage, block=block, mean_count=value,
                              theoretical_max=maxima[cls], n=5)
                summary.append(record)
                lookup[cls, stage, block] = value
    # Confirm agreement with all regions in the retained final trajectory figure.
    trajectories = read_tsv(ROOT / 'data/triplet_ABC_selected_region_trajectories_8signals_post_E11p5_n5_summary.tsv')
    assert len(trajectories) == len(STAGES) * 3 * len({r['block'] for r in trajectories})
    for r in trajectories:
        assert abs(lookup[r['triplet_class'], r['stage'], r['block']] - float(r['mean_count'])) < 1e-12
    # The sum of all mark occurrences must be three times the total triplet count.
    occurrences = read_tsv(ROOT / 'data/original_development_mark_participation_by_sample.tsv')
    for stage in STAGES:
        for tissue in tissues:
            sample = stage + '_' + tissue
            ref = [r for r in occurrences if r['sample'] == sample]
            assert len(ref) == 8
            assert sum(int(r['n_significant']) for r in ref) == 3 * sum(
                counts[sample, block, cls] for block in BLOCKS for cls in 'ABC')
    write_tsv(ROOT / 'data' / (PREFIX + '_by_sample.tsv'), per_sample)
    write_tsv(ROOT / 'data' / (PREFIX + '_summary.tsv'), summary)
    titles, figure_metadata = draw_heatmaps(
        lookup, stages=STAGES, stage_labels=[s + ' (n=5)' for s in STAGES], prefix=PREFIX,
        period_title='E11.5–P0', aggregation_title='mean count across five common tissues')
    validation = dict(status='complete', source=str(SOURCE), source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
                      signals=MARKS, repressive_marks=sorted(REPRESSIVE), theoretical_max=maxima,
                      class_definitions=titles, stages=STAGES, tissues=tissues, n_samples=35,
                      n_blocks=40, n_permutations=10000, BH_family_size=56,
                      criterion='observed_hoi < 0 and BH56 q <= 0.05', aggregation='Equal mean of five tissue counts per stage',
                      color_scale='Independent Reds scale for each class, matching reference',
                      **figure_metadata,
                      n_summary_rows=len(summary), n_sample_rows=len(per_sample),
                      eight_signal_families_complete=True, BH56_recomputed_and_matched=True,
                      ABC_partition_checked=True, H3K36me3_classification_checked=True,
                      trajectory_counts_match=True, mark_occurrence_totals_match=True)
    (ROOT / 'data' / (PREFIX + '_validation.json')).write_text(json.dumps(validation, indent=2) + '\n')
    print(json.dumps(validation, indent=2))


if __name__ == '__main__':
    main()
