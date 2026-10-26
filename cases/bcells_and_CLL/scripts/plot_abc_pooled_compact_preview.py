#!/usr/bin/env python3
"""Preview a compact ABC layout without changing counts or color normalization."""
import os
import hashlib
import json

from plot_abc_pooled_combined import PLOTS, BLOCKS, CAPACITIES, plt, font_manager, Path, pd, np
from matplotlib.text import Text

PREFIX = 'ABC_pooled_average_heatmap_compact_preview'


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    source = PLOTS / 'ABC_pooled_average_counts.tsv'
    retained = [source, PLOTS / 'ABC_pooled_average_heatmap.pdf', PLOTS / 'ABC_pooled_average_heatmap.png']
    before = {p: sha256(p) for p in retained}
    data = pd.read_csv(source, sep='\t')
    assert len(data) == 360
    font = Path(os.environ.get('BCELLS_ARIAL_FONT', str(Path.home() / '.local/share/fonts/arial.ttf')))
    if font.exists():
        font_manager.fontManager.addfont(str(font))
    plt.rcParams.update({'font.family': 'Arial' if font.exists() else 'DejaVu Sans',
                         'font.size': 8, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
                         'axes.labelsize': 8, 'pdf.fonttype': 42, 'ps.fonttype': 42})
    width_cm, height_cm = 18., 3.
    left_cm, bottom_cm, heatmap_width_cm, heatmap_height_cm = 1.85, .30, 14.6, 2.40
    panel_height_cm = heatmap_height_cm / 3
    fig = plt.figure(figsize=(width_cm/2.54, height_cm/2.54), dpi=300)
    axes, cell_texts, colorbars = [], [], []
    matrices = {}
    for index, cls in enumerate('ABC'):
        panel_bottom_cm = bottom_cm + (2-index)*panel_height_cm
        sub = data[data['class'].eq(cls)]
        matrix = sub.pivot(index='pooled_group', columns='block', values='mean_count').reindex(
            index=['R', 'S', 'pB'], columns=BLOCKS).to_numpy()
        assert matrix.shape == (3, 40) and np.isfinite(matrix).all()
        assert ((matrix >= 0) & (matrix <= CAPACITIES[cls])).all()
        matrices[cls] = matrix.tolist()
        ax = fig.add_axes([left_cm/width_cm, panel_bottom_cm/height_cm,
                           heatmap_width_cm/width_cm, panel_height_cm/height_cm])
        axes.append(ax)
        # Display the original mean counts using their original per-class limits.
        im = ax.imshow(matrix, cmap='YlOrRd', vmin=0, vmax=CAPACITIES[cls],
                       aspect='auto', interpolation='nearest')
        ax.set_yticks(range(3), ['R', 'S', 'pB'])
        ax.tick_params(axis='y', length=0, pad=2, labelsize=8)
        ax.set_xticks([])
        for k in range(1, 8):
            ax.axvline(k*5-.5, color='gray', lw=.5)
        for i in range(3):
            for j in range(40):
                cell_texts.append(ax.text(j, i, f'{matrix[i,j]:.1f}', ha='center', va='center',
                                          fontsize=3.5, color='white' if matrix[i,j] >= CAPACITIES[cls]*.62 else 'black'))
        fig.text(.67/width_cm, (panel_bottom_cm+panel_height_cm/2)/height_cm,
                 'Triplet '+cls, ha='center', va='center', fontsize=8)
        cax = fig.add_axes([16.65/width_cm, (panel_bottom_cm+.16)/height_cm,
                            .13/width_cm, .48/height_cm])
        cb = fig.colorbar(im, cax=cax, ticks=[0, CAPACITIES[cls]])
        cb.ax.tick_params(labelsize=8, length=1.5, pad=1.5)
        colorbars.append(cb)
    region_names = ['Promoter', 'pEnh+', 'pEnh−', 'dEnh+', 'dEnh−', 'Exon', 'Intron', 'None']
    for k, name in enumerate(region_names):
        x_cm = left_cm + (k*5+2.5)*heatmap_width_cm/40
        fig.text(x_cm/width_cm, .275/height_cm, name, ha='center', va='top', fontsize=8)
    feature_names = ['C', 'L', 'S', 'T', '–']
    for j in range(40):
        x_cm = left_cm + (j+.5)*heatmap_width_cm/40
        fig.text(x_cm/width_cm, (height_cm-.025)/height_cm, feature_names[j%5],
                 ha='center', va='top', fontsize=8)
    fig.text(17.70/width_cm, .5, 'Proportion', rotation=90, ha='center', va='center', fontsize=8)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    cells = set(cell_texts)
    for text in fig.findobj(Text):
        if text.get_visible() and text.get_text():
            assert text.get_fontsize() == (3.5 if text in cells else 8), text.get_text()
            box = text.get_window_extent(renderer)
            assert box.x0 >= 0 and box.y0 >= 0 and box.x1 <= fig.bbox.width and box.y1 <= fig.bbox.height, text.get_text()
            if text in cells:
                panel_box = text.axes.get_window_extent(renderer)
                assert box.width < panel_box.width/40 and box.height < panel_box.height/3
    for ax in axes:
        labels = ax.get_yticklabels()
        for first, second in zip(labels, labels[1:]):
            assert not first.get_window_extent(renderer).overlaps(second.get_window_extent(renderer))
    for ext in ['png', 'pdf']:
        path = PLOTS / (PREFIX+'.'+ext)
        fig.savefig(path, dpi=300)
        print(path)
    plt.close(fig)
    for path, old_digest in before.items():
        assert sha256(path) == old_digest
    metadata = dict(status='preview', figure_size_cm=[width_cm, height_cm], source=str(source),
                    source_sha256=before[source], authoritative_outputs_preserved=True,
                    original_mean_counts_unchanged=True, cell_annotation_font_size_pt=3.5,
                    other_font_sizes_pt=8, colorbar_display_label='Proportion',
                    color_normalization='Unchanged: original mean counts, no division by class capacity',
                    color_limits={c:[0, CAPACITIES[c]] for c in 'ABC'},
                    row_order={c:['R','S','pB'] for c in 'ABC'},
                    matrix_values=matrices, text_bounds_checked=True,
                    row_labels_do_not_overlap=True, cell_numbers_fit_cells=True,
                    abbreviated_axis_labels=True, full_axis_annotations_fit=False,
                    feature_key={'C':'CpG island','L':'LINE','S':'SINE','T':'LTR','–':'None'},
                    region_key={'pEnh+':'Proximal enhancer CTCF+', 'pEnh−':'Proximal enhancer CTCF−',
                                'dEnh+':'Distal enhancer CTCF+', 'dEnh−':'Distal enhancer CTCF−'},
                    note='Compact preview uses abbreviated annotations; full labels and FDR caption remain in the retained original figure.')
    (PLOTS/(PREFIX+'_layout.json')).write_text(json.dumps(metadata, indent=2)+'\n')
    print('Compact preview only: original count table, PDF and PNG are unchanged.')


if __name__ == '__main__':
    main()
