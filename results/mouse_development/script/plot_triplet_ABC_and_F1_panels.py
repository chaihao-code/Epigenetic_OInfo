#!/usr/bin/env python3
"""Panel a: ABC counts; panel b: the saved H3K9me3/H3K27me3/X F1 heatmap."""
import hashlib
import json
from pathlib import Path

from plot_triplet_counts_corrected_by_stage import (
    OUT_DIR, STAGE_ORDER, ORDERED_REGIONS, plt, np, pd,
)
from matplotlib.text import Text

OUTPUT = Path(OUT_DIR)
PREFIX = 'triplet_ABC_and_F1_panels_ab'
TITLES = {
    'A': 'Triplet A: DNA methylation + (H3K9me3 or H3K27me3) + other',
    'B': 'Triplet B: (H3K9me3 or H3K27me3) + other (no DNA methylation)',
    'C': 'Triplet C: neither A nor B',
    'F1': 'H3K9me3 + H3K27me3 + X: F1 (X = ATAC, H3K27ac, H3K4me1, H3K4me3)',
}
REGION_LABELS = ['Promoter', 'Proximal\nenh. CTCF+', 'Proximal\nenh. CTCF−',
                 'Distal\nenh. CTCF+', 'Distal\nenh. CTCF−', 'Exon', 'Intron', 'None']


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_saved_values():
    coverage_path = OUTPUT / 'stage_coverage.tsv'
    coverage = pd.read_csv(coverage_path, sep='\t')
    assert not coverage.stage.duplicated().any()
    stages = [s for s in STAGE_ORDER if s in set(coverage.stage)]
    assert len(stages) == len(coverage) == 17
    sizes = coverage.set_index('stage').n_samples.reindex(stages)
    assert sizes.sum() == 82
    matrices, sources = {}, [coverage_path]
    for cls in 'ABC':
        path = OUTPUT / f'triplet_{cls}_synergistic_count_by_stage.tsv'
        sources.append(path)
        data = pd.read_csv(path, sep='\t')
        assert len(data) == 17*40 and not data.duplicated(['stage','region']).any()
        assert set(data.region) == set(ORDERED_REGIONS)
        assert data.n_samples.eq(data.stage.map(sizes)).all()
        matrices[cls] = data.pivot(index='stage', columns='region',
                                  values='mean_fdr_significant_negative_count').reindex(
                                      index=stages, columns=ORDERED_REGIONS).to_numpy()
    path = OUTPUT / 'f1_by_stage_block.tsv'
    sources.append(path)
    data = pd.read_csv(path, sep='\t')
    assert len(data) == 17*40 and not data.duplicated(['stage','block']).any()
    assert set(data.block) == set(ORDERED_REGIONS)
    assert data.n_samples.eq(data.stage.map(sizes)).all()
    matrices['F1'] = data.pivot(index='stage', columns='block', values='mean_F1').reindex(
        index=stages, columns=ORDERED_REGIONS).to_numpy()
    assert ((matrices['F1'] >= 0) & (matrices['F1'] <= 1)).all()
    # Check the saved stage means against the original per-sample F1 table.
    source = OUTPUT / 'f1_by_sample_block.tsv'
    sources.append(source)
    samples = pd.read_csv(source, sep='\t')
    assert len(samples) == 82*40 and not samples.duplicated(['sample','block']).any()
    assert (samples.TP + samples.FP).eq(4).all()
    expected = 2*samples.TP/(2*samples.TP+samples.FP+samples.FN)
    assert np.allclose(expected, samples.F1, rtol=0, atol=1e-14)
    means = samples.groupby(['stage','block']).F1.mean().unstack().reindex(
        index=stages, columns=ORDERED_REGIONS).to_numpy()
    assert np.allclose(means, matrices['F1'], rtol=0, atol=1e-14)
    for values in matrices.values():
        assert values.shape == (17,40) and np.isfinite(values).all()
    return stages, sizes, matrices, sources


def draw_panels(stages, sizes, matrices):
    plt.rcParams.update({'font.size':8, 'axes.titlesize':8, 'axes.labelsize':8,
                         'xtick.labelsize':8, 'ytick.labelsize':8, 'legend.fontsize':8})
    width_cm, height_cm = 18., 22.5
    left_cm, matrix_width_cm = 2.20, 14.5
    top_margin_cm, bottom_margin_cm = .35, 2.40
    gaps_cm = [.32, .32, 1.00]
    panel_height_cm = (height_cm-top_margin_cm-bottom_margin_cm-sum(gaps_cm))/4
    fig = plt.figure(figsize=(width_cm/2.54,height_cm/2.54), dpi=300)
    row_labels = list(stages)
    axes, title_texts, matrices_meta = [], [], {}
    top_cm = height_cm-top_margin_cm
    for index, key in enumerate(['A','B','C','F1']):
        bottom_cm = top_cm-panel_height_cm
        ax = fig.add_axes([left_cm/width_cm,bottom_cm/height_cm,
                           matrix_width_cm/width_cm,panel_height_cm/height_cm],
                          sharex=axes[0] if axes else None)
        axes.append(ax)
        values = matrices[key]
        vmax = 1. if key == 'F1' else max(1.,float(values.max()))
        matrices_meta[key] = {'values':values.tolist(), 'color_limits':[0,vmax]}
        image = ax.imshow(values,cmap='Reds',vmin=0,vmax=vmax,aspect='auto',interpolation='nearest')
        ax.set_xlim(-.5,39.5)
        ax.set_yticks(range(17),row_labels)
        ax.tick_params(axis='y',labelsize=8,length=2,pad=2,width=.6)
        ax.set_xticks([])
        title_texts.append(fig.text((left_cm+matrix_width_cm/2)/width_cm,
                                   (top_cm+.035)/height_cm,TITLES[key],
                                   fontsize=8,ha='center',va='bottom'))
        for boundary in range(5,40,5):
            ax.axvline(boundary-.5,color='black',lw=.8,alpha=.4)
        cax = fig.add_axes([16.84/width_cm,bottom_cm/height_cm,.15/width_cm,panel_height_cm/height_cm])
        cb = fig.colorbar(image,cax=cax)
        cb.ax.tick_params(labelsize=8,length=2,pad=2)
        cb.set_label('Mean F1' if key == 'F1' else 'Mean count',fontsize=8,labelpad=3)
        if key == 'F1':
            cb.set_ticks([0,.2,.4,.6,.8,1])
        else:
            cb.set_ticks([tick for tick in cb.get_ticks() if 0 <= tick <= vmax])
        if index in [0,3]:
            fig.text(.18/width_cm,(top_cm+.035)/height_cm,'a' if index == 0 else 'b',
                     ha='left',va='bottom',fontsize=8,fontweight='bold')
        if index < 3:
            top_cm = bottom_cm-gaps_cm[index]
    ax = axes[-1]
    region_texts,feature_texts = [],[]
    for index,label in enumerate(REGION_LABELS):
        arrow_y = -.12/panel_height_cm
        ax.annotate('',xy=(5*index+4.3,arrow_y),xytext=(5*index-.3,arrow_y),
                    xycoords=ax.get_xaxis_transform(),
                    arrowprops=dict(arrowstyle='<->',color='gray',lw=1.0),annotation_clip=False)
        region_texts.append(ax.text(5*index+2,-.22/panel_height_cm,label,
            transform=ax.get_xaxis_transform(),ha='center',va='top',fontsize=8,clip_on=False))
    for index,block in enumerate(ORDERED_REGIONS):
        feature_texts.append(ax.text(index,-.88/panel_height_cm,
            block.split('_AND_')[1].replace('_',' '),transform=ax.get_xaxis_transform(),
            ha='center',va='top',rotation=90,fontsize=8,clip_on=False))
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for text in fig.findobj(Text):
        if text.get_visible() and text.get_text():
            assert text.get_fontsize() == 8,text.get_text()
            box = text.get_window_extent(renderer)
            assert box.x0 >= 0 and box.y0 >= 0 and box.x1 <= fig.bbox.width and box.y1 <= fig.bbox.height,text.get_text()
    for ax,title in zip(axes,title_texts):
        labels = ax.get_yticklabels()
        for first,second in zip(labels,labels[1:]):
            assert not first.get_window_extent(renderer).overlaps(second.get_window_extent(renderer)), 'stage labels overlap'
        assert not title.get_window_extent(renderer).overlaps(ax.get_window_extent(renderer)),title.get_text()
    for previous,title in zip(axes,title_texts[1:]):
        assert not previous.get_window_extent(renderer).overlaps(title.get_window_extent(renderer)),title.get_text()
    for first,second in zip(region_texts,region_texts[1:]):
        assert not first.get_window_extent(renderer).overlaps(second.get_window_extent(renderer))
    assert min(t.get_window_extent(renderer).y0 for t in region_texts) > max(
        t.get_window_extent(renderer).y1 for t in feature_texts)
    reference_x = axes[0].transData.transform(np.column_stack([np.arange(40),np.zeros(40)]))[:,0]
    for ax in axes[1:]:
        assert np.allclose(reference_x,ax.transData.transform(np.column_stack([np.arange(40),np.zeros(40)]))[:,0],rtol=0,atol=1e-8)
    for extension in ['pdf','png']:
        output = OUTPUT/(PREFIX+'.'+extension)
        fig.savefig(output,dpi=300)
        print(output,flush=True)
    plt.close(fig)
    return dict(figure_size_cm=[width_cm,height_cm],panels={'a':['A','B','C'],'b':['F1']},
                titles=TITLES,label_font_size_pt=8,cell_annotations_shown=False,
                sample_sizes_shown=False,stage_labels=row_labels,region_labels=REGION_LABELS,
                method_caption_shown=False,
                panel_height_cm=panel_height_cm,panel_gaps_cm=gaps_cm,
                matrices=matrices_meta,shared_x_labels_shown_once=True,
                shared_x_positions_verified=True,text_bounds_checked=True,
                stage_labels_do_not_overlap=True,
                region_and_feature_labels_do_not_overlap=True)


def main():
    stages,sizes,matrices,sources = load_saved_values()
    retained = list(OUTPUT.glob('*.tsv'))+[OUTPUT/'all_triplets_fdr_10000.tsv.gz']
    retained += [OUTPUT/(stem+'.'+extension) for stem in ['triplet_ABC_stacked','H3K9me3_H3K27me3_X_F1_heatmap']
                 for extension in ['pdf','png']]
    before = {path:sha256(path) for path in retained}
    metadata = draw_panels(stages,sizes,matrices)
    assert all(sha256(path) == digest for path,digest in before.items())
    metadata.update(source_sha256={path.name:sha256(path) for path in sources},
                    original_data_and_standalone_figures_unchanged=True,
                    stages=stages,regions=ORDERED_REGIONS,sample_counts=sizes.to_dict(),
                    ABC_selection='observed_hoi < 0 and q_down_sample_block_bh <= .05; BH35',
                    ABC_aggregation='Original per-stage equal-weight mean across available tissues',
                    F1_universe='Triplet B; no DNA methylation and at least one of H3K9me3 or H3K27me3',
                    F1_prediction='H3K9me3 + H3K27me3 + X; X = ATAC,H3K27ac,H3K4me1,H3K4me3',
                    F1_aggregation='Original arithmetic mean of per-sample F1 within stage, not pooled F1',
                    F1_saved_values_independently_checked=True)
    (OUTPUT/(PREFIX+'_layout.json')).write_text(json.dumps(metadata,indent=2)+'\n')
    print('Verified saved values, 40-column alignment, 22.5 cm height, fonts and text bounds.',flush=True)


if __name__ == '__main__':
    main()
