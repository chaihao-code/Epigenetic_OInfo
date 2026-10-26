#!/usr/bin/env python3
"""Restyle the original TC/DTC trajectories and stack TC above DTC."""

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
import csv
import hashlib
import json
import os
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

os.environ.setdefault('MPLCONFIGDIR', '/tmp/mouse-development-mpl')
os.environ.setdefault('XDG_CACHE_HOME', '/tmp/mouse-development-cache')
import plot_combined_hoi_trajectories_2x4 as original
import pandas as pd
import numpy as np
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.text import Text
from matplotlib.ticker import MaxNLocator, StrMethodFormatter
from matplotlib.patches import Patch

plt = original.plt
OUTPUT = Path(__file__).resolve().parent
PROJECT = OUTPUT.parent
WORKSPACE_PROJECTS = PROJECT.parent
SOURCE = WORKSPACE_PROJECTS / 'permutation_test/2026-09-28/results_mouse_development'
ARCHIVE = WORKSPACE_PROJECTS / 'mouse_development_fine_te/results_oinfo_n_bins_5'
WIDTH_CM = 18.
ORIGINAL_ASPECT = {'tc':.6115615887418628, 'dtc':.6200183311852661}
COMBINED_HEIGHT_CM = 20.
METRIC_GAP_CM = .5
METRIC_HEIGHTS = {key:(COMBINED_HEIGHT_CM-METRIC_GAP_CM)/2 for key in ['tc','dtc']}
EARLY_TIME_STEP = 1.25
PREFIX = 'combined_tc_dtc_trajectories_6x2_n5'
TITLES = ['Promoter / CpG island', 'Promoter / LINE',
          'Proximal enhancer CTCF- / CpG island', 'Proximal enhancer CTCF- / LINE',
          'Distal enhancer CTCF- / CpG island', 'Distal enhancer CTCF- / LINE']


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_tree():
    tree = original.nx.DiGraph()
    for line in original.TREE_PATH.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        parent,children = line.split('\t')
        for child in children.split(','):
            tree.add_edge(parent,child.strip())
    assert original.nx.is_arborescence(tree)
    assert original.nx.descendants(tree,'zygote')|{'zygote'} == set(tree)
    return tree


def load_original_values():
    """Read the frozen overall-system table retained with this release."""
    tree = read_tree()
    report = OUTPUT / 'tc_dtc_missing_inputs.csv'
    issues = pd.read_csv(report)
    assert issues.reason.eq('sample directory missing').all()
    samples = set(tree) - set(issues.cell_type)
    path = OUTPUT / 'tc_dtc_overall_values_for_plot.tsv'
    table = pd.read_csv(path, sep='\t', float_precision='round_trip')
    blocks = [item[2] for item in original.BLOCKS]
    assert len(samples) == 82 and set(table['sample']) == samples
    assert len(table) == 82 * 6 and not table.duplicated(['sample', 'block']).any()
    assert set(table.block) == set(blocks) and table.groupby('sample').size().eq(6).all()
    assert table.order.eq(7).all() and table.n_missing.eq(0).all()
    assert np.isfinite(table[['hoi', 'tc', 'dtc']].to_numpy()).all()
    assert np.allclose(table.tc - table.dtc, table.hoi, rtol=0, atol=1e-10)
    values = {(r.sample, r.block, metric): getattr(r, metric)
              for r in table.itertuples() for metric in ['tc', 'dtc']}
    info = dict(source_directory=str(OUTPUT), source_file_sha256={str(path):sha256(path)},
                n_samples=82, n_blocks=6, n_values=len(values), n_signals=7, n_bins=5,
                order=7, n_missing=0, no_recomputation=True, input_kind='retained overall-system intermediate TSV',
                differentiation_tree=str(original.TREE_PATH), tree_sha256=sha256(original.TREE_PATH),
                original_missing_inputs_report_sha256=sha256(report),
                excluded_missing_tree_samples=sorted(set(issues.cell_type)), original_sample_coverage_preserved=True)
    print(f'Loaded 82 samples x 6 blocks from retained intermediate table {path}', flush=True)
    return tree, values, table, info


def configure_fonts():
    path = Path(_release_path('@ROOT@/data/fonts/arial.ttf'))
    if path.exists():
        font_manager.fontManager.addfont(str(path))
    plt.rcParams.update({'font.family':'Arial' if path.exists() else 'DejaVu Sans',
                         'font.size':8,'axes.titlesize':8,'axes.labelsize':8,
                         'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8,
                         'pdf.fonttype':42,'ps.fonttype':42})


def line_signature(ax):
    digest = hashlib.sha256()
    for line in ax.lines:
        for coords in [line.get_xdata(),line.get_ydata()]:
            digest.update(np.asarray(coords,dtype=float).tobytes())
        digest.update(str((line.get_color(),line.get_linestyle(),line.get_linewidth(),line.get_marker())).encode())
    return digest.hexdigest()


def draw_metric(fig,metric,bottom_cm,height_cm,tree,values):
    figure_height = fig.get_figheight()*2.54
    left_positions = [1.,7.9]
    panel_width = 6.3
    top_margin,bottom_margin,row_gap = .86,1.34,.48
    panel_height = (height_cm-top_margin-bottom_margin-2*row_gap)/3
    heading = fig.text(7.6/WIDTH_CM,(bottom_cm+height_cm-.07)/figure_height,
                       metric.upper(),fontsize=10,fontweight='bold',ha='center',va='top')
    axes,titles,signatures = [],[],[]
    for index,(_,_,block) in enumerate(original.BLOCKS):
        row,col = divmod(index,2)
        top = bottom_cm+height_cm-top_margin-row*(panel_height+row_gap)
        ax = fig.add_axes([left_positions[col]/WIDTH_CM,(top-panel_height)/figure_height,
                           panel_width/WIDTH_CM,panel_height/figure_height])
        trajectories = original.collect_trajectories(block,metric,tree,values)
        assert trajectories
        # Widen categorical early-stage intervals to keep the 8 pt labels readable
        # after narrowing the axes; preserve days, values and the late daily scale.
        before = {(d['cell_type'],d['e_day'],d['value']) for points in trajectories.values() for d in points}
        anchor = min(d['x_pos'] for points in trajectories.values() for d in points if d['stage']=='organogenesis')
        for points in trajectories.values():
            for d in points:
                old_x = d['x_pos']
                d['x_pos'] = old_x*EARLY_TIME_STEP if old_x < anchor else anchor*EARLY_TIME_STEP+(old_x-anchor)
        assert before == {(d['cell_type'],d['e_day'],d['value']) for points in trajectories.values() for d in points}
        original.plot_on_ax(ax,trajectories)
        unchanged = line_signature(ax)
        ax.margins(x=.01)
        ax.tick_params(axis='both',pad=2,length=1.5,labelsize=8)
        locator = MaxNLocator(nbins=5,steps=[1,2,2.5,5,10])
        lo,hi = ax.get_ylim()
        ax.set_yticks([tick for tick in locator.tick_values(lo,hi) if lo <= tick <= hi])
        ax.yaxis.set_major_formatter(StrMethodFormatter('{x:g}'))
        if col == 0:
            ax.set_ylabel(metric.upper(),fontsize=8,labelpad=2)
        if row == 2:
            for text in ax.get_xticklabels():
                text.set_rotation(90)
                text.set_horizontalalignment('center')
                text.set_verticalalignment('top')
            ax.set_xlabel('Stage',fontsize=8,labelpad=2)
        else:
            ax.tick_params(axis='x',labelbottom=False)
        for text in ax.findobj(Text):
            text.set_fontsize(8)
        titles.append(fig.text((left_positions[col]+panel_width/2)/WIDTH_CM,
                               (top+.06)/figure_height,TITLES[index],fontsize=8,ha='center',va='bottom'))
        assert unchanged == line_signature(ax)
        axes.append(ax)
        signatures.append(unchanged)
    return axes,titles,heading,dict(metric=metric,height_cm=height_cm,panel_size_cm=[panel_width,panel_height],
                                  y_axis_limits=[list(ax.get_ylim()) for ax in axes],
                                  x_axis_limits=[list(ax.get_xlim()) for ax in axes],
                                  line_signature_sha256=signatures,
                                  display_coordinates_unchanged_during_styling=True,
                                  trajectory_days_values_colors_and_styles_preserved=True)


def add_legend(fig,center_y_cm):
    handles = [Patch(facecolor='white',edgecolor='none',label='Preimplantation'),
               Line2D([],[],color=original.preimpl_color,lw=1,alpha=.6,label='Preimpl.'),
               Patch(facecolor='white',edgecolor='none',label='Gastrulation')]
    for node in original.gastrulation_nodes:
        handles.append(Line2D([],[],color=original.gastrulation_color_map[node],lw=1,alpha=.6,label=node))
    handles.append(Patch(facecolor='white',edgecolor='none',label='Organogenesis'))
    for tissue in original.organ_tissues:
        label = 'Embryonic facial\nprominence' if tissue == 'Embryonic_Facial_Prominence' else tissue.replace('_',' ')
        handles.append(Line2D([],[],color=original.organ_color_map[tissue],lw=1,alpha=.6,label=label))
    legend = fig.legend(handles=handles,loc='center left',bbox_to_anchor=(14.55/WIDTH_CM,center_y_cm/(fig.get_figheight()*2.54)),
                        frameon=False,fontsize=8,handlelength=1.,handletextpad=.35,labelspacing=.25,
                        borderpad=0,borderaxespad=0)
    for text in legend.get_texts():
        if text.get_text() in ['Preimplantation','Gastrulation','Organogenesis']:
            text.set_fontweight('bold')
    return legend


def validate_and_save(fig,axes,titles,headings,legend,prefix):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    heading_set = set(headings)
    for text in fig.findobj(Text):
        if text.get_visible() and text.get_text():
            assert text.get_fontsize() == (10 if text in heading_set else 8),text.get_text()
            box = text.get_window_extent(renderer)
            assert box.x0 >= 0 and box.y0 >= 0 and box.x1 <= fig.bbox.width and box.y1 <= fig.bbox.height,text.get_text()
    box = legend.get_window_extent(renderer)
    assert box.x0 >= 0 and box.y0 >= 0 and box.x1 <= fig.bbox.width and box.y1 <= fig.bbox.height
    assert box.x1 <= fig.bbox.width-.4/2.54*fig.dpi, 'Legend must leave at least 0.4 cm at the right edge'
    for ax,title in zip(axes,titles):
        labels = [text for text in ax.get_xticklabels() if text.get_visible()]
        for first,second in zip(labels,labels[1:]):
            assert not first.get_window_extent(renderer).overlaps(second.get_window_extent(renderer)),('overlapping time labels',first.get_text(),second.get_text())
        assert not title.get_window_extent(renderer).overlaps(ax.get_window_extent(renderer))
    for heading,title in zip(headings,titles[::6]):
        assert not heading.get_window_extent(renderer).overlaps(title.get_window_extent(renderer))
    for ext in ['pdf','png']:
        path = OUTPUT/(prefix+'.'+ext)
        fig.savefig(path,dpi=300)
        print(path,flush=True)
    plt.close(fig)


def main():
    configure_fonts()
    tree,values,table,source_info = load_original_values()
    original_report_hash = sha256(OUTPUT/'tc_dtc_missing_inputs.csv')
    standalone = {}
    for metric in ['tc','dtc']:
        height = METRIC_HEIGHTS[metric]
        fig = plt.figure(figsize=(WIDTH_CM/2.54,height/2.54),dpi=300)
        axes,titles,heading,layout = draw_metric(fig,metric,0,height,tree,values)
        legend = add_legend(fig,height/2)
        validate_and_save(fig,axes,titles,[heading],legend,f'combined_{metric}_trajectories_3x2_n5')
        standalone[metric] = layout
    height = sum(METRIC_HEIGHTS.values())+METRIC_GAP_CM
    fig = plt.figure(figsize=(WIDTH_CM/2.54,height/2.54),dpi=300)
    top = draw_metric(fig,'tc',METRIC_HEIGHTS['dtc']+METRIC_GAP_CM,METRIC_HEIGHTS['tc'],tree,values)
    bottom = draw_metric(fig,'dtc',0,METRIC_HEIGHTS['dtc'],tree,values)
    axes,titles = top[0]+bottom[0],top[1]+bottom[1]
    legend = add_legend(fig,height/2)
    validate_and_save(fig,axes,titles,[top[2],bottom[2]],legend,PREFIX)
    assert original_report_hash == sha256(OUTPUT/'tc_dtc_missing_inputs.csv')
    table.to_csv(OUTPUT/'tc_dtc_overall_values_for_plot.tsv',sep='\t',index=False)
    metadata = dict(figure_size_cm=[WIDTH_CM,height],label_font_size_pt=8,stack_order=['TC','DTC'],
                    metric_heading_font_size_pt=10,legend_right_margin_at_least_cm=.4,
                    metric_gap_cm=METRIC_GAP_CM,original_metric_aspect_ratios=ORIGINAL_ASPECT,
                    standalone_layouts=standalone,combined_layouts={'tc':top[3],'dtc':bottom[3]},
                    shared_legend=True,region_titles=TITLES,x_tick_rotation_degrees=90,
                    x_coordinates='Categorical preimplantation/gastrulation: 1.25 units per time point; organogenesis: original 2 units per day; ICM and TE share E3.5',
                    early_categorical_time_step=EARLY_TIME_STEP,
                    time_spacing_adjusted_for_8pt_labels=True,
                    text_bounds_checked=True,time_labels_do_not_overlap=True,
                    original_missing_input_report_unchanged=True,**source_info)
    (OUTPUT/(PREFIX+'_validation.json')).write_text(json.dumps(metadata,indent=2)+'\n')
    print(f'Combined figure: {WIDTH_CM:g} x {height:.4f} cm; TC/DTC headings 10 pt, other text 8 pt.',flush=True)


if __name__ == '__main__':
    main()
