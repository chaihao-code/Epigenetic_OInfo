#!/usr/bin/env python3
"""Plot total seven-signal O-information in the TC/DTC layout, 18 x 8 cm."""
import hashlib
import json
from pathlib import Path

import plot_tc_dtc_stacked_18cm as shared

pd, np, plt = shared.pd, shared.np, shared.plt
original = shared.original
OUTPUT = Path(__file__).resolve().parent
WIDTH_CM = 18.
HEIGHT_CM = 8.
PREFIX = 'combined_oinfo_trajectories_3x2_n5'
CACHE = OUTPUT / 'tc_dtc_overall_values_for_plot.tsv'
TC_DTC_METADATA = OUTPUT / 'combined_tc_dtc_trajectories_6x2_n5_validation.json'


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_total_oinfo():
    """Reuse the exact overall-system table saved for the original TC/DTC plots."""
    table = pd.read_csv(CACHE, sep='\t', float_precision='round_trip')
    provenance = json.loads(TC_DTC_METADATA.read_text())
    tree = shared.read_tree()
    missing_report = OUTPUT / 'tc_dtc_missing_inputs.csv'
    missing = pd.read_csv(missing_report)
    assert missing.reason.eq('sample directory missing').all()
    expected_samples = set(tree) - set(missing.cell_type)
    expected_blocks = [block for _, _, block in original.BLOCKS]
    assert len(expected_samples) == provenance['n_samples'] == 82
    assert set(table['sample']) == expected_samples
    assert set(table['block']) == set(expected_blocks)
    assert len(table) == 82 * 6
    assert not table.duplicated(['sample', 'block']).any()
    assert table.groupby('sample').size().eq(6).all()
    assert table['order'].eq(7).all() and table.n_missing.eq(0).all()
    assert provenance['n_signals'] == 7 and provenance['n_bins'] == 5
    assert provenance['tree_sha256'] == sha256(original.TREE_PATH)
    assert provenance['original_missing_inputs_report_sha256'] == sha256(missing_report)
    assert np.isfinite(table[['hoi', 'tc', 'dtc']].to_numpy()).all()
    error = np.abs(table.hoi.to_numpy() - (table.tc - table.dtc).to_numpy())
    assert error.max() < 1e-10
    values = {(row.sample, row.block, 'oinfo'): row.hoi for row in table.itertuples()}
    coverage = {}
    for block in expected_blocks:
        trajectories = original.collect_trajectories(block, 'oinfo', tree, values)
        actual_samples = {point['cell_type'] for points in trajectories.values() for point in points}
        assert actual_samples == expected_samples, block
        coverage[block] = len(actual_samples)
    sources = {
        'plot_input_table': str(CACHE), 'plot_input_table_sha256': sha256(CACHE),
        'input_column': 'hoi', 'n_samples': 82, 'n_blocks': 6, 'n_values': len(values),
        'n_signals': 7, 'n_bins': 5, 'order': 7, 'n_missing': 0,
        'metric_definition': 'Total O-information = TC - DTC',
        'max_abs_error_hoi_vs_tc_minus_dtc': float(error.max()),
        'no_recomputation': True, 'per_block_trajectory_sample_coverage': coverage,
        'tc_dtc_validation_file': str(TC_DTC_METADATA),
        'tc_dtc_validation_file_sha256': sha256(TC_DTC_METADATA),
        'original_archived_source_directory': provenance['source_directory'],
        'original_archived_source_file_sha256': provenance['source_file_sha256'],
        'archived_regional_comparisons_in_tc_dtc_run': provenance['archived_region_comparisons'],
        'differentiation_tree': str(original.TREE_PATH), 'tree_sha256': sha256(original.TREE_PATH),
        'excluded_missing_tree_samples': sorted(set(missing.cell_type)),
        'missing_inputs_report_sha256': sha256(missing_report),
    }
    print(f'Loaded {len(values)} total O-information values (82 samples x 6 regions).', flush=True)
    return tree, values, table, sources


def compact_legend(fig):
    """Preserve all tissue colors while fitting the legend within an 8 cm page."""
    handles = [shared.Line2D([], [], color=original.preimpl_color, lw=1, alpha=.6,
                             label='Preimplantation'),
               shared.Patch(facecolor='white', edgecolor='none', label='Gastrulation')]
    for node in original.gastrulation_nodes:
        handles.append(shared.Line2D([], [], color=original.gastrulation_color_map[node],
                                     lw=1, alpha=.6, label=node))
    handles.append(shared.Patch(facecolor='white', edgecolor='none', label='Organogenesis'))
    for tissue in original.organ_tissues:
        label = ('Embryonic facial\nprominence' if tissue == 'Embryonic_Facial_Prominence'
                 else tissue.replace('_', ' '))
        handles.append(shared.Line2D([], [], color=original.organ_color_map[tissue],
                                     lw=1, alpha=.6, label=label))
    legend = fig.legend(handles=handles, loc='center left',
                        bbox_to_anchor=(14.55 / WIDTH_CM, .5), frameon=False,
                        fontsize=8, handlelength=1., handletextpad=.35,
                        labelspacing=.03, borderpad=0, borderaxespad=0)
    for text in legend.get_texts():
        if text.get_text() in ['Preimplantation', 'Gastrulation', 'Organogenesis']:
            text.set_fontweight('bold')
    return legend


def main():
    shared.configure_fonts()
    tree, values, table, sources = load_total_oinfo()
    fig = plt.figure(figsize=(WIDTH_CM / 2.54, HEIGHT_CM / 2.54), dpi=300)
    axes, titles, heading, layout = shared.draw_metric(
        fig, 'oinfo', 0, HEIGHT_CM, tree, values)
    heading.set_text('O-information')
    # Negative tick labels and the longer metric name need a little more room
    # than TC/DTC. Keep each panel's right edge and the legend position fixed.
    heading.set_x(heading.get_position()[0] + .1 / WIDTH_CM)
    for ax, title in zip(axes, titles):
        position = ax.get_position()
        ax.set_position([position.x0 + .2 / WIDTH_CM, position.y0,
                         position.width - .2 / WIDTH_CM, position.height])
        title.set_x(title.get_position()[0] + .1 / WIDTH_CM)
    layout['panel_size_cm'][0] = 6.1
    layout['extra_left_margin_for_negative_ticks_cm'] = .2
    for i, ax in enumerate(axes):
        if i % 2 == 0:
            ax.set_ylabel('O-information', fontsize=8, labelpad=2)
        assert shared.line_signature(ax) == layout['line_signature_sha256'][i]
    legend = compact_legend(fig)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for ax in axes:
        tick_boxes = [text.get_window_extent(renderer) for text in ax.get_yticklabels()
                      if text.get_visible() and text.get_text()]
        assert not any(a.overlaps(b) for i, a in enumerate(tick_boxes)
                       for b in tick_boxes[i + 1:]), 'Overlapping vertical tick labels'
        annotations = [t for t in ax.texts if t.get_text() in ['ICM', 'TE']]
        assert len(annotations) == 2
        assert not annotations[0].get_window_extent(renderer).overlaps(
            annotations[1].get_window_extent(renderer)), 'Overlapping ICM/TE labels'
    legend_box = legend.get_window_extent(renderer)
    assert legend_box.x0 > max(ax.bbox.x1 for ax in axes), 'Legend overlaps a panel'
    metadata = {
        'figure_size_cm': [WIDTH_CM, HEIGHT_CM],
        'reference_dtc_figure_size_cm': [18, shared.METRIC_HEIGHTS['dtc']],
        'height_ratio_vs_dtc': HEIGHT_CM / shared.METRIC_HEIGHTS['dtc'],
        'metric': 'total O-information', 'label_font_size_pt': 8,
        'metric_heading_font_size_pt': 10, 'layout': layout,
        'region_titles': shared.TITLES,
        'zero_reference_line': True,
        'tissue_colors_line_styles_and_time_coordinates_match_tc_dtc': True,
        'early_categorical_time_step': shared.EARLY_TIME_STEP,
        'x_coordinates': 'Same TC/DTC display scale: early categorical step 1.25; organogenesis 2 units/day; ICM and TE share E3.5',
        'legend_labelspacing': .03, 'legend_n_entries': len(legend.get_texts()),
        'legend_right_margin_at_least_cm': .4,
        'legend_bbox_cm': [v / fig.dpi * 2.54 for v in legend_box.bounds],
        'validation': {'text_bounds_checked': True, 'time_labels_do_not_overlap': True,
                       'vertical_tick_labels_do_not_overlap': True,
                       'ICM_TE_labels_do_not_overlap': True, 'legend_clear_of_panels': True,
                       'all_other_text_8pt_heading_10pt': True},
        **sources,
    }
    # Shared validation also verifies all text bounds/fonts and the legend margin.
    shared.validate_and_save(fig, axes, titles, [heading], legend, PREFIX)
    export = table[['sample', 'block', 'order', 'n_missing', 'hoi', 'tc', 'dtc']].rename(
        columns={'hoi': 'O_information', 'tc': 'TC', 'dtc': 'DTC'})
    export.to_csv(OUTPUT / 'oinfo_overall_values_for_plot.tsv', sep='\t', index=False)
    (OUTPUT / f'{PREFIX}_validation.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print(f'Total O-information: {WIDTH_CM:g} x {HEIGHT_CM:g} cm; heading 10 pt, other text 8 pt.', flush=True)


if __name__ == '__main__':
    main()
