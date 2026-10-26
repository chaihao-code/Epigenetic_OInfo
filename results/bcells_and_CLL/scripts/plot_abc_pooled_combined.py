#!/usr/bin/env python3
"""Stack A/B/C donor means after verified global BH correction of all triplet tests."""
import hashlib
import json
from plot_abc_triplet_classes import (
    PROJECT, OUT, PLOTS, BLOCKS, REGION_LABELS, SIGNALS, ANCHORS,
    canonical, classify, plt, font_manager, Path, pd, np,
)
from matplotlib.text import Text

CRITERION = 'observed_hoi < 0 AND q_down_global_BH < 0.05'
GROUP_SIZES = {'R': 8, 'S': 5, 'pB': 4}
CAPACITIES = {'A': 12, 'B': 19, 'C': 4}


def corrected_group_means():
    source = PROJECT / 'results_oinfo_triplet_n_bins_5/all_samples_triplet_significance.tsv'
    df = pd.read_csv(source, sep='\t')
    assert len(df) == 23800 and df['sample'].nunique() == 17
    assert df.order.eq(3).all() and df.n_perm.eq(1000).all()
    assert set(df.block) == set(BLOCKS)
    df['triplet_signals'] = df.subset_names.map(canonical)
    df['class'] = df.triplet_signals.map(classify)
    assert not df.duplicated(['sample', 'block', 'triplet_signals']).any()
    assert df.groupby(['sample', 'block']).size().eq(35).all()
    # Recompute one global lower-tail BH family before any sign/class/group filtering.
    p = df.p_down.to_numpy(dtype=float)
    assert np.isfinite(p).all() and ((p >= 0) & (p <= 1)).all()
    order = np.argsort(p, kind='stable')
    ranked_q = np.minimum(1., np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1])
    q = np.empty(len(p), dtype=float)
    q[order] = ranked_q
    discrepancy = float(np.max(np.abs(q - df.q_down_global_BH.to_numpy(dtype=float))))
    assert discrepancy < 1e-12
    df['synergistic_FDR005'] = df.observed_hoi.lt(0) & (q < .05)
    df['synergistic_nominal'] = df.observed_hoi.lt(0) & df.p_down.lt(.05)
    assert not (df.synergistic_FDR005 & ~df.synergistic_nominal).any()
    mapping = pd.read_csv(OUT / 'sample_labels.tsv', sep='\t')
    assert not mapping['sample'].duplicated().any() and set(mapping['sample']) == set(df['sample'])
    assert mapping.pooled_group.value_counts().to_dict() == GROUP_SIZES
    counts = df.groupby(['sample', 'block', 'class']).agg(
        n_tested=('class', 'size'), count_FDR005=('synergistic_FDR005', 'sum'),
        count_nominal=('synergistic_nominal', 'sum')).reset_index()
    assert len(counts) == 17 * 40 * 3
    assert counts.n_tested.eq(counts['class'].map(CAPACITIES)).all()
    counts = counts.merge(mapping[['sample', 'pooled_group']], on='sample', validate='many_to_one')
    means = counts.groupby(['class', 'pooled_group', 'block']).agg(
        n_samples=('sample', 'nunique'), mean_count=('count_FDR005', 'mean'),
        std_count=('count_FDR005', 'std'), mean_count_nominal=('count_nominal', 'mean')).reset_index()
    assert len(means) == 3 * 3 * 40
    assert means.n_samples.eq(means.pooled_group.map(GROUP_SIZES)).all()
    # Existing sample-level FDR summaries provide an independent aggregation check.
    reference = pd.read_csv(OUT / 'counts_by_sample_region_class.tsv', sep='\t')
    matched = counts.merge(reference[['sample', 'block', 'class', 'count_FDR005']],
                           on=['sample', 'block', 'class'], validate='one_to_one', suffixes=('', '_reference'))
    assert len(matched) == len(counts) and matched.count_FDR005.eq(matched.count_FDR005_reference).all()
    group_order = {g:i for i,g in enumerate(GROUP_SIZES)}
    block_order = {b:i for i,b in enumerate(BLOCKS)}
    means = means.sort_values(['class', 'pooled_group', 'block'], key=lambda s:
                             s.map(group_order) if s.name == 'pooled_group' else
                             s.map(block_order) if s.name == 'block' else s).reset_index(drop=True)
    validation = dict(status='complete', source=str(source), source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                      criterion=CRITERION, BH_family='All 23,800 sample × annotation block × triplet tests',
                      BH_family_size=len(df), BH_applied_before_sign_or_class_filtering=True,
                      global_BH_independently_recomputed=True, maximum_q_discrepancy=discrepancy,
                      n_permutations=1000, signals=SIGNALS, anchor_signals=sorted(ANCHORS),
                      class_capacities=CAPACITIES, sample_groups=GROUP_SIZES,
                      aggregation='Equal-weight donor mean within R, S and pB separately; no group-difference test',
                      n_nominal_synergistic=int(df.synergistic_nominal.sum()),
                      n_FDR_synergistic=int(df.synergistic_FDR005.sum()),
                      n_pooled_cells=len(means),
                      n_pooled_cells_changed=int((means.mean_count != means.mean_count_nominal).sum()),
                      all_sample_level_FDR_counts_match_existing_summaries=True)
    return means, validation


def main():
    means, validation = corrected_group_means()
    font = Path(os.environ.get("BCELLS_ARIAL_FONT", str(Path.home() / ".local/share/fonts/arial.ttf")))
    if font.exists():
        font_manager.fontManager.addfont(str(font))
    plt.rcParams.update({"font.family": "Arial" if font.exists() else "DejaVu Sans",
        "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
        "xtick.labelsize": 8, "ytick.labelsize": 8,
        "pdf.fonttype": 42, "ps.fonttype": 42})
    width_cm = 18.0
    heatmap_width_cm = 14.6
    cell_cm = heatmap_width_cm / 40
    panel_height_cm = cell_cm * 3
    top_margin_cm, bottom_margin_cm, panel_gap_cm = .60, 2.95, .60
    height_cm = top_margin_cm + 3*panel_height_cm + 2*panel_gap_cm + bottom_margin_cm
    fig = plt.figure(figsize=(width_cm/2.54, height_cm/2.54))
    capacities = CAPACITIES
    export = []
    axes = []
    cell_texts = []
    for panel_index, cl in enumerate("ABC"):
        bottom_cm = height_cm - top_margin_cm - panel_height_cm - panel_index*(panel_height_cm+panel_gap_cm)
        data = means[means['class'].eq(cl)].drop(columns='class').copy()
        assert len(data) == 120
        assert data.n_samples.eq(data.pooled_group.map({"R": 8, "S": 5, "pB": 4})).all()
        matrix = data.pivot(index="pooled_group", columns="block", values="mean_count").reindex(
            index=["R", "S", "pB"], columns=BLOCKS)
        assert matrix.shape == (3, 40) and matrix.notna().all().all()
        assert np.all((matrix.to_numpy() >= 0) & (matrix.to_numpy() <= capacities[cl]))
        data.insert(0, "class", cl)
        export.append(data)
        ax = fig.add_axes([1.8/width_cm, bottom_cm/height_cm,
            heatmap_width_cm/width_cm, panel_height_cm/height_cm])
        axes.append(ax)
        im = ax.imshow(matrix.to_numpy(), cmap="YlOrRd", vmin=0, vmax=capacities[cl],
            aspect="equal", interpolation="nearest")
        ax.set_yticks(range(3), ["R", "S", "pB"])
        ax.tick_params(axis="y", length=0, labelsize=8)
        ax.set_xticks([])
        for k in range(1, 8):
            ax.axvline(k*5-.5, color="gray", lw=.5)
        for i in range(3):
            for j in range(40):
                value = matrix.iloc[i, j]
                cell_texts.append(ax.text(j, i, f"{value:.1f}", ha="center", va="center", fontsize=3.5,
                    color="white" if value >= capacities[cl]*.62 else "black"))
        ax.set_title(f"Triplet {cl}", loc="left", fontsize=8, pad=4)
        cax = fig.add_axes([16.65/width_cm, bottom_cm/height_cm, .13/width_cm, panel_height_cm/height_cm])
        ticks = sorted(set([0, round(capacities[cl]/3), round(2*capacities[cl]/3), capacities[cl]]))
        cb = fig.colorbar(im, cax=cax, ticks=ticks)
        cb.ax.tick_params(labelsize=8, length=2, pad=2)
        if cl == "C":
            for k, label in enumerate(REGION_LABELS):
                ax.annotate("", xy=(k*5+4.4, -.08), xytext=(k*5-.4, -.08),
                    xycoords=ax.get_xaxis_transform(), arrowprops=dict(arrowstyle="<->", color="gray", lw=.5))
                ax.text(k*5+2, -.23, label, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=8)
            for j in range(40):
                ax.text(j, -1.05/panel_height_cm, ["CpG island", "LINE", "SINE", "LTR", "None"][j % 5],
                    transform=ax.get_xaxis_transform(), ha="center", va="top", rotation=90, fontsize=8)
    shared_color_label_y = (height_cm - top_margin_cm - (3*panel_height_cm + 2*panel_gap_cm)/2) / height_cm
    fig.text(17.65/width_cm, shared_color_label_y, 'Mean triplet count',
             rotation=90, ha='center', va='center', fontsize=8)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for ax in axes:
        box = ax.get_window_extent()
        assert np.isclose(box.width/40, box.height/3, atol=1e-9), "Cells must be square"
    fig.text(.5, .10/height_cm,
             'O-information < 0; global BH FDR < 0.05 (23,800 tests); 1,000 permutations; equal donor means.',
             ha='center', va='bottom', fontsize=8)
    fig.canvas.draw()
    cells = set(cell_texts)
    for text in fig.findobj(Text):
        if text.get_visible() and text.get_text():
            assert text.get_fontsize() == (3.5 if text in cells else 8), text.get_text()
            bbox = text.get_window_extent(renderer)
            assert bbox.x0 >= 0 and bbox.y0 >= 0 and bbox.x1 <= fig.bbox.width and bbox.y1 <= fig.bbox.height, text.get_text()
            if text in cells:
                panel_box = text.axes.get_window_extent(renderer)
                assert bbox.width < panel_box.width/40 and bbox.height < panel_box.height/3
    layout = dict(width_cm=width_cm, height_cm=height_cm,
        cell_width_cm=cell_cm, cell_height_cm=cell_cm, panel_rows=3, panel_columns=40,
        titles=["Triplet A", "Triplet B", "Triplet C"],
        sample_groups={"R": 8, "S": 5, "pB": 4},
        displayed_row_labels=['R', 'S', 'pB'], sample_sizes_shown=False,
        label_font_size_pt=8, cell_number_font_size_pt=3.5, shared_colorbar_label=True,
        criterion=CRITERION, BH_family_size=23800, n_permutations=1000,
        color_limits={"A": [0, 12], "B": [0, 19], "C": [0, 4]})
    (PLOTS / "ABC_pooled_average_layout.json").write_text(json.dumps(layout, indent=2))
    print(f"Figure: {width_cm:.3f} x {height_cm:.3f} cm; square cells: {cell_cm:.3f} x {cell_cm:.3f} cm")
    stem = "ABC_pooled_average_heatmap"
    for ext in ["png", "pdf"]:
        path = PLOTS / f"{stem}.{ext}"
        fig.savefig(path, dpi=300)
        print(path)
    plt.close(fig)
    pd.concat(export, ignore_index=True).to_csv(PLOTS / "ABC_pooled_average_counts.tsv", sep="\t", index=False)
    validation.update(figure_size_cm=[width_cm, height_cm], cells_are_square=True, text_bounds_checked=True,
                      label_font_size_pt=8, cell_number_font_size_pt=3.5,
                      displayed_row_labels=['R', 'S', 'pB'], sample_sizes_shown=False,
                      cell_numbers_fit_cells=True, shared_colorbar_label=True)
    (PLOTS / 'ABC_pooled_average_FDR005_validation.json').write_text(json.dumps(validation, indent=2) + '\n')
    config_path = OUT / 'plot_config.json'
    config = json.loads(config_path.read_text())
    config['pooled_combined_criterion'] = CRITERION + '; global BH family of 23800 tests'
    config['pooled_combined_script'] = '../scripts/plot_abc_pooled_combined.py'
    config_path.write_text(json.dumps(config, indent=2) + '\n')
    print('Global BH verified: nominal discoveries =', validation['n_nominal_synergistic'],
          '; FDR discoveries =', validation['n_FDR_synergistic'],
          '; changed pooled cells =', validation['n_pooled_cells_changed'])


if __name__ == "__main__":
    main()
