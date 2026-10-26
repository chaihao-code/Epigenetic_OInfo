#!/usr/bin/env python3
"""Four full-signal occurrence plots, preserving the supplied sign-only triplet rule.

Preimplantation: all nine signals; post-E11.5: all ten signals. Cores are
recomputed with the shared replacement-background engine, in separate outputs.
"""

# Resolve repository-relative paths in retained operational JSON manifests.
_RELEASE_ROOT = next(p for p in __import__('pathlib').Path(__file__).resolve().parents
                     if (p / '.epigenetic_oinfo_root').is_file())
def _release_manifest(value):
    if isinstance(value, dict):
        return {k: _release_manifest(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_release_manifest(v) for v in value]
    if isinstance(value, str) and value.startswith(('cases/', 'common/', 'data/')):
        return str(_RELEASE_ROOT / value)
    return value
import argparse
import csv
import hashlib
import itertools
import json
import math
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/full-signal-occurrence-mpl')
os.environ.setdefault('XDG_CACHE_HOME', '/tmp/full-signal-occurrence-cache')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1')
ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[2]
OUTPUT = ROOT / 'full_signal_occurrence'
ENGINE = REPO / 'common/03_triplets_and_redundancy/compute_redundancy_core.py'
MANIFEST = OUTPUT / 'inputs/analysis_manifest.json'


def canonical_signal(token):
    # Same canonical identifiers as the core engine, without importing its
    # numerical dependencies during data validation or plotting.
    match = re.search(r'H3K\d+(?:me\d|ac)|H2A\.Z|ATAC|DNA_methylation|CpG_methyl|WGBS', token)
    name = match[0] if match else token.strip()
    return {'CpG_methyl': 'DNA_methylation', 'WGBS': 'DNA_methylation'}.get(name, name)

COMMON = ['ATAC', 'H3K4me1', 'H3K4me3', 'H3K27ac', 'H3K9me3',
          'H3K27me3', 'H3K36me3', 'DNA_methylation']
SIGNALS = dict(preimplantation=COMMON + ['H2A.Z'],
               post_E11p5=COMMON + ['H3K9ac', 'H3K4me2'])
DISPLAY_ORDER = ['DNA_methylation', 'H3K9me3', 'H3K27me3', 'ATAC', 'H2A.Z',
                 'H3K4me1', 'H3K4me2', 'H3K4me3', 'H3K27ac', 'H3K9ac',
                 'H3K36me3']
PERIOD_LABEL = dict(preimplantation='Preimplantation', post_E11p5='Post-E11.5')
TOP_N = 25
WIDTH_CM, HEIGHT_CM = 9., 7.5


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def table(path, rows, columns=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns or list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def read_table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def prepare():
    for folder in ['inputs', 'data', 'results', 'plots', 'logs']:
        (OUTPUT / folder).mkdir(parents=True, exist_ok=True)
    source_manifest = ROOT / 'inputs/analysis_manifest.json'
    original = _release_manifest(json.loads(source_manifest.read_text()))
    tasks, triplet_counts, block_reference = [], [], None
    for index, old in enumerate(original['tasks']):
        period = 'preimplantation' if old['period'] == 'preimplantation' else 'post_E11p5'
        sample, path = old['sample'], Path(old['input'])
        meta_path = path.with_name(sample + '_metadata.json')
        meta = _release_manifest(json.loads(meta_path.read_text()))
        signals = frozenset(SIGNALS[period])
        assert meta['n_features'] == len(signals)
        assert {canonical_signal(s) for s in meta['feature_names']} == signals
        assert meta['n_bins'] == 5 and len(meta['blocks']) == 40
        if block_reference is None:
            block_reference = set(meta['blocks'])
        assert set(meta['blocks']) == block_reference
        expected = {frozenset(combo) for k in range(3, len(signals) + 1)
                    for combo in itertools.combinations(sorted(signals), k)}
        by_block = defaultdict(set)
        counter = Counter()
        with path.open() as handle:
            for row in csv.DictReader(handle):
                order = int(row['order'])
                if order < 3:
                    continue
                subset = frozenset(canonical_signal(s) for s in row['remaining_names'].split('/'))
                assert len(subset) == order and subset <= signals
                assert subset not in by_block[row['block']], (sample, row['block'], subset)
                by_block[row['block']].add(subset)
                value = float(row['hoi'])
                assert math.isfinite(value)
                if order == 3:
                    assert (row['interpretation'] == 'Synergy') == (value < 0)
                    if value < 0:
                        counter['/'.join(sorted(subset))] += 1
        assert set(by_block) == block_reference
        assert all(subsets == expected for subsets in by_block.values()), sample
        task = dict(index=index, sample=sample, stage=old['stage'], period=period,
                    input=str(path), input_sha256=sha(path),
                    input_metadata=str(meta_path), input_metadata_sha256=sha(meta_path),
                    signals=SIGNALS[period], n_blocks=40,
                    n_subsets_per_block=len(expected), n_triplets_per_block=math.comb(len(signals), 3),
                    n_negative_triplet_occurrences=sum(counter.values()))
        tasks.append(task)
        for combo in itertools.combinations(sorted(signals), 3):
            name = '/'.join(combo)
            triplet_counts.append(dict(period=period, stage=old['stage'], sample=sample,
                                       combination=name, count=counter[name], n_blocks=40))
        print(f'Validated {sample}: {len(signals)} signals; {sum(counter.values())} negative triplets', flush=True)
    assert Counter(t['period'] for t in tasks) == {'preimplantation': 7, 'post_E11p5': 66}
    counts_path = OUTPUT / 'data/triplet_combinations_by_sample.tsv'
    table(counts_path, triplet_counts)
    manifest = dict(status='inputs_validated', tasks=tasks, signals=SIGNALS,
                    source_manifest=str(source_manifest), source_manifest_sha256=sha(source_manifest),
                    blocks=sorted(block_reference), n_blocks=40,
                    triplet_criterion='observed O-information < 0 (interpretation == Synergy); no P-value/FDR filter, matching supplied script',
                    triplet_FDR_applied=False,
                    core_criterion='positive O-information; modified z > 2.5; GMM BIC or ANOVA eta_new > 0.5; promotion and maximal-set pruning',
                    core_FDR_applied=False, core_min_order=3,
                    core_engine=str(ENGINE), core_engine_sha256=sha(ENGINE),
                    triplet_counts_sha256=sha(counts_path), top_n=TOP_N,
                    pooling='All 40 regions and all available samples within each period; each sample x block x combination is one occurrence',
                    note='Full signal universes differ between periods. Existing shared-eight-signal analyses are unchanged.')
    MANIFEST.write_text(json.dumps(manifest, indent=2) + '\n')


def core_task(index):
    manifest = _release_manifest(json.loads(MANIFEST.read_text()))
    task = manifest['tasks'][index]
    assert sha(ENGINE) == manifest['core_engine_sha256']
    assert sha(task['input']) == task['input_sha256']
    path = OUTPUT / 'results' / task['period'] / task['sample'] / 'redundancy_cores.csv'
    if path.exists() and path.with_suffix('.metadata.json').exists():
        meta = _release_manifest(json.loads(path.with_suffix('.metadata.json').read_text()))
        if (meta.get('status') == 'complete' and meta.get('input_sha256') == task['input_sha256']
                and meta.get('program_sha256') == manifest['core_engine_sha256']
                and set(meta.get('signals', [])) == set(task['signals'])):
            print(f'Using validated cached cores: {task["sample"]}', flush=True)
            return
    command = [sys.executable, str(ENGINE), '--input', task['input'], '--output', str(path),
               '--sample', task['sample'], '--signals', *task['signals'], '--min-order', '3',
               '--z-threshold', '2.5', '--eta-threshold', '0.5', '--expected-blocks', '40']
    subprocess.run(command, check=True)


def collect():
    manifest = _release_manifest(json.loads(MANIFEST.read_text()))
    assert sha(ENGINE) == manifest['core_engine_sha256']
    triplet_file = OUTPUT / 'data/triplet_combinations_by_sample.tsv'
    assert sha(triplet_file) == manifest['triplet_counts_sha256']
    counters = {(period, kind): Counter() for period in SIGNALS for kind in ['synergistic_triplet', 'redundancy_core']}
    task_map = {t['sample']: t for t in manifest['tasks']}
    sample_triplets = Counter()
    for row in read_table(triplet_file):
        task = task_map[row['sample']]
        assert row['period'] == task['period'] and int(row['n_blocks']) == 40
        value = int(row['count'])
        assert 0 <= value <= 40
        if value:
            counters[row['period'], 'synergistic_triplet'][row['combination']] += value
        sample_triplets[row['sample']] += value
    core_rows, sample_counts = [], []
    for task in manifest['tasks']:
        path = OUTPUT / 'results' / task['period'] / task['sample'] / 'redundancy_cores.csv'
        metadata = _release_manifest(json.loads(path.with_suffix('.metadata.json').read_text()))
        assert metadata['status'] == 'complete' and metadata['n_blocks'] == 40
        assert metadata['input_sha256'] == task['input_sha256']
        assert metadata['program_sha256'] == manifest['core_engine_sha256']
        assert set(metadata['signals']) == set(task['signals'])
        assert metadata['z_threshold'] == 2.5 and metadata['eta_threshold'] == .5
        with path.open() as handle:
            rows = list(csv.DictReader(handle))
        assert len(rows) == metadata['n_cores']
        seen, blocks = set(), defaultdict(list)
        for row in rows:
            subset = frozenset(row['signals'].split('/'))
            assert len(subset) == int(row['order']) and subset <= set(task['signals'])
            assert float(row['hoi']) > 0
            if row['promoted'] != 'True':
                assert float(row['mod_z']) > 2.5 and row['passes_modality'] == 'True'
            key = row['block'], subset
            assert key not in seen
            seen.add(key)
            blocks[row['block']].append(subset)
            counters[task['period'], 'redundancy_core']['/'.join(sorted(subset))] += 1
            core_rows.append(dict(period=task['period'], stage=task['stage'], **row))
        for subsets in blocks.values():
            assert not any(a < b for a in subsets for b in subsets)
        assert sample_triplets[task['sample']] == task['n_negative_triplet_occurrences']
        sample_counts.append(dict(period=task['period'], stage=task['stage'], sample=task['sample'],
                                  n_signals=len(task['signals']), n_blocks=40,
                                  n_negative_triplets=sample_triplets[task['sample']], n_cores=len(rows)))
    table(OUTPUT / 'data/redundancy_cores.tsv', core_rows)
    table(OUTPUT / 'data/counts_by_sample.tsv', sample_counts)
    summaries = []
    for (period, kind), counter in counters.items():
        total = sum(counter.values())
        assert total > 0
        ranked = sorted(counter.items(), key=lambda item: (-item[1], item[0]))
        cumulative, rows = 0, []
        for rank, (combo, count) in enumerate(ranked, 1):
            cumulative += count
            rows.append(dict(rank=rank, combination=combo, count=count,
                             proportion_pct=count / total * 100,
                             cumulative_pct=cumulative / total * 100,
                             total_occurrences=total, selected_top25=rank <= TOP_N))
        table(OUTPUT / 'data' / f'{kind}_combination_occurrence_{period}_all.tsv', rows)
        table(OUTPUT / 'data' / f'{kind}_combination_occurrence_{period}_top25.tsv', rows[:TOP_N])
        summaries.append(dict(period=period, kind=kind, n_signals=len(SIGNALS[period]),
                              n_samples=sum(t['period'] == period for t in manifest['tasks']),
                              total_occurrences=total, n_distinct_combinations=len(rows),
                              n_plotted=min(len(rows), TOP_N),
                              top25_coverage_pct=rows[min(TOP_N, len(rows)) - 1]['cumulative_pct']))
    table(OUTPUT / 'data/occurrence_summary.tsv', summaries)
    return manifest, summaries


def style():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    font = font_manager.findfont('Arial', fallback_to_default=False)
    plt.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'axes.labelsize': 8,
                         'axes.titlesize': 8, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
                         'legend.fontsize': 8, 'pdf.fonttype': 42, 'ps.fonttype': 42,
                         'axes.linewidth': .5, 'savefig.facecolor': 'white'})
    return plt, font


def draw_panel(fig, bounds, period, kind, letter=None):
    import numpy as np
    from matplotlib.ticker import MaxNLocator
    left, bottom, width, height = bounds
    def rect(x, y, w, h):
        return [left + x / WIDTH_CM * width, bottom + y / HEIGHT_CM * height,
                w / WIDTH_CM * width, h / HEIGHT_CM * height]
    rows = read_table(OUTPUT / 'data' / f'{kind}_combination_occurrence_{period}_top25.tsv')
    n = len(rows)
    marks = [s for s in DISPLAY_ORDER if s in SIGNALS[period]]
    assert marks[:3] == ['DNA_methylation', 'H3K9me3', 'H3K27me3']
    # Physical margins reserve space for 8 pt mark names and both core axes.
    ax = fig.add_axes(rect(2.35, 3.65, 5.40, 2.80))
    dot = fig.add_axes(rect(2.35, .48, 5.40, 2.93), sharex=ax)
    x = np.arange(n)
    core = kind == 'redundancy_core'
    bar_values = [float(row['proportion_pct']) if core else int(row['count']) for row in rows]
    ax.bar(x, bar_values, width=.72, color='#78B5D3', edgecolor='#303030', linewidth=.35)
    ax.set_ylabel('Proportion (%)' if core else 'Occurrence count', labelpad=2)
    ax.set_ylim(0, max(bar_values) * 1.18)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4, integer=not core))
    ax.set_yticks([v for v in ax.get_yticks() if ax.get_ylim()[0] <= v <= ax.get_ylim()[1]])
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(axis='x', bottom=False, labelbottom=False)
    ax.tick_params(axis='y', width=.5, length=2, pad=2)
    if core:
        cumulative = ax.twinx()
        cumulative.plot(x, [float(row['cumulative_pct']) for row in rows],
                        color='#BF3851', linewidth=.85, marker='o', markersize=2.1)
        cumulative.set_ylim(0, 105)
        cumulative.set_yticks([0, 50, 100])
        cumulative.set_ylabel('Cumulative (%)', color='#BF3851', labelpad=2)
        cumulative.tick_params(axis='y', colors='#BF3851', width=.5, length=2, pad=2)
        cumulative.spines['top'].set_visible(False)
        cumulative.spines['right'].set_color('#BF3851')
        cumulative.axhline(100, color='#AAAAAA', linewidth=.45, linestyle='--', zorder=0)
    for i, mark in enumerate(marks):
        y = len(marks) - 1 - i
        dot.axhline(y, color='#DEDEDE', linewidth=.4, zorder=0)
        dot.scatter(x, np.full(n, y), color='#E0E0E0', s=5, linewidths=0, zorder=1)
    for j, row in enumerate(rows):
        subset = set(row['combination'].split('/'))
        ys = [len(marks) - 1 - i for i, mark in enumerate(marks) if mark in subset]
        assert len(ys) == len(subset) >= 3
        dot.plot([j, j], [min(ys), max(ys)], color='#303030', linewidth=.65, zorder=2)
        dot.scatter(np.full(len(ys), j), ys, color='#151515', s=8, linewidths=0, zorder=3)
    dot.set_xlim(-.6, n - .4)
    dot.set_ylim(-.55, len(marks) - .45)
    dot.set_yticks(range(len(marks)), [m.replace('_', ' ') for m in marks[::-1]])
    dot.set_xticks([])
    dot.tick_params(axis='y', length=0, pad=4)
    dot.spines[:].set_visible(False)
    title = PERIOD_LABEL[period] + '\n' + ('Redundancy cores' if core else 'Synergistic triplets')
    # Align titles with the data, keeping them independent of the twins' y labels.
    tx = left + (2.35 + 5.40 / 2) / WIDTH_CM * width
    ty = bottom + 7.22 / HEIGHT_CM * height
    fig.text(tx, ty, title, ha='center', va='top', fontsize=8, linespacing=1.15)
    if letter:
        fig.text(left + .12 / WIDTH_CM * width, bottom + 7.22 / HEIGHT_CM * height,
                 letter, fontsize=8, fontweight='bold', va='top')
    return dict(period=period, kind=kind, n_columns=n, n_signal_rows=len(marks),
                signal_order_top_to_bottom=[m.replace('_', ' ') for m in marks],
                dots_match_all_combination_members=True)


def validate_figure(fig):
    from matplotlib.text import Text
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    errors = []
    for item in fig.findobj(Text):
        if item.get_visible() and item.get_text():
            assert item.get_fontsize() == 8, item.get_text()
            # Tick locators also generate invisible, out-of-range ticks.
            if item in [label for ax in fig.axes for label in ax.get_yticklabels()]:
                ax = item.axes
                if ax is not None and not ax.get_ylim()[0] <= item.get_position()[1] <= ax.get_ylim()[1]:
                    continue
            bb = item.get_window_extent(renderer)
            if bb.x0 < -.5 or bb.y0 < -.5 or bb.x1 > fig.bbox.width + .5 or bb.y1 > fig.bbox.height + .5:
                errors.append((item.get_text(), tuple(bb.bounds)))
    if errors:
        raise ValueError(f'Figure text outside canvas: {errors}')
    for ax in fig.axes:
        labels = [item for item in ax.get_yticklabels() if item.get_visible() and item.get_text()]
        for first, second in itertools.combinations(labels, 2):
            assert not first.get_window_extent(renderer).overlaps(second.get_window_extent(renderer)), (
                first.get_text(), second.get_text())


def plot(plots_only=False):
    (OUTPUT / 'plots').mkdir(parents=True, exist_ok=True)
    if plots_only:
        manifest = _release_manifest(json.loads(MANIFEST.read_text()))
        summaries = read_table(OUTPUT / 'data/occurrence_summary.tsv')
    else:
        manifest, summaries = collect()
    plt, font = style()
    panels = []
    for period in SIGNALS:
        for kind in ['synergistic_triplet', 'redundancy_core']:
            fig = plt.figure(figsize=(WIDTH_CM / 2.54, HEIGHT_CM / 2.54))
            panels.append(draw_panel(fig, (0, 0, 1, 1), period, kind))
            validate_figure(fig)
            prefix = OUTPUT / 'plots' / f'{kind}_combination_occurrence_{period}_top25'
            for extension in ['pdf', 'png']:
                fig.savefig(prefix.with_suffix('.' + extension), dpi=300)
            plt.close(fig)
    # Rows compare core/triplet types; columns compare developmental periods.
    fig = plt.figure(figsize=(2 * WIDTH_CM / 2.54, 2 * HEIGHT_CM / 2.54))
    combined_panels = []
    for kind, y, letters in [('redundancy_core', .5, ['a', 'b']), ('synergistic_triplet', 0, ['c', 'd'])]:
        for i, period in enumerate(['preimplantation', 'post_E11p5']):
            panel = draw_panel(fig, (i * .5, y, .5, .5), period, kind, letters[i])
            combined_panels.append(dict(letter=letters[i], row='top' if y else 'bottom',
                                        column='left' if i == 0 else 'right', **panel))
    validate_figure(fig)
    for extension in ['pdf', 'png']:
        fig.savefig(OUTPUT / 'plots' / ('full_signal_triplets_and_cores_2x2.' + extension), dpi=300)
    plt.close(fig)
    validation = dict(status='complete', summaries=summaries, panels=panels,
                      combined_panels=combined_panels,
                      combined_rows=['Redundancy cores', 'Synergistic triplets'],
                      combined_columns=['Preimplantation', 'Post-E11.5'],
                      individual_size_cm=[WIDTH_CM, HEIGHT_CM], combined_size_cm=[2 * WIDTH_CM, 2 * HEIGHT_CM],
                      all_font_sizes_pt=8, font_family='Arial', font_file=font,
                      text_bounds_checked=True, y_tick_overlap_checked=True,
                      full_signal_cores_recomputed=not plots_only,
                      plotting_from_retained_tables=plots_only,
                      triplet_FDR_applied=False, triplet_criterion=manifest['triplet_criterion'],
                      top_n=TOP_N, n_core_samples=73,
                      script_sha256=sha(Path(__file__)), core_engine_sha256=sha(ENGINE))
    (OUTPUT / 'data/plot_validation.json').write_text(json.dumps(validation, indent=2) + '\n')
    manifest['status'] = 'complete'
    MANIFEST.write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(validation, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--prepare', action='store_true')
    group.add_argument('--core-task', type=int)
    group.add_argument('--plot', action='store_true')
    group.add_argument('--plots-only', action='store_true',
                       help='Redraw retained combination-occurrence tables without recollecting per-sample cores')
    args = parser.parse_args()
    if args.prepare:
        prepare()
    elif args.core_task is not None:
        core_task(args.core_task)
    else:
        plot(plots_only=args.plots_only)


if __name__ == '__main__':
    main()
