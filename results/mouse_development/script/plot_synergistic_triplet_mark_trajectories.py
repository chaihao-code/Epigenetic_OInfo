#!/usr/bin/env python3
"""Count each mark's significant triplet occurrences across full development.

Seven signals, all 40 genomic blocks, all 82 samples and 17 stages. Each
sample/block retains the existing 35-test BH correction; negative O-information
and q <= 0.05 define a significant synergistic triplet. No permutations rerun.
"""
import argparse
import csv
import gzip
import hashlib
import itertools
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean, stdev

from plot_triplet_A_B_C_counts_trajectories import OUTPUT, STAGES

SOURCE = OUTPUT / 'all_triplets_fdr_10000.tsv.gz'
PREFIX = 'synergistic_triplet_mark_trajectories_7signals'
PALETTE_SOURCE = OUTPUT.parents[1] / 'H3K36me3_developmental_switch/redundancy_core/inputs/analysis_manifest.json'
MARKS = ['ATAC', 'H3K4me1', 'H3K4me3', 'H3K27ac', 'H3K9me3', 'H3K27me3', 'CpG_methyl']
DISPLAY = {m: 'DNA methylation' if m == 'CpG_methyl' else m for m in MARKS}
FIG_W_CM, FIG_H_CM = 8.5, 5


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_tsv(path, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def positions():
    # Early stages are categories; later spacing is one unit per day, with
    # a two-unit gap between E7.5 and E11.5. ICM and TE stay separate.
    xs = {stage: i for i, stage in enumerate(STAGES[:10])}
    xs.update({stage: 11 + (19.5 if stage == 'P0' else float(stage[1:])) - 11.5
               for stage in STAGES[10:]})
    return xs


def collect_counts():
    universe = {frozenset(t) for t in itertools.combinations(MARKS, 3)}
    assert len(universe) == 35
    families = defaultdict(list)
    counts, total_triplets = Counter(), Counter()
    sample_stage, triplet_cache = {}, {}
    input_rows = 0
    for_row = gzip.open(SOURCE, 'rt')
    with for_row as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            input_rows += 1
            sample, stage, block = row['sample'], row['stage'], row['block']
            assert stage == sample.split('_')[0] and stage in STAGES
            assert sample_stage.setdefault(sample, stage) == stage
            assert int(row['order']) == 3 and int(row['n_perm']) == 10000
            if row['triplet'] not in triplet_cache:
                marks = frozenset(row['triplet'].split('/'))
                assert len(marks) == 3 and marks in universe
                triplet_cache[row['triplet']] = marks
            marks = triplet_cache[row['triplet']]
            p, q, observed = float(row['p_down']), float(row['q_down_sample_block_bh']), float(row['observed_hoi'])
            assert 0 <= p <= 1 and 0 <= q <= 1 and math.isfinite(observed)
            families[sample, block].append((marks, p, q))
            significant = observed < 0 and q <= .05
            assert (row['significant_negative'] == 'True') == significant
            if significant:
                total_triplets[sample] += 1
                for mark in marks:
                    counts[sample, mark] += 1
    assert len(sample_stage) == 82 and set(sample_stage.values()) == set(STAGES)
    all_blocks = {block for sample, block in families}
    assert len(all_blocks) == 40 and input_rows == 82 * 40 * 35
    for sample in sample_stage:
        assert {b for s, b in families if s == sample} == all_blocks
    max_q_error = 0.
    for key, tests in families.items():
        assert len(tests) == 35 and {t[0] for t in tests} == universe, key
        ordered = sorted(tests, key=lambda t: t[1])
        running_q = 1.
        for i in range(34, -1, -1):
            running_q = min(running_q, ordered[i][1] * 35 / (i + 1))
            error = abs(running_q - ordered[i][2])
            max_q_error = max(max_q_error, error)
            assert error < 1e-12, (key, i, error)
    per_sample = []
    stage_samples = {stage: sorted(s for s, value in sample_stage.items() if value == stage)
                     for stage in STAGES}
    for sample in sorted(sample_stage):
        assert sum(counts[sample, mark] for mark in MARKS) == 3 * total_triplets[sample]
        for mark in MARKS:
            count = counts[sample, mark]
            assert 0 <= count <= 40 * math.comb(6, 2)
            per_sample.append(dict(sample=sample, stage=sample_stage[sample], mark=mark,
                                   display_mark=DISPLAY[mark], count=count, n_blocks=40,
                                   total_significant_triplets=total_triplets[sample]))
    xs = positions()
    summary = []
    for stage in STAGES:
        for mark in MARKS:
            values = [counts[sample, mark] for sample in stage_samples[stage]]
            n = len(values)
            sd = stdev(values) if n > 1 else 0.
            summary.append(dict(stage=stage, mark=mark, display_mark=DISPLAY[mark],
                                x_pos=xs[stage], mean_count=mean(values), std_count=sd,
                                sem_count=sd / math.sqrt(n), n=n))
    # Independent agreement with the preserved A+B+C total count trajectories.
    reference = OUTPUT / 'total_synergistic_triplet_count_by_stage.tsv'
    stage_reference = {r['stage']: r for r in read_tsv(reference)}
    assert set(stage_reference) == set(STAGES)
    for stage in STAGES:
        ref = stage_reference[stage]
        assert int(ref['n_samples']) == len(stage_samples[stage])
        assert abs(mean(total_triplets[s] for s in stage_samples[stage]) - float(ref['mean_total'])) < 1e-10
        mark_total = sum(float(r['mean_count']) for r in summary if r['stage'] == stage)
        assert abs(mark_total - 3 * float(ref['mean_total'])) < 1e-10
    metadata = dict(
        source=str(SOURCE), source_sha256=sha256(SOURCE), source_rows=input_rows,
        signals=MARKS, n_signals=7, H3K36me3_included=False,
        n_blocks=40, blocks=sorted(all_blocks), n_samples=82, n_stages=17, stages=STAGES,
        samples_per_stage={stage: len(stage_samples[stage]) for stage in STAGES},
        sample_selection='All available samples per stage; ICM and TE separate; gastrulation included',
        criterion='observed_hoi < 0 and q_down_sample_block_bh <= 0.05',
        BH_family='35 triplets separately within each sample x genomic block',
        n_permutations=10000, n_complete_test_families=len(families),
        BH35_recomputed_and_matched=True, max_abs_BH_q_error=max_q_error,
        counting='For each sample and mark, sum significant triplets containing that mark across all 40 blocks; the same combination in different blocks counts separately',
        normalization='None; raw occurrence counts per sample',
        aggregation='Arithmetic mean across available samples/tissues at each stage',
        error_bars='SEM from sample SD (ddof=1); one-sample SEM=0',
        zero_samples_included=True,
        mass_balance_mark_counts_equal_three_times_triplet_count=True,
        preserved_total_count_table=str(reference), preserved_total_count_sha256=sha256(reference),
        stage_mean_totals_match_preserved_ABC_counts=True,
        n_sample_rows=len(per_sample), n_summary_rows=len(summary),
        no_permutations_rerun=True,
    )
    print(f'Counted {input_rows} tests: 82 samples, 40 blocks, 7 marks, 17 stages; BH35 and A+B+C totals matched.', flush=True)
    return per_sample, summary, metadata


def plot(summary, metadata):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    from matplotlib.text import Text
    from matplotlib.ticker import MaxNLocator

    font = Path('/lustre/home/2101110354/.local/share/fonts/arial.ttf')
    if font.exists():
        font_manager.fontManager.addfont(str(font))
    plt.rcParams.update({'font.family': 'Arial' if font.exists() else 'DejaVu Sans',
                         'font.size': 8, 'axes.labelsize': 8, 'axes.titlesize': 8,
                         'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 8,
                         'pdf.fonttype': 42, 'ps.fonttype': 42})
    reference_palette = json.loads(PALETTE_SOURCE.read_text())['palette']
    palette = {mark: reference_palette['DNA_methylation' if mark == 'CpG_methyl' else mark]
               for mark in MARKS}
    lookup = {(r['stage'], r['mark']): r for r in summary}
    assert len(lookup) == 17 * 7
    fig, ax = plt.subplots(figsize=(FIG_W_CM / 2.54, FIG_H_CM / 2.54), dpi=300)
    handles = []
    for mark in MARKS:
        records = [lookup[stage, mark] for stage in STAGES]
        handles.append(ax.errorbar(
            [float(r['x_pos']) for r in records], [float(r['mean_count']) for r in records],
            yerr=[float(r['sem_count']) for r in records],
            color=palette[mark], marker='o', markersize=2.8,
            markeredgecolor='#606060', markeredgewidth=.3,
            linewidth=.9, elinewidth=.5, capsize=1.5, label=DISPLAY[mark]))
    xs = positions()
    labels = ['Zygote', '2-cell', '4-cell', '8-cell', 'Morula', 'ICM', 'TE'] + STAGES[7:]
    ax.set_xticks([xs[stage] for stage in STAGES], labels, rotation=90, ha='center', va='top')
    ax.set_xlim(-.55, xs['P0'] + .55)
    ax.set_ylim(bottom=0)
    for boundary in [6.5, 10]:
        ax.axvline(boundary, color='#999999', linestyle='--', linewidth=.6)
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(width=.5, length=2)
    ax.set_ylabel('Mean count\nper sample', fontsize=8, labelpad=2)
    ax.set_title('Synergistic triplets containing each signal', fontsize=8, pad=4)
    legend = fig.legend(handles, [DISPLAY[m] for m in MARKS], loc='lower center',
                        bbox_to_anchor=(.5, .05 / FIG_H_CM), ncol=4, frameon=False,
                        handlelength=.8, columnspacing=.5, handletextpad=.25,
                        labelspacing=.15, borderpad=0, borderaxespad=0)
    fig.subplots_adjust(left=1.42 / FIG_W_CM, right=(FIG_W_CM - .2) / FIG_W_CM,
                        top=(FIG_H_CM - .55) / FIG_H_CM, bottom=1.75 / FIG_H_CM)
    lo, hi = ax.get_ylim()
    locator = MaxNLocator(nbins=5, steps=[1, 2, 2.5, 5, 10])
    ax.set_yticks([tick for tick in locator.tick_values(lo, hi) if lo <= tick <= hi])
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    text_layout = []
    for text in fig.findobj(Text):
        if not text.get_visible() or not text.get_text():
            continue
        assert text.get_fontsize() == 8, text.get_text()
        box = text.get_window_extent(renderer)
        assert (box.x0 >= 0 and box.y0 >= 0 and box.x1 <= fig.bbox.width and box.y1 <= fig.bbox.height), (
                    'Text outside canvas', text.get_text(), box.bounds)
        text_layout.append({'text': text.get_text(), 'font_pt': text.get_fontsize(),
                            'bbox_cm': [v / fig.dpi * 2.54 for v in box.bounds]})
    ticks = [t.get_window_extent(renderer) for t in ax.get_xticklabels()]
    assert not any(a.overlaps(b) for i, a in enumerate(ticks) for b in ticks[i + 1:]), 'Stage labels overlap'
    legend_box = legend.get_window_extent(renderer)
    assert min(b.y0 for b in ticks) > legend_box.y1, 'Stage labels overlap legend'
    assert legend_box.x0 >= 0 and legend_box.x1 <= fig.bbox.width
    for ext in ['pdf', 'png']:
        fig.savefig(OUTPUT / f'{PREFIX}.{ext}', dpi=300)
    plt.close(fig)
    metadata.update(
        figure_size_cm=[FIG_W_CM, FIG_H_CM], font_size_pt=8, title_font_size_pt=8,
        palette=palette, palette_source=str(PALETTE_SOURCE),
        x_spacing='Early categorical positions 0-9; E11.5 at 11; later one unit/day, P0=19.5 days',
        ICM_TE='Separate stages', gastrulation='E6.5, E7.0 and E7.5 included',
        legend_columns=4, legend_rows=2, stage_label_rotation_deg=90,
        y_axis_label='Mean count\nper sample',
        validation={'text_bounds_checked': True, 'stage_labels_do_not_overlap': True,
                    'legend_clear_of_stage_labels': True, 'all_marks_same_line_weight': True},
        text=text_layout,
    )
    (OUTPUT / f'{PREFIX}_metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print(OUTPUT / f'{PREFIX}.pdf', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plot-only', action='store_true', help='Restyle the saved mark-count summary.')
    args = parser.parse_args()
    if args.plot_only:
        summary = read_tsv(OUTPUT / f'{PREFIX}_summary.tsv')
        metadata = json.loads((OUTPUT / f'{PREFIX}_metadata.json').read_text())
        assert metadata['source_sha256'] == sha256(SOURCE)
    else:
        per_sample, summary, metadata = collect_counts()
        write_tsv(OUTPUT / f'{PREFIX}_by_sample.tsv', per_sample)
        write_tsv(OUTPUT / f'{PREFIX}_summary.tsv', summary)
        (OUTPUT / f'{PREFIX}_metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    plot(summary, metadata)


if __name__ == '__main__':
    main()
