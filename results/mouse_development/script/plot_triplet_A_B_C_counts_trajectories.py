#!/usr/bin/env python3
"""A/B/C relative-count trajectories, adapted from plot_triplet_A_C_D_counts_trajectories.py."""
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/mouse-development-mpl')
os.environ.setdefault('XDG_CACHE_HOME', '/tmp/mouse-development-cache')
import argparse
import csv
import gzip
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
OUTPUT = BASE / 'synergistic_triplets'
STAGES = ['zygote', '2cell', '4cell', '8cell', 'morula', 'ICM', 'TE',
          'E6.5', 'E7.0', 'E7.5', 'E11.5', 'E12.5', 'E13.5', 'E14.5', 'E15.5', 'E16.5', 'P0']
COLORS = {'A': '#ca2a28', 'B': '#349939', 'C': '#ef7c21'}
PANELS = [
    ('exon / SINE', [('exon_purified_bins_AND_SINE', '', '-')]),
    ('proximal enh. CTCF− / SINE', [('proximal_not_ctcf_bound_enhancer_purified_bins_AND_SINE', '', '-')]),
    ('proximal enh. CTCF− / CpG island', [('proximal_not_ctcf_bound_enhancer_purified_bins_AND_CpG_island', '', '-')]),
]
DISTAL = {
    'positive': ('distal_ctcf_bound_enhancer_purified_bins_AND_LINE', 'CTCF+', '-'),
    'negative': ('distal_not_ctcf_bound_enhancer_purified_bins_AND_LINE', 'CTCF−', '--'),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--distal-ctcf', choices=['both', 'positive', 'negative'], default='negative')
    parser.add_argument('--prepare-only', action='store_true')
    args = parser.parse_args()
    selected = [DISTAL[k] for k in ('positive', 'negative')] if args.distal_ctcf == 'both' else [DISTAL[args.distal_ctcf]]
    distal_title = 'distal enh. / LINE' if args.distal_ctcf == 'both' else 'distal enh. CTCF' + ('+' if args.distal_ctcf == 'positive' else '−') + ' / LINE'
    if args.distal_ctcf != 'both':
        selected = [(selected[0][0], selected[0][1], '-')]
    panels = PANELS + [(distal_title, selected)]
    blocks = {entry[0] for _, entries in panels for entry in entries}
    grouped = defaultdict(list)
    with gzip.open(OUTPUT / 'all_triplets_fdr_10000.tsv.gz', 'rt') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            assert int(row['order']) == 3 and int(row['n_perm']) == 10000
            if row['block'] in blocks:
                grouped[row['sample'], row['block'], row['triplet_class']].append(row)
    assert len({key[0] for key in grouped}) == 82
    per_sample = []
    maxima = defaultdict(set)
    for (sample, block, cls), rows in grouped.items():
        denominator = len(rows)
        maxima[cls].add(denominator)
        count = sum(float(r['observed_hoi']) < 0 and float(r['q_down_sample_block_bh']) <= .05 for r in rows)
        per_sample.append(dict(sample=sample, stage=sample.split('_')[0], block=block,
                               triplet_class=cls, count=count, theoretical_max=denominator,
                               relative=count/denominator))
    assert {c: sorted(v) for c, v in maxima.items()} == {'A': [9], 'B': [16], 'C': [10]}
    stage_groups = defaultdict(list)
    for row in per_sample:
        stage_groups[row['triplet_class'], row['block'], row['stage']].append(row)
    summaries = []
    lookup = {}
    # Distinct positions for all seven preimplantation stages, followed by gastrulation.
    positions = {s: i for i, s in enumerate(STAGES[:10])}
    positions.update({s: 10 + 2*((19.5 if s == 'P0' else float(s[1:]))-11.5) for s in STAGES[10:]})
    for (cls, block, stage), rows in stage_groups.items():
        relative = [r['relative'] for r in rows]
        n = len(relative)
        std = statistics.pstdev(relative)
        record = dict(triplet_class=cls, block=block, stage=stage, x_pos=positions[stage],
                      mean_count=statistics.mean(r['count'] for r in rows),
                      theoretical_max=rows[0]['theoretical_max'], mean_relative=statistics.mean(relative),
                      std_relative=std, sem_relative=std/math.sqrt(n) if n > 1 else 0., n=n)
        summaries.append(record)
        lookup[cls, block, stage] = record
    summaries.sort(key=lambda r: (r['triplet_class'], r['block'], STAGES.index(r['stage'])))
    path = OUTPUT / 'triplet_A_B_C_relative_counts_summary.tsv'
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summaries[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(summaries)
    metadata = dict(n_permutations=10000, n_samples=82, n_stages=17, distal_ctcf=args.distal_ctcf,
                    panels=[dict(title=title, blocks=[e[0] for e in entries]) for title, entries in panels],
                    theoretical_max={c: next(iter(v)) for c, v in maxima.items()},
                    selection='observed_hoi < 0 and q_down_sample_block_bh <= .05; BH35 per sample x block',
                    aggregation='Mean of count / class theoretical maximum across tissues per stage; population SD / sqrt(n) SEM',
                    distal_handling='CTCF categories kept separate; no pooling of loci or permutation tests')
    (OUTPUT / 'triplet_A_B_C_trajectories_metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
    print(json.dumps(metadata, indent=2), flush=True)
    if args.prepare_only:
        return
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    arial = Path.home() / '.local/share/fonts/arial.ttf'
    if arial.exists():
        font_manager.fontManager.addfont(str(arial))
        plt.rcParams['font.family'] = 'Arial'
    plt.rcParams.update({'font.size': 7, 'pdf.fonttype': 42, 'ps.fonttype': 42})
    fig, axes = plt.subplots(3, 4, figsize=(25/2.54, 16/2.54), sharex=True, sharey=True)
    fig.subplots_adjust(left=.075, right=.99, top=.87, bottom=.15, wspace=.10, hspace=.19)
    for row, cls in enumerate('ABC'):
        color = COLORS[cls]
        for col, (title, entries) in enumerate(panels):
            ax = axes[row, col]
            for block, label, linestyle in entries:
                values = [lookup[cls, block, stage] for stage in STAGES]
                xs = [v['x_pos'] for v in values]
                ax.errorbar(xs, [v['mean_relative'] for v in values],
                            yerr=[v['sem_relative'] for v in values], color=color,
                            linestyle=linestyle, marker='o', markersize=2.4,
                            markerfacecolor='white' if linestyle == '--' else color,
                            markeredgewidth=.4, linewidth=.8, elinewidth=.5,
                            capsize=1.3, label=label)
            ax.set_ylim(-.04, 1.06)
            ax.set_yticks([0, .25, .5, .75, 1])
            ax.set_xticks([positions[s] for s in STAGES], labels=STAGES, rotation=60, ha='right', fontsize=5.3)
            ax.set_xlim(-.6, positions['P0']+.6)
            ax.tick_params(length=2, pad=1, labelsize=6)
            ax.axhline(0, color='gray', linewidth=.4, alpha=.3)
            for split in (6.5, 9.5):
                ax.axvline(split, color='gray', linestyle='--', linewidth=.5, alpha=.5)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            if row == 0:
                ax.set_title(title, fontsize=7.5, pad=7)
                if col == 3 and args.distal_ctcf == 'both':
                    ax.legend(frameon=False, fontsize=6, loc='upper left')
            if col == 0:
                ax.set_ylabel('Triplet '+cls+'\noccurrence proportion', color=color, fontsize=7)
                ax.tick_params(axis='y', colors=color)
            if row == 2:
                ax.set_xlabel('Stage', fontsize=7)
    fig.suptitle('Significant synergistic triplets: relative-count trajectories\n'
                 '10,000 permutations; BH35 q ≤ 0.05 and O-info < 0; mean ± SEM across tissues\n'
                 'Seven-signal panel: A / 9, B / 16, C / 10', fontsize=9, y=.98)
    for extension in ('pdf', 'png'):
        path = OUTPUT / ('triplet_A_B_C_counts_trajectories_n5.'+extension)
        fig.savefig(path, dpi=300, bbox_inches='tight', pad_inches=.03)
        print('Saved', path, flush=True)
    plt.close(fig)


if __name__ == '__main__':
    main()
