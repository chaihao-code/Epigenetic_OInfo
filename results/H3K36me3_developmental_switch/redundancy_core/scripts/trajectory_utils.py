"""Shared input validation, aggregation and plotting for developmental trajectories."""
import csv
import hashlib
import json
import math
import os
from collections import Counter
from pathlib import Path
from statistics import mean, stdev

os.environ.setdefault('MPLCONFIGDIR', '/tmp/bulk16-support-mpl')
os.environ.setdefault('XDG_CACHE_HOME', '/tmp/bulk16-support-cache')
ROOT = Path(__file__).resolve().parents[1]
MARKS = ['ATAC', 'H3K4me1', 'H3K4me3', 'H3K27ac', 'H3K9me3', 'H3K27me3', 'DNA_methylation', 'H3K36me3']


def read(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f, delimiter='\t' if str(path).endswith('.tsv') else ','))


def write(path, rows, fields=None):
    with Path(path).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields or list(rows[0]), delimiter='\t' if str(path).endswith('.tsv') else ',')
        w.writeheader()
        w.writerows(rows)


def config():
    return json.loads((ROOT / 'inputs/analysis_manifest.json').read_text())


def selected_tasks(manifest):
    return [t for t in manifest['tasks'] if t['period'] == 'preimplantation'
            or t['sample'].removeprefix(t['stage'] + '_') in manifest['late_plot_tissues']]


def collect_cores(manifest):
    all_cores = []
    program_hash = hashlib.sha256(Path(manifest['core_program']).read_bytes()).hexdigest()
    for task in manifest['tasks']:
        path = ROOT / 'results' / task['sample'] / 'redundancy_cores.csv'
        metadata = json.loads(path.with_suffix('.metadata.json').read_text())
        assert metadata['status'] == 'complete' and metadata['sample'] == task['sample']
        assert metadata['n_blocks'] == 40 and metadata['program_sha256'] == program_hash
        assert set(metadata['signals']) == set(manifest['signals'])
        assert metadata['z_threshold'] == 2.5 and metadata['eta_threshold'] == .5
        rows = read(path)
        assert len(rows) == metadata['n_cores']
        assert len({(r['block'], r['signals']) for r in rows}) == len(rows)
        for r in rows:
            assert float(r['hoi']) > 0 and int(r['order']) == len(r['signals'].split('/'))
            if r['promoted'] != 'True':
                assert float(r['mod_z']) > 2.5 and r['passes_modality'] == 'True'
        for block in {r['block'] for r in rows}:
            subsets = [set(r['signals'].split('/')) for r in rows if r['block'] == block]
            assert not any(a < b for a in subsets for b in subsets)
        all_cores.extend(dict(stage=task['stage'], period=task['period'], **r) for r in rows)
    fields = ['stage', 'period', 'sample', 'block', 'order', 'signals', 'hoi', 'mod_z', 'n_background',
              'n_components', 'delta_bic', 'eta_new', 'positive', 'passes_z', 'passes_modality',
              'selected_before_pruning', 'promoted']
    write(ROOT / 'results/redundancy_cores.csv', all_cores, fields)
    counts = Counter(r['sample'] for r in all_cores)
    sample_summary = [dict(sample=t['sample'], stage=t['stage'], n_blocks=40, n_cores=counts[t['sample']]) for t in manifest['tasks']]
    write(ROOT / 'results/core_counts_by_sample.tsv', sample_summary)
    return all_cores


def sample_counts(kind, manifest):
    tasks = selected_tasks(manifest)
    assert len(tasks) == 42
    values = {(t['sample'], mark): 0 for t in tasks for mark in MARKS}
    if kind == 'redundancy_core':
        rows = collect_cores(manifest)
        for r in rows:
            for mark in r['signals'].split('/'):
                key = r['sample'], mark
                if key in values:
                    values[key] += 1
    else:
        rows = read(manifest['synergistic_triplet_counts_source'])
        assert len(rows) == len(values)
        seen = set()
        for r in rows:
            sample = r['stage'] if r['stage'] in manifest['early_stages'] else r['sample']
            key = sample, r['mark']
            assert key in values and key not in seen
            seen.add(key)
            values[key] = int(r['n_significant'])
        assert seen == set(values)
    records = [dict(sample=t['sample'], stage=t['stage'], period=t['period'], mark=m,
                    count=values[t['sample'], m]) for t in tasks for m in MARKS]
    write(ROOT / 'results' / (kind + '_mark_counts_by_sample.tsv'), records)
    return records


def plot_trajectories(kind):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    manifest = config()
    records = sample_counts(kind, manifest)
    early = ['zygote', '2cell', '4cell', '8cell', 'morula', 'E3.5']
    stages = early + manifest['late_stages']
    positions = dict(zip(early, range(6)))
    positions.update({s: 7 + 2 * ((19.5 if s == 'P0' else float(s[1:])) - 11.5) for s in manifest['late_stages']})
    summary = []
    for stage in stages:
        for mark in MARKS:
            z = [r['count'] for r in records if r['mark'] == mark
                 and (r['stage'] in ['ICM', 'TE'] if stage == 'E3.5' else r['stage'] == stage)]
            expected = 2 if stage == 'E3.5' else (1 if stage in early else 5)
            assert len(z) == expected
            sd = stdev(z) if len(z) > 1 else 0.
            summary.append(dict(stage=stage, mark=mark, x_pos=positions[stage], mean_count=mean(z),
                                std_count=sd, sem_count=sd / math.sqrt(len(z)), n=len(z)))
    font_path = font_manager.findfont('Arial', fallback_to_default=False)
    plt.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'axes.labelsize': 8,
                         'axes.titlesize': 10, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
                         'legend.fontsize': 8, 'pdf.fonttype': 42, 'ps.fonttype': 42})
    fig, ax = plt.subplots(figsize=(18 / 2.54, 6 / 2.54))
    handles = []
    for mark in MARKS:
        highlight = mark == 'H3K36me3'
        color = '#000000' if highlight else manifest['palette'][mark]
        for i, period in enumerate([stages]):
            z = [next(r for r in summary if r['stage'] == s and r['mark'] == mark) for s in period]
            line = ax.errorbar([r['x_pos'] for r in z], [r['mean_count'] for r in z],
                              yerr=[r['sem_count'] for r in z], color=color, marker='o',
                              markersize=4.2 if highlight else 2.8,
                              markeredgecolor='#000000' if highlight else '#606060', markeredgewidth=.3,
                              linewidth=2.6 if highlight else .9,
                              elinewidth=1.1 if highlight else .5, capsize=1.5,
                              zorder=10 if highlight else 3,
                              label=mark.replace('_', ' ') if i == 0 else None)
            if i == 0:
                handles.append(line)
    ax.set_xticks([positions[s] for s in stages], ['Zygote', '2-cell', '4-cell', '8-cell', 'Morula', 'E3.5'] + manifest['late_stages'], rotation=45, ha='right')
    ax.axvline(6, color='#999999', linestyle='--', linewidth=.6)
    ax.set_ylim(bottom=0)
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(width=.5, length=2)
    label = 'Redundancy cores' if kind == 'redundancy_core' else 'Synergistic triplets'
    ax.set_title(label + ' containing each signal')
    ax.set_ylabel('Mean count per sample')
    fig.legend(handles, [m.replace('_', ' ') for m in MARKS], loc='lower center',
               bbox_to_anchor=(.5, .005), ncol=4, frameon=False, handlelength=1.4,
               columnspacing=1., handletextpad=.4, labelspacing=.35)
    fig.subplots_adjust(left=.09, right=.99, top=.85, bottom=.36)
    for ext in ['pdf', 'png']:
        fig.savefig(ROOT / 'plots' / (kind + '_mark_trajectories.' + ext), dpi=300)
    plt.close(fig)
    write(ROOT / 'results' / (kind + '_mark_trajectories.tsv'), summary)
    (ROOT / 'results' / (kind + '_plot_validation.json')).write_text(json.dumps(dict(
        status='complete', kind=kind, n_samples=42, n_stage_points=13,
        ICM_TE='Combined under E3.5 (mean of both lineages)',
        late_tissues=manifest['late_plot_tissues'], zero_samples_included=True,
        error_bars='SEM across samples/tissues; sample SD ddof=1; single-sample SEM=0',
        font_family='Arial', font_file=font_path, font_size_pt=8,
        title_font_size_pt=10, figure_size_cm=[18, 6],
        connect_E3p5_to_E11p5=True,
        H3K36me3_style=dict(color='#000000', linewidth_pt=2.6, markersize_pt=4.2)), indent=2) + '\n')
    print('Saved ' + kind + ' trajectory figure and tables.', flush=True)
