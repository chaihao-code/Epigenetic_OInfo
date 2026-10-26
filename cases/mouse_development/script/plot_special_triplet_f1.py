#!/usr/bin/env python3
"""F1 heatmap for H3K9me3 + H3K27me3 + a non-methylation signal, within B."""
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/mouse-development-mpl')
os.environ.setdefault('XDG_CACHE_HOME', '/tmp/mouse-development-cache')
from pathlib import Path
import csv
import gzip
import hashlib
import json
from collections import defaultdict
from datetime import datetime, timezone

BASE = Path(__file__).resolve().parents[1]
ROOT = BASE / 'result_oinfo_n_bins_5_triplet_perm10000'
OUTPUT = BASE / 'synergistic_triplets'
SOURCE = OUTPUT / 'all_triplets_fdr_10000.tsv.gz'
STAGES = ['zygote', '2cell', '4cell', '8cell', 'morula', 'ICM', 'TE',
          'E6.5', 'E7.0', 'E7.5', 'E11.5', 'E12.5', 'E13.5', 'E14.5', 'E15.5', 'E16.5', 'P0']
REGIONS = ['promoter', 'proximal_ctcf_bound_enhancer_purified_bins',
           'proximal_not_ctcf_bound_enhancer_purified_bins',
           'distal_ctcf_bound_enhancer_purified_bins',
           'distal_not_ctcf_bound_enhancer_purified_bins',
           'exon_purified_bins', 'intron_purified_bins', 'None']
LABELS = ['promoter', 'proximal\nCTCF+', 'proximal\nCTCF−',
          'distal\nCTCF+', 'distal\nCTCF−', 'exon', 'intron', 'None']
FEATURES = ['CpG_island', 'LINE', 'SINE', 'LTR', 'None']
BLOCKS = [r + '_AND_' + f for r in REGIONS for f in FEATURES]


def bh(p):
    q = [1.] * len(p)
    current = 1.
    for rank, i in reversed(list(enumerate(sorted(range(len(p)), key=p.__getitem__), 1))):
        current = min(current, p[i] * len(p) / rank)
        q[i] = current
    return q


def table(path, rows):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def main():
    OUTPUT.mkdir(exist_ok=True)
    manifest = list(csv.DictReader((ROOT / 'input_manifest.tsv').open(), delimiter='\t'))
    results, inputs, missing = [], [], []
    with gzip.open(SOURCE, 'rt') as handle:
        corrected = list(csv.DictReader(handle, delimiter='\t'))
    sample_groups = defaultdict(lambda: defaultdict(list))
    for row in corrected:
        assert int(row['order']) == 3 and int(row['n_perm']) == 10000
        sample_groups[row['sample']][row['block']].append(row)
    table_hash = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    for entry in manifest:
        sample = entry['sample']
        groups = sample_groups[sample]
        assert set(groups) == set(BLOCKS), sample
        for block in BLOCKS:
            rows = groups[block]
            assert len(rows) == 35
            tp = fp = fn = tn = 0
            for row in rows:
                names = row['triplet']
                if any(x in names for x in ('WGBS', 'CpG_methyl')):
                    continue
                if not any(x in names for x in ('H3K9me3', 'H3K27me3', 'H3K36me3')):
                    continue
                predicted = 'H3K9me3' in names and 'H3K27me3' in names
                positive = float(row['observed_hoi']) < 0 and float(row['q_down_sample_block_bh']) <= .05
                tp += predicted and positive
                fp += predicted and not positive
                fn += not predicted and positive
                tn += not predicted and not positive
            assert tp + fp == 4
            results.append(dict(sample=sample, stage=sample.split('_')[0], block=block,
                                TP=tp, FP=fp, FN=fn, TN=tn,
                                precision=tp / (tp + fp), recall=tp / (tp + fn) if tp + fn else 0.,
                                F1=2 * tp / (2 * tp + fp + fn)))
        inputs.append(dict(sample=sample, source=str(SOURCE), sha256=table_hash))
    assert results
    grouped = defaultdict(list)
    for row in results:
        grouped[row['stage'], row['block']].append(row)
    stages = [s for s in STAGES if any(k[0] == s for k in grouped)]
    averages = []
    for stage in stages:
        for block in BLOCKS:
            rows = grouped[stage, block]
            assert rows
            averages.append(dict(stage=stage, block=block, n_samples=len(rows),
                                 mean_F1=sum(x['F1'] for x in rows) / len(rows)))
    table(OUTPUT / 'f1_by_sample_block.tsv', results)
    table(OUTPUT / 'f1_by_stage_block.tsv', averages)
    table(OUTPUT / 'special_triplet_f1_inputs.tsv', inputs)
    metadata = dict(created_at_utc=datetime.now(timezone.utc).isoformat(), n_samples=len(inputs),
                    missing_samples=missing, n_permutations=10000,
                    prediction='H3K9me3 + H3K27me3 + X; X = ATAC,H3K27ac,H3K4me1,H3K4me3',
                    universe='Triplet B only; no DNA methylation and at least one repressive mark',
                    reference='observed_hoi < 0 and q_down <= .05; BH of all 35 triplets per sample x block',
                    aggregation='Arithmetic mean of per-sample F1 across available tissues within stage; not pooled F1',
                    formula='F1=2TP/(2TP+FP+FN); four predicted positives per sample/block')
    (OUTPUT / 'special_triplet_f1_metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    arial = Path.home() / '.local/share/fonts/arial.ttf'
    if arial.exists():
        font_manager.fontManager.addfont(str(arial))
        plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['pdf.fonttype'] = 42
    values = {(x['stage'], x['block']): x['mean_F1'] for x in averages}
    matrix = np.array([[values[s, b] for b in BLOCKS] for s in stages])
    fig, ax = plt.subplots(figsize=(17, max(6, .36 * len(stages) + 2.5)))
    fig.subplots_adjust(left=.12, right=.94, top=.84, bottom=.24)
    im = ax.imshow(matrix, cmap='Reds', vmin=0, vmax=1, aspect='auto', interpolation='nearest')
    counts = {s: len(grouped[s, BLOCKS[0]]) for s in stages}
    ax.set_yticks(range(len(stages)), labels=[f'{s} (n={counts[s]})' for s in stages], fontsize=10)
    ax.set_xticks([])
    for i in range(len(stages)):
        for j in range(40):
            value = matrix[i, j]
            ax.text(j, i, f'{value:.2f}', ha='center', va='center', fontsize=6,
                    color='white' if value > .55 else 'black')
    for j in range(5, 40, 5):
        ax.axvline(j - .5, color='black', lw=.8, alpha=.4)
    for i, label in enumerate(LABELS):
        ax.annotate('', xy=(5*i+4.3, -.045), xytext=(5*i-.3, -.045),
                    xycoords=ax.get_xaxis_transform(), arrowprops=dict(arrowstyle='<->', color='gray', lw=1.5))
        ax.text(5*i+2, -.06, label, transform=ax.get_xaxis_transform(), ha='center', va='top', fontsize=10)
    for j, block in enumerate(BLOCKS):
        ax.text(j, -.20, block.split('_AND_')[1].replace('_', ' '), transform=ax.get_xaxis_transform(),
                ha='center', va='top', rotation=90, fontsize=9)
    cbar = fig.colorbar(im, ax=ax, fraction=.02, pad=.01)
    cbar.set_label('Mean F1 across tissues', fontsize=10)
    cbar.set_ticks([0, .2, .4, .6, .8, 1])
    ax.set_title('H3K9me3 + H3K27me3 + X: F1 for significant synergistic Triplet B\n'
                 'X = ATAC, H3K27ac, H3K4me1, H3K4me3; 10,000 permutations; BH within each sample × region (35 tests)\n'
                 f'{len(inputs)}/{len(manifest)} samples; synergy: O-info < 0 and q ≤ 0.05', fontsize=12, pad=15)
    for ext in ('png', 'pdf'):
        path = OUTPUT / ('H3K9me3_H3K27me3_X_F1_heatmap.' + ext)
        fig.savefig(path, dpi=300, bbox_inches='tight')
        print('Saved', path, flush=True)
    plt.close(fig)
    print(json.dumps(dict(n_samples=len(inputs), n_stages=len(stages), n_blocks=40)), flush=True)


if __name__ == '__main__':
    main()
