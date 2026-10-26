#!/usr/bin/env python3
"""Rebuild the retained descriptive C-triplet recurrence table (no new tests)."""
import argparse
from pathlib import Path

import pandas as pd

PROJECT = Path(__file__).resolve().parents[1]


def recurrence_table():
    data = pd.read_csv(PROJECT / 'ABC_triplet_classes/all_triplet_classes_and_tests.tsv', sep='\t')
    labels = pd.read_csv(PROJECT / 'ABC_triplet_classes/sample_labels.tsv', sep='\t')
    data = data[data['class'].eq('C')].merge(labels[['sample', 'sample_id', 'pooled_group']],
                                           on='sample', validate='many_to_one')
    data['nominal'] = data.observed_hoi.lt(0) & data.p_down.lt(.05)
    data['fdr'] = data.observed_hoi.lt(0) & data.q_down_global_BH.lt(.05)
    keys = ['block', 'triplet']
    table = None
    for group in ['R', 'S', 'pB']:
        sub = data[data.pooled_group.eq(group)]
        grouped = sub.groupby(keys).agg(n_samples=('sample', 'nunique'),
            n_synergy_nominal=('nominal', 'sum'), n_synergy_FDR005=('fdr', 'sum'),
            mean_omega=('observed_hoi', 'mean')).add_prefix(group + '_')
        table = grouped if table is None else table.join(grouped, validate='one_to_one')
    for sample in sorted(data.loc[data.pooled_group.eq('S'), 'sample_id'].unique()):
        sub = data[data.sample_id.eq(sample)].set_index(keys)
        table = table.join(sub[['observed_hoi', 'p_down', 'q_down_global_BH']].add_prefix(sample + '_'),
                           validate='one_to_one')
    table = table.reset_index().sort_values(
        ['S_n_synergy_nominal', 'S_n_synergy_FDR005', 'block', 'triplet'],
        ascending=[False, False, True, True])
    assert len(table) == 160 and table.S_n_samples.eq(5).all()
    return table, data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=PROJECT / 'ABC_triplet_classes/C_pattern')
    args = parser.parse_args()
    table, data = recurrence_table()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output_dir / 'C_triplet_recurrence_by_region.tsv', sep='\t', index=False)
    sub = data[data.pooled_group.eq('S')]
    lines = ['# C 类在 S 样本中的描述性模式', '',
        '名义协同：observed_hoi<0 且 p_down<0.05。FDR 协同：observed_hoi<0 且 q_down_global_BH<0.05。',
        '沿用全部 17 个样本、23,800 项检验的全局下尾 BH 校正；此处不进行组间差异或富集检验。', '',
        '| 三信号组合 | S 名义显著 / 200 | S FDR显著 / 200 |', '|---|---:|---:|']
    counts = sub.groupby('triplet')[['nominal', 'fdr']].sum().sort_values('nominal', ascending=False)
    for triplet, row in counts.iterrows():
        lines.append(f'| {triplet} | {int(row.nominal)} | {int(row.fdr)} |')
    lines += ['', f'S 组名义协同 {int(sub.nominal.sum())} 次，FDR 协同 {int(sub.fdr.sum())} 次。', '']
    for name, column in [('名义', 'S_n_synergy_nominal'), ('FDR', 'S_n_synergy_FDR005')]:
        counts = [int(table[column].ge(n).sum()) for n in [1, 4, 5]]
        lines.append(f'{name}：至少 1/5、4/5、5/5 样本复现的区域×组合数量分别为 {counts}。')
    lines += ['', '各区域并非独立生物重复；本表仅描述跨供者复现情况，不能据此证明生物机制。',
              '完整逐样本效应、P 值和 FDR 见 C_triplet_recurrence_by_region.tsv。', '']
    (args.output_dir / 'summary_zh.md').write_text('\n'.join(lines))


if __name__ == '__main__':
    main()
