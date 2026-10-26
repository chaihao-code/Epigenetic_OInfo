#!/usr/bin/env python3
"""
Generic script to analyze local O-information on genomic regions defined by feature combinations.

通用脚本：对基因组特征组合定义的目标区域，批量计算 local O-information 和年龄回归。

Usage:
  python3 run_analysis_generic.py --name <analysis_name> --features <feat1> <feat2> ... --op AND|OR [--exclude feat3 feat4 ...] [--out-dir <custom_output_dir>]

计算流程
========
1. compute_target_bins()    → 用 AND/OR/NOT 组合 feature masks 得到目标 bins
2. extract_signals()        → 从离散 BED 文件读取每个 sample 的 7 个 mark 状态
3. oi_local(X)              → 对每组 7 个状态计算 local O-info
4. run_regression()         → 对每个 bin，O-info ~ age 做 Pearson r + 线性回归
5. BH FDR correction        → qvalue_bh, sig_fdr_0.05, sig_fdr_0.10

输入数据
--------
- 离散信号: data/aging_thyroid_input_normed_binned_discreted/
  7 marks × 9 samples，每文件是 1kb bin 的离散状态（0-4）
- 基因组特征 masks: result/genome_features/hg38/
  LINE, SINE, LTR, CpG_island, proximal_enhancer 等二进制 mask

输出文件
--------
- target_bins.bed              → 目标区域坐标
- extracted_signals.tsv        → 每个 bin 的 7 个 mark 离散状态
- oinfo_per_sample.bed         → 原始 O-info 值（9 samples × N bins）+ 回归统计
- oinfo_per_sample_fdr.bed     → 加上 qvalue, sig_fdr_0.05, sig_fdr_0.10
"""

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

import sys
import os
from pathlib import Path
import os
import glob
import argparse
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

# O-info core: hoi_code exact counting (compute_hoi_oinfo.py) + permutation null
from oi_local_hoi import oi_local_hoi, permutation_null

# ===================== Config =====================
# Machine-specific paths (cluster)
# 元数据路径：包含 sample_id 和 age 的对应关系
META_PATH = _release_path('@ROOT@/cases/aging/thyroid/metadata/thyroid_healthy_metadata.tsv')

# 离散信号目录：7 个 marks × 9 个 samples 的离散 BED 文件
DISCRETE_DIR = _release_path('@ROOT@/data/aging/thyroid/input_normed_binned_discreted')

# 基因组特征 mask 目录：二进制 mask 文件（含 *_purified_bins.tsv 与 None_func/None_rep）
FEATURE_DIR = _release_path('@ROOT@/data/genome_features/hg38')

# 输出基础目录（可通过 --out-dir 参数覆盖）
OUT_BASE = _release_path('@ROOT@/cases/aging/thyroid/outputs/local_oinfo_age_correlation')

# 7 个表观遗传标记
MARKS = ['H3K4me1', 'H3K4me3', 'H3K9me3', 'H3K27me3', 'H3K36me3', 'H3K27ac', 'WGB-Seq']

# 只分析常染色体（chr1-chr22），排除性染色体和线粒体
MAIN_CHRS = [f'chr{i}' for i in range(1, 23)]

# 排除样本
EXCLUDED_SAMPLES = {'CEMT0086', 'CEMT0088', 'CEMT0090'}

# 排除样本
EXCLUDED_SAMPLES = {'CEMT0086', 'CEMT0088', 'CEMT0090'}


def get_discrete_file(sample, mark):
    """
    根据 sample ID 和 mark 名称，找到对应的离散信号文件路径。

    文件名模式：
    - WGB-Seq: *.{sample}.WGB-Seq.methylation_profile.bin1000.step1000.discrete.bed
    - 其他 mark: *.{sample}.{mark}.signal_unstranded.bin1000.step1000.divInput.discrete.bed
    """
    if mark == 'WGB-Seq':
        pattern = os.path.join(DISCRETE_DIR, f'*.{sample}.WGB-Seq.methylation_profile.bin1000.step1000.discrete.bed')
    else:
        pattern = os.path.join(DISCRETE_DIR, f'*.{sample}.{mark}.signal_unstranded.bin1000.step1000.divInput.discrete.bed')
    files = glob.glob(pattern)
    if not files:
        raise FileNotFoundError(f"Missing file for {sample} {mark}: {pattern}")
    return files[0]


def load_feature(feature_name):
    """
    加载基因组特征 mask 文件。

    文件格式：TSV，两列
      bin_id    value
      chr1:0-1000    1
      chr1:1000-2000  0
      ...

    返回 DataFrame，包含 bin_id, value, chr, start, end 列。
    只保留常染色体（chr1-chr22）。
    """
    path = os.path.join(FEATURE_DIR, f'{feature_name}_bins.tsv')
    df = pd.read_csv(path, sep='\t')
    # 从 bin_id (如 "chr1:0-1000") 解析染色体和坐标
    df['chr'] = df['bin_id'].str.split(':').str[0]
    df['start'] = df['bin_id'].str.split(':').str[1].str.split('-').str[0].astype(int)
    df['end'] = df['bin_id'].str.split(':').str[1].str.split('-').str[1].astype(int)
    df = df[df['chr'].isin(MAIN_CHRS)].copy()
    return df


def compute_target_bins(features, op, exclude_features=None):
    """
    根据特征组合计算目标 bin ID 列表。

    参数
    ----
    features : list[str]
        要组合的特征名称列表。
    op : str
        'AND' 或 'OR'。AND 表示交集，OR 表示并集。
    exclude_features : list[str] or None
        要排除的特征名称列表（NOT 操作）。

    逻辑
    ----
    1. 加载第一个特征作为基础
    2. 对其他特征应用 op（AND/OR）
    3. 对 exclude_features 应用 NOT（排除这些特征为 1 的 bins）
    4. 返回最终 value == 1 的 bin_id 列表
    """
    # 加载第一个特征
    base_df = load_feature(features[0])
    base_col = base_df.columns[1]  # 第二列是 mask 值（0 或 1）
    result = base_df[['bin_id', base_col]].copy()
    result.columns = ['bin_id', 'val']

    if op == 'AND':
        # AND: 只保留所有特征都为 1 的 bins
        for feat in features[1:]:
            df = load_feature(feat)
            col = df.columns[1]
            df = df[['bin_id', col]].copy()
            df.columns = ['bin_id', 'val2']
            result = result.merge(df, on='bin_id', how='inner')
            result['val'] = result['val'] & result['val2']
            result = result[['bin_id', 'val']]
    elif op == 'OR':
        # OR: 保留任一特征为 1 的 bins
        for feat in features[1:]:
            df = load_feature(feat)
            col = df.columns[1]
            df = df[['bin_id', col]].copy()
            df.columns = ['bin_id', 'val2']
            result = result.merge(df, on='bin_id', how='outer')
            result['val'] = result['val'].fillna(0).astype(int) | result['val2'].fillna(0).astype(int)
            result = result[['bin_id', 'val']]

    # 应用 exclusions（NOT）
    if exclude_features:
        for feat in exclude_features:
            df = load_feature(feat)
            col = df.columns[1]
            df = df[['bin_id', col]].copy()
            df.columns = ['bin_id', 'val_exc']
            result = result.merge(df, on='bin_id', how='left')
            result['val_exc'] = result['val_exc'].fillna(0).astype(int)
            # val & (1 - val_exc): 如果排除特征为 1，则结果为 0
            result['val'] = result['val'] & (1 - result['val_exc'])
            result = result[['bin_id', 'val']]

    target_bins = result[result['val'] == 1]['bin_id'].tolist()
    return target_bins


def run_analysis(analysis_name, features, op='AND', exclude_features=None,
                 n_perm=10000, perm_seed=1234):
    """
    主分析函数：计算指定区域的 local O-information 和年龄回归。

    完整流程：
    1. 加载元数据（sample -> age 映射）
    2. 计算目标 bins（特征组合）
    3. 构建坐标索引（从参考文件建立 bin_id -> 行索引 映射）
    4. 提取信号（从离散 BED 文件读取目标 bins 的 7-mark 状态）
    5. 计算 local O-information（每个 sample 独立计算）
    6. 线性回归（每个 bin 的 O-info ~ age）
    7. FDR 校正（Benjamini-Hochberg）
    8. 输出结果文件
    """
    out_dir = os.path.join(OUT_BASE, analysis_name)
    os.makedirs(out_dir, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"Analysis: {analysis_name}")
    print(f"Features: {features}, op={op}, exclude={exclude_features}")
    print(f"{'='*60}")

    # ===================== 1. 加载元数据 =====================
    meta = pd.read_csv(META_PATH, sep='\t')
    # 按 sample_id 去重，获取每个 sample 的年龄
    age_map = meta.groupby('sample_id')['age'].first().to_dict()
    # 按年龄排序 samples
    samples = sorted([s for s in age_map.keys() if s not in EXCLUDED_SAMPLES],
                     key=lambda s: age_map[s])
    print(f"Samples: {samples}")
    print(f"Ages: {[age_map[s] for s in samples]}")

    # ===================== 2. 计算目标 bins =====================
    print("\n=== Computing target bins ===")
    target_bin_ids = compute_target_bins(features, op, exclude_features)
    print(f"Target bins: {len(target_bin_ids)}")
    if len(target_bin_ids) == 0:
        print("No target bins found! Skipping.")
        return

    # 解析 bin_id 为坐标（chr, start, end）
    target_df = pd.DataFrame({
        'bin_id': target_bin_ids,
        'chr': [b.split(':')[0] for b in target_bin_ids],
        'start': [int(b.split(':')[1].split('-')[0]) for b in target_bin_ids],
        'end': [int(b.split(':')[1].split('-')[1]) for b in target_bin_ids],
    })
    target_df = target_df.sort_values(['chr', 'start']).reset_index(drop=True)
    target_df[['chr', 'start', 'end']].to_csv(os.path.join(out_dir, 'target_bins.bed'), sep='\t', index=False, header=False)
    print(f"Saved target bins to {out_dir}/target_bins.bed")

    # ===================== 3. 提取信号（按 bin_id 对齐）=====================
    print("\n=== Extracting signals (aligned by bin_id) ===")
    signal_dict = {}  # (sample, mark) -> numpy array of target bin signals
    target_bin_ids = target_df['bin_id'].values
    for s in samples:
        for m in MARKS:
            fpath = get_discrete_file(s, m)
            # 读取完整文件（包括坐标和值），确保按 bin_id 对齐
            df_sig = pd.read_csv(fpath, sep='\t', header=None, names=['chr', 'start', 'end', 'val'])
            df_sig = df_sig[df_sig['chr'].isin(MAIN_CHRS)].copy()
            df_sig['bin_id'] = df_sig['chr'] + ':' + df_sig['start'].astype(str) + '-' + df_sig['end'].astype(str)
            # 按 bin_id 对齐提取目标 bins
            df_sig = df_sig.set_index('bin_id')
            target_vals = df_sig.loc[target_bin_ids, 'val'].values.astype(np.int8)
            # 检查缺失
            n_missing = pd.isna(target_vals).sum()
            if n_missing > 0:
                print(f"  Warning: {n_missing} missing bins for {s} {m}")
            signal_dict[(s, m)] = target_vals

    # 构建信号 DataFrame：每行是一个 bin，每列是一个 sample × mark 的组合
    signal_df = target_df[['chr', 'start', 'end', 'bin_id']].copy()
    for s in samples:
        for m in MARKS:
            signal_df[f'{s}_{m}'] = signal_dict[(s, m)]
    signal_df.to_csv(os.path.join(out_dir, 'extracted_signals.tsv'), sep='\t', index=False)

    # ===================== 5. 计算 local O-information =====================
    print("\n=== Computing local O-information ===")
    oinfo_dict = {}  # sample -> o_local array
    for s in samples:
        # 提取该 sample 的 7 个 mark 的列
        cols = [f'{s}_{m}' for m in MARKS]
        X = signal_df[cols].values.astype(np.int8)
        # 计算该 sample 的 local O-information（hoi_code 精确计数实现）
        # 注意：每个 sample 独立计算，即每个 sample 有自己的概率分布
        o_local = oi_local_hoi(X)
        oinfo_dict[s] = o_local
        print(f"  {s} ({age_map[s]:.0f}y): mean={o_local.mean():.4f}")

    # ===================== 6. 线性回归 =====================
    print("\n=== Running linear regression ===")
    ages_arr = np.array([age_map[s] for s in samples])
    regression_results = []
    n_bins = len(target_df)
    for i in range(n_bins):
        # 收集该 bin 在所有 samples 中的 O-info 值
        oi_vals = np.array([oinfo_dict[s][i] for s in samples])
        # Pearson 相关系数
        r, p = stats.pearsonr(ages_arr, oi_vals)
        # 线性回归：斜率、截距、标准误
        slope, intercept, r_value, p_value, std_err = stats.linregress(ages_arr, oi_vals)
        regression_results.append({
            'r': r, 'p': p, 'slope': slope, 'intercept': intercept,
            'std_err': std_err, 'mean_oinfo': oi_vals.mean(), 'std_oinfo': oi_vals.std(ddof=1),
        })
    reg_df = pd.DataFrame(regression_results)
    print(f"Significant (p < 0.05): {(reg_df['p'] < 0.05).sum()}")

    # ===================== 7. 构建输出 =====================
    out_df = target_df[['chr', 'start', 'end', 'bin_id']].copy()
    for s in samples:
        out_df[f'{s}_{age_map[s]:.0f}y'] = oinfo_dict[s]
    out_df['r'] = reg_df['r'].values
    out_df['p'] = reg_df['p'].values
    out_df['slope'] = reg_df['slope'].values
    out_df['intercept'] = reg_df['intercept'].values
    out_df['std_err'] = reg_df['std_err'].values
    out_df['mean_oinfo'] = reg_df['mean_oinfo'].values
    out_df['std_oinfo'] = reg_df['std_oinfo'].values

    out_bed = os.path.join(out_dir, 'oinfo_per_sample.bed')
    out_df.to_csv(out_bed, sep='\t', index=False)
    print(f"Saved: {out_bed}")

    # ===================== 8. FDR 校正 =====================
    print("\n=== FDR correction ===")
    mask = out_df['p'].notna()
    p_values = out_df.loc[mask, 'p'].values
    # Benjamini-Hochberg FDR 校正
    reject_05, q_values, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
    reject_10, _, _, _ = multipletests(p_values, alpha=0.10, method='fdr_bh')
    out_df['qvalue_bh'] = np.nan
    out_df.loc[mask, 'qvalue_bh'] = q_values
    out_df['sig_fdr_0.05'] = False
    out_df['sig_fdr_0.10'] = False
    out_df.loc[mask, 'sig_fdr_0.05'] = reject_05
    out_df.loc[mask, 'sig_fdr_0.10'] = reject_10

    print(f"FDR < 0.05: {reject_05.sum()}")
    print(f"FDR < 0.10: {reject_10.sum()}")
    print(f"Bonferroni: {(out_df['p'] < 0.05/len(out_df)).sum()}")

    out_df.to_csv(os.path.join(out_dir, 'oinfo_per_sample_fdr.bed'), sep='\t', index=False)

    # ===================== 8b. 置换零检验（每 sample 的 global O-info）=====================
    # 零假设：7 个 mark 在 region 内相互独立（逐 mark 跨 bin 置换）
    print("\n=== Permutation null (global O-information per sample) ===")
    perm_rows = []
    for s in samples:
        cols = [f'{s}_{m}' for m in MARKS]
        X = signal_df[cols].values.astype(np.int8)
        res = permutation_null(X, n_perm=n_perm, n_states=5,
                                seed=perm_seed)
        perm_rows.append({'sample': s, 'age': age_map[s], **res})
        print(f"  {s}: hoi={res['observed_hoi']:.4f}, null_mean={res['null_mean']:.4f}, "
              f"p_dev={res['p_dev']:.4g}")
    perm_df = pd.DataFrame(perm_rows)
    perm_df.to_csv(os.path.join(out_dir, 'region_oinfo_summary.tsv'),
                   sep='\t', index=False)

    # Summary
    print("\n=== Summary ===")
    print(f"Negative r: {(out_df['r'] < 0).sum()}")
    print(f"Positive r: {(out_df['r'] > 0).sum()}")
    for s in samples:
        col = f'{s}_{age_map[s]:.0f}y'
        vals = out_df[col]
        print(f"{col}: mean={vals.mean():.4f}, pos={(vals > 0).mean():.2%}, neg={(vals < 0).mean():.2%}")

    print(f"\nDone: {analysis_name}\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--name', required=True, help='Analysis name (used as output subdirectory)')
    parser.add_argument('--features', nargs='+', required=True, help='Feature names to combine')
    parser.add_argument('--op', default='AND', choices=['AND', 'OR'], help='Set operation: AND or OR')
    parser.add_argument('--exclude', nargs='*', default=None, help='Features to exclude (NOT operation)')
    parser.add_argument('--out-dir', default=None, help='Custom output directory base (overrides default)')
    parser.add_argument('--n-permutations', type=int, default=10000,
                        help='Locus permutations for the per-sample permutation null')
    parser.add_argument('--perm-seed', type=int, default=1234)
    args = parser.parse_args()

    # 如果指定了 --out-dir，覆盖默认输出目录
    if args.out_dir:
        OUT_BASE = args.out_dir

    run_analysis(args.name, args.features, args.op, args.exclude,
                 n_perm=args.n_permutations, perm_seed=args.perm_seed)
