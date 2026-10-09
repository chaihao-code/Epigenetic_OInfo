#!/usr/bin/env python3
"""deconv_utils.py — 去卷积共享工具: R 矩阵加载, 约束 NNLS, block bootstrap.

约束: p >= 0, sum(p) = 1 (若模型细胞类型不完整, 估计为相对比例, 需另行报告)。
"""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import nnls

import yaml

ROOT = Path(__file__).resolve().parents[1]

# 命名约定 (2026-10-09 起): 主分析为 Panel A / Model A
# Panel A = 独立筛选 marker (主); Panel B = UXM U250 文献 marker (敏感性)
# Model A = 9 类细粒度 (主);   Model B = 4 类粗粒度 (敏感性)
MODEL_A = ["Thyroid_epithelial", "T_CD4", "T_CD8", "B_cell", "Monocyte", "NK",
           "Granulocyte", "Endothelial", "Fibroblast"]
MODEL_B = ["Thyroid", "Lymphoid", "Myeloid", "Stromal"]
MERGE_B = {  # Model A 细类 -> Model B 粗类
    "Thyroid_epithelial": "Thyroid",
    "T_CD4": "Lymphoid", "T_CD8": "Lymphoid", "B_cell": "Lymphoid",
    "Monocyte": "Myeloid", "NK": "Myeloid", "Granulocyte": "Myeloid",
    "Endothelial": "Stromal", "Fibroblast": "Stromal",
}


def load_R(panel="A", model="A"):
    """返回 (R: (n_marker, k) ndarray, marker_ids, cell_types, marker_chroms)."""
    df = pd.read_csv(ROOT / f"data/processed/reference_matrix_{panel}.tsv", sep="\t")
    wide = df.pivot_table(index="marker_id", columns="cell_type", values="methylation")
    if model == "B":
        wide = wide.groupby(MERGE_B, axis=1).mean()
        cell_types = MODEL_B
    else:
        cell_types = [c for c in MODEL_A if c in wide.columns]
    wide = wide[cell_types].dropna()
    # 每条 marker 的染色体 (block bootstrap 用)
    bed = ROOT / f"data/processed/selected_markers_{panel}.bed"
    chrom = [None] * len(wide)
    if bed.exists():
        if panel == "A":
            b = pd.read_csv(bed, sep="\t")
            lut = dict(zip(b["marker_id"], b["chrom"]))
        else:
            b = pd.read_csv(bed, sep="\t", header=None, names=["chrom", "start", "end", "name"])
            lut = dict(zip(b["name"], b["chrom"]))
        chrom = [lut.get(m) for m in wide.index]
    return wide.values, list(wide.index), cell_types, chrom


def fit_nnls(y, R, ridge=0.0):
    """约束 NNLS: min ||y - R p||^2 + ridge||p||^2, s.t. p>=0, sum p=1.

    用增广矩阵法: 在 R 追加一行全 1 (权重 lam), y 追加 lam, 再解 NNLS,
    lam 逐步增大直至和为 1 (迭代投影, 通常 1-3 次收敛)。
    """
    y = np.asarray(y, dtype=np.float64)
    R = np.asarray(R, dtype=np.float64)
    ok = ~np.isnan(y)
    if ok.sum() < R.shape[1]:
        return None
    Rk, yk = R[ok], y[ok]
    if ridge > 0:
        Rk = np.vstack([Rk, np.sqrt(ridge) * np.eye(R.shape[1])])
        yk = np.concatenate([yk, np.zeros(R.shape[1])])
    lam = 10.0
    p, _ = nnls(np.vstack([Rk, lam * np.ones((1, R.shape[1]))]),
                np.concatenate([yk, [lam]]))
    for _ in range(50):
        s = p.sum()
        if abs(s - 1) < 1e-6:
            break
        lam *= max(2.0, 1.0 / max(s, 1e-9))
        p, _ = nnls(np.vstack([Rk, lam * np.ones((1, R.shape[1]))]),
                    np.concatenate([yk, [lam]]))
    return p / max(p.sum(), 1e-12)


def fit_bootstrap(y, R, marker_chroms, n_boot=200, seed=42):
    """染色体 block bootstrap: 重抽样整染色体上的 marker, 重拟合。"""
    rng = np.random.default_rng(seed)
    chroms = np.array([c if c else "?" for c in marker_chroms])
    unique = np.unique(chroms)
    boots = []
    for _ in range(n_boot):
        pick = rng.choice(unique, size=len(unique), replace=True)
        idx = np.concatenate([np.where(chroms == c)[0] for c in pick])
        p = fit_nnls(y[idx], R[idx])
        if p is not None:
            boots.append(p)
    if not boots:
        return None, None
    boots = np.array(boots)
    return boots.mean(axis=0), np.percentile(boots, [2.5, 97.5], axis=1)


def aggregate_B(p_a):
    """Model A (9 类细) 比例 -> Model B (4 类粗) 比例。"""
    out = np.zeros(4)
    for i, ct in enumerate(MODEL_A):
        out[MODEL_B.index(MERGE_B[ct])] += p_a[i]
    return out
