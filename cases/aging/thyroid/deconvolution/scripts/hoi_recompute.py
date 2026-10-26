#!/usr/bin/env python3
"""hoi_recompute.py — 从离散化矩阵精确重算 7 变量 O-information.

复现原始管线 (code/hoi_analysis/compute_hoi_oinfo.py):
  - 输入: data/aging_thyroid_input_normed_binned_discreted/*.discrete.bed (0-4 状态,
    行序 = bins_1kb 常染色体序, 与原 signal_loader 输出一致)
  - 区域定义: block_processor.cartesian_product_blocks 逻辑
    group1 = [None, promoter, proximal_ctcf_bound_enhancer_purified_bins,
              proximal_not_ctcf_bound_enhancer_purified_bins,
              distal_ctcf_bound_enhancer_purified_bins,
              distal_not_ctcf_bound_enhancer_purified_bins,
              exon_purified_bins, intron_purified_bins]
    group2 = [CpG_island, LINE, LTR, SINE, None]
    None = all_bins - union(group)
  - 熵: 离散频率计数, H = -sum(p log2 p) (与原 binning 方法一致)
  - O-info: hoi = tc - dtc;  tc = sum(H_i) - H_all;  dtc = sum(H(all\j)) - (n-1)*H_all
    (shared information_from_codes API)
  - n_max=100000 子采样: pandas df.sample(random_state=42), 行序 = 常染色体 bin 序

本模块同时提供 validate 模式: 对全部 40 区域 × 样本与磁盘上原始
results_oinfo_n_bins_5 结果逐值比对, 作为公式/区域/状态正确性验证。
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
import logging
import sys
sys.path.insert(0, str(_RELEASE_ROOT / "common/02_oinformation"))
from compute_hoi_oinfo import information_from_codes, entropy_of_states
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import pyBigWig  # noqa: F401  (环境检查)
import yaml

log = logging.getLogger("hoi_recompute")

N_FEATURES = 7
N_MAX = 100000
RANDOM_SEED = 42

GROUP1_MASKS = [
    ("promoter", "promoter_bins.tsv"),
    ("proximal_ctcf_bound_enhancer_purified_bins", "proximal_ctcf_bound_enhancer_purified_bins.tsv"),
    ("proximal_not_ctcf_bound_enhancer_purified_bins", "proximal_not_ctcf_bound_enhancer_purified_bins.tsv"),
    ("distal_ctcf_bound_enhancer_purified_bins", "distal_ctcf_bound_enhancer_purified_bins.tsv"),
    ("distal_not_ctcf_bound_enhancer_purified_bins", "distal_not_ctcf_bound_enhancer_purified_bins.tsv"),
    ("exon_purified_bins", "exon_purified_bins.tsv"),
    ("intron_purified_bins", "intron_purified_bins.tsv"),
]
GROUP2_MASKS = [
    ("CpG_island", "CpG_island_bins.tsv"),
    ("LINE", "LINE_bins.tsv"),
    ("LTR", "LTR_bins.tsv"),
    ("SINE", "SINE_bins.tsv"),
]


def load_config():
    return _release_config(yaml.safe_load(open(Path(__file__).resolve().parents[1] / "config.yaml")))


def autosomal_bin_ids(cfg):
    """bins_1kb 中常染色体 bin 的 bin_id 列表 (chrN:start-end)。"""
    bins = pd.read_csv(cfg["data"]["bins_1kb"], sep="\t", header=None,
                       names=["chrom", "start", "end"])
    bins = bins[bins["chrom"].str.match(r"^chr\d+$")].reset_index(drop=True)
    return bins, (bins["chrom"] + ":" + bins["start"].astype(str) + "-" + bins["end"].astype(str)).tolist()


def load_mask(cfg, fname, bin_ids):
    """读 mask tsv -> 与 bin_ids 对齐的 bool 数组。缺失 bin 视为 0。"""
    df = pd.read_csv(Path(cfg["data"]["genome_features_dir"]) / fname, sep="\t")
    col = df.columns[1]
    s = df.set_index("bin_id")[col]
    return s.reindex(bin_ids).fillna(0).astype(bool).values


def region_definitions(cfg):
    """按 cartesian_product_blocks 逻辑构造全部 40 个区域的 bin 索引数组。

    返回 dict: region_name -> np.ndarray (bool, len = n_autosomal_bins)
    """
    _, bin_ids = autosomal_bin_ids(cfg)
    n = len(bin_ids)
    all_bins = np.ones(n, dtype=bool)

    def build(group_masks):
        group = {}
        union = np.zeros(n, dtype=bool)
        for name, fname in group_masks:
            m = load_mask(cfg, fname, bin_ids)
            group[name] = m
            union |= m
        group["None"] = all_bins & ~union
        return group

    g1 = build(GROUP1_MASKS)
    g2 = build(GROUP2_MASKS)
    regions = {}
    for n1, m1 in g1.items():
        for n2, m2 in g2.items():
            inter = m1 & m2
            if inter.any():
                regions[f"{n1}_AND_{n2}"] = inter
    return regions


def find_discrete_files(cfg, sample):
    """定位样本的 7 个离散化 bed。返回 mark -> path (特征序无关, O-info 对称)。"""
    ddir = Path(cfg["data"]["discretized_dir"])
    files = {}
    for mark in cfg["marks"]:
        if mark == "DNA_methylation":
            pat = f"*CEEHRC.{sample}.WGB-Seq.methylation_profile.*discrete.bed*"
        else:
            pat = f"*CEEHRC.{sample}.{mark}.signal_unstranded.*discrete.bed*"
        hits = sorted(ddir.glob(pat))
        hits = [h for h in hits if not h.name.endswith(".tbi")]
        if not hits:
            raise FileNotFoundError(f"{sample} {mark}: {pat}")
        files[mark] = hits[0]
    return files


def load_sample_states(cfg, sample, n_expected=None):
    """加载样本 7 轨状态矩阵, 行序 = bins_1kb 常染色体序。"""
    files = find_discrete_files(cfg, sample)
    cols = []
    for mark in cfg["marks"]:
        df = pd.read_csv(files[mark], sep="\t", header=None,
                         names=["chrom", "start", "end", "state"],
                         usecols=[3], dtype={"state": np.int8})
        cols.append(df["state"].values)
    X = np.column_stack(cols)
    if n_expected is not None and len(X) != n_expected:
        raise ValueError(f"{sample}: 行数 {len(X)} != 期望 {n_expected}")
    return X


def subsample_indices(n_rows):
    """复现 df.sample(n=N_MAX, random_state=42) 的位置索引 (pandas 行为)。"""
    if n_rows <= N_MAX:
        return np.arange(n_rows)
    df = pd.DataFrame(index=np.arange(n_rows))
    return df.sample(n=N_MAX, random_state=RANDOM_SEED).index.values


def joint_entropy(states):
    return entropy_of_states(states, n_states=5)


def oinfo_from_states(X):
    return information_from_codes(X, n_states=5)


def compute_region_oinfo(X, region_mask):
    idx = np.where(region_mask)[0]
    sel = subsample_indices(len(idx))
    Xb = X[idx[sel]]
    return oinfo_from_states(Xb)


def validate(cfg, samples, out_tsv):
    """与原始 results_oinfo_n_bins_5 逐值比对。"""
    _, bin_ids = autosomal_bin_ids(cfg)
    regions = region_definitions(cfg)
    odir = Path(cfg["data"]["oinfo_results_dir"])
    rows = []
    for sample in samples:
        X = load_sample_states(cfg, sample, n_expected=len(bin_ids))
        for rname, mask in sorted(regions.items()):
            f = odir / f"CEEHRC.{sample}/{rname}/CEEHRC.{sample}_{rname}_oinfo_overall.csv"
            if not f.exists():
                continue
            orig = pd.read_csv(f)
            orig = orig[(orig["order"] == 7) & (orig["n_missing"] == 0)].iloc[0]
            rec = compute_region_oinfo(X, mask)
            rows.append(dict(
                sample=sample, region=rname,
                n_recomputed=rec["n"],
                hoi_orig=float(orig["hoi"]), hoi_recomp=rec["hoi"],
                tc_orig=float(orig["tc"]), tc_recomp=rec["tc"],
                dtc_orig=float(orig["dtc"]), dtc_recomp=rec["dtc"],
                abs_diff_hoi=abs(float(orig["hoi"]) - rec["hoi"]),
            ))
        log.info("validated %s (%d regions)", sample, len(rows))
    df = pd.DataFrame(rows)
    df.to_csv(out_tsv, sep="\t", index=False)
    return df


if __name__ == "__main__":
    import sys
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    cfg = load_config()
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else \
        Path(__file__).resolve().parents[1] / "results/qc/hoi_recompute_validation.tsv"
    samples = cfg["samples"]["included"]
    df = validate(cfg, samples, out)
    print(df.groupby("sample")["abs_diff_hoi"].agg(["max", "mean"]))
    print(f"\nwrote {out}")
