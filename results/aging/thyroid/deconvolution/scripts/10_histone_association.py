#!/usr/bin/env python3
"""10_histone_association.py — Step 10b: 组蛋白修饰变化与免疫比例的关联.

使用原始连续信号 (input 归一化 bigWig)。为避免跨样本 ChIP-seq
绝对定量不可靠的问题, 同时报告三种度量:
  raw: 区域内 bins 的中位数 (原值)
  locus_relative: log2(区域中位 / 该样本全基因组中位) — 样本内自校准
  rank: 6 个样本间的排名 (1-6, 完全稳健于单调变换)

对 40 个区域 × {DNA_methylation, H3K4me1, H3K4me3, H3K27ac, H3K27me3, H3K9me3}
分别计算与 age / estimated immune fraction 的 Spearman 相关及 OLS 斜率
(locus_relative 度量)。 n=6, 报告 effect size, 不做过度显著性解释。

输出:
  results/tables/histone_immune_association.tsv
  results/figures/fig5_histone_immune.{pdf,png}
"""
import logging
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig
import statsmodels.api as sm
import yaml
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from hoi_recompute import autosomal_bin_ids, region_definitions

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(Path(__file__).resolve().parents[1] / "logs/11_histone.log"),
        logging.StreamHandler(),
    ],
)
log = logging.getLogger("histone")

ROOT = Path(__file__).resolve().parents[1]
TBL = ROOT / "results/tables"
FIG = ROOT / "results/figures"
SAMPLES = ["CEMT0098", "CEMT0130", "CEMT0132", "CEMT0042", "CEMT0044", "CEMT0040"]
AGES = [30, 37, 37, 46, 55, 67]
HISTONE_MARKS = ["H3K4me1", "H3K4me3", "H3K27ac", "H3K27me3", "H3K9me3", "DNA_methylation"]
FOCUS = "proximal_ctcf_bound_enhancer_purified_bins_AND_None"


def find_bw(bw_dir, sample, mark):
    if mark == "DNA_methylation":
        pat = f"*CEEHRC.{sample}.WGB-Seq.methylation_profile.bigWig"
    else:
        pat = f"*CEEHRC.{sample}.{mark}.signal_unstranded.bigWig"
    hits = sorted(Path(bw_dir).glob(pat))
    hits = [h for h in hits if not h.name.endswith(".tbi")]
    if not hits:
        raise FileNotFoundError(pat)
    # 同一信号可能有多个下载版本, 取最大(最完整)的文件
    return max(hits, key=lambda h: h.stat().st_size)


def load_bin_means(bw, bins_df):
    """按 bins_df 行序返回全常染色体 1kb bin 均值数组。

    values(numpy=True) 一次 C 传递取整条染色体 (float32, NaN=无覆盖),
    再 reshape 成 (n_kb, 1000) 做向量化 1kb 均值 (忽略染色体尾端余数)。
    bins_df 必须只含常染色体, 1kb 对齐 (断言)。
    """
    out = np.full(len(bins_df), np.nan, dtype=float)
    pos = 0
    for chrom, grp in bins_df.groupby("chrom", sort=False):
        s = grp.start.values
        e = grp.end.values
        assert (np.diff(s) == 1000).all() and s[0] % 1000 == 0, f"{chrom}: bins 非 1kb 对齐"
        clen = bw.chroms(chrom)
        nb = len(grp)
        assert s[0] == 0 and (np.diff(s) == 1000).all(), f"{chrom}: bins 非 1kb 对齐"
        nb_full = min(nb, clen // 1000)   # 完整 1kb 槽位数
        if nb_full:
            v = np.asarray(bw.values(chrom, 0, nb_full * 1000, numpy=True),
                           dtype=np.float32).reshape(nb_full, 1000)
            with np.errstate(invalid="ignore"):
                out[pos:pos + nb_full] = np.nanmean(v, axis=1)
        if nb > nb_full:  # 尾端不足 1kb 的 bin
            v = np.asarray(bw.values(chrom, nb_full * 1000, clen, numpy=True),
                           dtype=np.float32)
            with np.errstate(invalid="ignore"):
                out[pos + nb_full] = np.nanmean(v) if v.size else np.nan
        pos += nb
    return out


def main():
    cfg = yaml.safe_load(open(ROOT / "config.yaml"))
    bw_dir = cfg["data"]["bigwig_dir"]
    bins_df, _ = autosomal_bin_ids(cfg)
    regions = region_definitions(cfg)
    imm = pd.read_csv(TBL / "thyroid_cell_proportions.tsv", sep="\t")
    imm_frac = (imm[(imm.panel == "A") & (imm.model == "A")
                    & imm.cell_type.isin(["T_CD4", "T_CD8", "B_cell", "Monocyte", "NK", "Granulocyte"])]
                .groupby("sample_id")["fraction"].sum())

    rows = []
    for mark in HISTONE_MARKS:
        # 每样本一次性载入全常染色体 1kb 均值 (~23MB/样本)
        sample_means = {}
        cache_dir = ROOT / "data/processed/bin_means_cache"
        cache_dir.mkdir(exist_ok=True)
        for s in SAMPLES:
            cf = cache_dir / f"{mark}__{s}.npy"
            if cf.exists():
                sample_means[s] = np.load(cf)
                continue
            try:
                bw = pyBigWig.open(str(find_bw(bw_dir, s, mark)))
                sample_means[s] = load_bin_means(bw, bins_df)
                bw.close()
            except Exception as e:
                # 文件损坏/缺失: 该样本该 mark 记 NaN (如 CEMT0044 H3K27me3 bigWig 损坏)
                log.warning("mark %s sample %s unavailable: %s", mark, s, e)
                sample_means[s] = np.full(len(bins_df), np.nan)
            if np.isfinite(sample_means[s]).any():
                np.save(cf, sample_means[s])
            else:
                log.error("mark %s sample %s: all-NaN, not cached", mark, s)
        # 全基因组中位 (样本自校准基准)
        genome_med = {s: np.nanmedian(sample_means[s]) for s in SAMPLES}
        log.info("%s genome medians: %s", mark,
                 {k: round(v, 3) for k, v in genome_med.items()})

        for rname, mask in sorted(regions.items()):
            idx = np.where(mask)[0]
            meds = {}
            for s in SAMPLES:
                v = sample_means[s][idx]
                meds[s] = np.nanmedian(v) if np.isfinite(v).any() else np.nan
            raw = np.array([meds[s] for s in SAMPLES], dtype=float)
            lr = np.log2(raw / np.array([genome_med[s] for s in SAMPLES]))
            rk = pd.Series(raw).rank().values  # NaN 保持 NaN (scipy rankdata 会传播 NaN)
            for metric, y in (("raw", raw), ("locus_relative", lr), ("rank", rk)):
                # 成对完整观测 (个别样本×mark 缺失, 如 CEMT0044 H3K27me3)
                ys = pd.Series(y, index=SAMPLES)
                age_s = pd.Series(AGES, index=SAMPLES, dtype=float)
                imm_s = imm_frac.loc[SAMPLES].astype(float)
                ra, pa = (ys.corr(age_s, method="spearman"),) * 2
                from scipy import stats as _st
                msk = ys.notna() & age_s.notna()
                ra, pa = _st.spearmanr(age_s[msk], ys[msk]) if msk.sum() >= 4 else (np.nan, np.nan)
                msk_i = ys.notna() & imm_s.notna()
                ri, pi = _st.spearmanr(imm_s[msk_i], ys[msk_i]) if msk_i.sum() >= 4 else (np.nan, np.nan)
                if msk.sum() >= 4:
                    f = sm.OLS(ys[msk].values.astype(float),
                               sm.add_constant(age_s[msk].values.astype(float))).fit()
                    slope, p_age = float(f.params[1]), float(f.pvalues[1])
                else:
                    slope, p_age = np.nan, np.nan
                rows.append(dict(mark=mark, region=rname, metric=metric,
                                 spearman_r_age=ra, spearman_p_age=pa,
                                 spearman_r_immune=ri, spearman_p_immune=pi,
                                 ols_slope_per_year=slope,
                                 ols_p_age=p_age,
                                 n_complete=int(msk.sum()),
                                 values_by_sample=";".join(f"{s}={v:.4g}" for s, v in zip(SAMPLES, y))))
        log.info("done %s", mark)

    out = pd.DataFrame(rows)
    out.to_csv(TBL / "histone_immune_association.tsv", sep="\t", index=False)

    # 重点: 变化方向是否符合审稿人描述 (H3K4me1/3, H3K27ac ↓; H3K27me3, methylation ↑)
    main = out[out.metric == "locus_relative"]
    foc = main[main.region == FOCUS].set_index("mark")
    log.info("focus region locus-relative slopes:\n%s",
             foc[["ols_slope_per_year", "ols_p_age", "spearman_r_age", "spearman_r_immune"]])

    # ---------------- Figure 5 ----------------
    fig, axes = plt.subplots(2, 3, figsize=(14, 7))
    for k, mark in enumerate(HISTONE_MARKS):
        ax = axes[k // 3, k % 3]
        d = main[main.mark == mark]
        key = d[d.region.str.contains("proximal_ctcf|proximal_not_ctcf")].sort_values("ols_slope_per_year")
        x = np.arange(len(key))
        ax.bar(x, key.ols_slope_per_year, color=np.where(key.ols_slope_per_year < 0, "steelblue", "indianred"))
        ax.axhline(0, color="k", lw=0.6)
        ax.set_xticks(x)
        ax.set_xticklabels([r.replace("_bins_AND_", "\n") for r in key.region], rotation=90, fontsize=5)
        ax.set_ylabel("slope (locus-relative / y)")
        ax.set_title(mark, fontsize=9)
    fig.suptitle(f"locus-relative signal change with age — enhancer regions (focus: {FOCUS.replace('_bins_AND_','/')})", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(FIG / "fig5_histone_immune.pdf")
    fig.savefig(FIG / "fig5_histone_immune.png", dpi=200)
    plt.close(fig)

    # 散点: 重点区每个 mark 的 locus-relative vs immune fraction
    fig, axes = plt.subplots(2, 3, figsize=(14, 7))
    for k, mark in enumerate(HISTONE_MARKS):
        ax = axes[k // 3, k % 3]
        d = main[(main.mark == mark) & (main.region == FOCUS)].iloc[0]
        y = np.array([float(v.split("=")[1]) for v in d.values_by_sample.split(";")])
        xi = imm_frac.loc[SAMPLES].values
        ax.scatter(xi, y, s=50, color="seagreen", edgecolor="k")
        r, p = stats.spearmanr(xi, y)
        ax.set_xlabel("immune fraction")
        ax.set_ylabel("locus-relative signal")
        ax.set_title(f"{mark}: r_immune={r:.2f} p={p:.3g}", fontsize=9)
    fig.suptitle(f"focus region signal vs immune fraction ({FOCUS.replace('_bins_AND_','/')})", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(FIG / "fig5b_histone_vs_immune.pdf")
    fig.savefig(FIG / "fig5b_histone_vs_immune.png", dpi=200)
    plt.close(fig)
    log.info("saved fig5")
    return 0


if __name__ == "__main__":
    sys.exit(main())
