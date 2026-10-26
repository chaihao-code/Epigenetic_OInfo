#!/usr/bin/env python3
"""08_region_filtering.py — Step 9: 免疫特异性区域排除敏感性分析.

三臂统一使用 hoi_recompute 重算估计器 (与原始分析相同的五分箱状态矩阵,
不重新离散化, 阈值沿用全基因组五分位):

  arm A: all original bins (重算)
  arm B: excluding immune-specific bins (|d_immune - other| >= 0.4, cov>=10, rep_sd<=0.08)
  arm C: excluding strongly immune-specific bins (|d| >= 0.5)
  arm M: matched downsampling control — 对每个区域从原始集合随机抽取与
         过滤后相同数量的 bins (构成按区域天然匹配), 100 次重复,
         报告年龄斜率的零分布.

注意: 该分析为 immune-marker-exclusion sensitivity analysis,
不等于重建纯甲状腺细胞的 O-information。

输出:
  results/tables/oinfo_region_filtering.tsv    (每区域×臂的 age 回归)
  results/tables/oinfo_filtering_values.tsv    (每区域×样本×臂的 O-info 值)
  results/figures/fig3_filtering_sensitivity.{pdf,png}
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
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from hoi_recompute import (autosomal_bin_ids, compute_region_oinfo,
                           load_sample_states, region_definitions)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(Path(__file__).resolve().parents[1] / "logs/09_region_filtering.log"),
        logging.StreamHandler(),
    ],
)
log = logging.getLogger("region_filtering")

ROOT = Path(__file__).resolve().parents[1]
TBL = ROOT / "results/tables"
FIG = ROOT / "results/figures"
N_REPS = 100
FOCUS_FAMILY = "proximal_ctcf_bound_enhancer_purified_bins"


def age_slope(hoi_by_sample):
    """hoi_by_sample: Series indexed by sample_id -> slope, p, n"""
    meta = pd.read_csv(ROOT / "data/sample_metadata/analysis_samples.tsv", sep="\t")
    d = hoi_by_sample.reset_index().merge(meta, on="sample_id")
    if len(d) < 4:
        return np.nan, np.nan
    f = sm.OLS(d.hoi.values.astype(float),
               sm.add_constant(d.age.values.astype(float))).fit()
    return float(f.params[1]), float(f.pvalues[1])


def main():
    cfg = _release_config(yaml.safe_load(open(ROOT / "config.yaml")))
    samples = cfg["samples"]["included"]
    bins_df, bin_ids = autosomal_bin_ids(cfg)
    log.info("autosomal bins: %d", len(bin_ids))
    bid_to_idx = {b: i for i, b in enumerate(bin_ids)}

    # --- 免疫特异 bin mask ---
    imm = pd.read_csv(ROOT / "data/processed/immune_specific_regions.bed", sep="\t")
    imm_idx = [bid_to_idx[f"{c}:{s}-{e}"] for c, s, e
               in zip(imm.chrom, imm.start, imm.end) if f"{c}:{s}-{e}" in bid_to_idx]
    immune_mask = np.zeros(len(bin_ids), dtype=bool)
    immune_mask[np.array(imm_idx)] = True
    strong_mask = np.zeros(len(bin_ids), dtype=bool)
    strong_idx = [bid_to_idx[f"{c}:{s}-{e}"] for c, s, e, d in
                  zip(imm.chrom, imm.start, imm.end, imm.d_immune_minus_other)
                  if f"{c}:{s}-{e}" in bid_to_idx and abs(d) >= 0.5]
    strong_mask[np.array(strong_idx)] = True
    log.info("immune-specific bins: %d (all), %d (strong)", immune_mask.sum(), strong_mask.sum())

    regions = region_definitions(cfg)
    log.info("regions: %d", len(regions))

    # --- 加载样本状态矩阵 ---
    Xs = {}
    for s in samples:
        Xs[s] = load_sample_states(cfg, s, n_expected=len(bin_ids))
        log.info("loaded %s", s)

    rng = np.random.default_rng(20261009)
    value_rows = []
    filt_rows = []
    for rname, mask in sorted(regions.items()):
        n_orig = int(mask.sum())
        arms = {
            "A_all": mask,
            "B_excl_immune": mask & ~immune_mask,
            "C_excl_strong": mask & ~strong_mask,
        }
        for arm, amask in arms.items():
            filt_rows.append(dict(region=rname, arm=arm, n_bins=int(amask.sum()),
                                  frac_of_original=round(amask.sum() / n_orig, 4)))
        # 每样本: 3 臂 + matched downsampling (B/C 各 N_REPS)
        for s in samples:
            for arm in ("A_all", "B_excl_immune", "C_excl_strong"):
                rec = compute_region_oinfo(Xs[s], arms[arm])
                value_rows.append(dict(region=rname, sample_id=s, arm=arm,
                                       rep=-1, hoi=rec["hoi"], n=rec["n"]))
            for base_arm, base_mask in (("B", arms["B_excl_immune"]), ("C", arms["C_excl_strong"])):
                idx_all = np.where(mask)[0]
                n_keep = int(base_mask.sum())
                if n_keep < 100 or n_keep >= n_orig:
                    continue
                for rep in range(N_REPS):
                    # matched control: 从原始区域集合随机抽同样数量 bins
                    # (不按免疫状态排除, 构成按区域天然匹配)
                    sub = rng.choice(idx_all, size=n_keep, replace=False)
                    rec = oinfo_from_states_sub(Xs[s], sub)
                    value_rows.append(dict(region=rname, sample_id=s,
                                           arm=f"M_down_{base_arm}", rep=rep,
                                           hoi=rec["hoi"], n=rec["n"]))
        log.info("done %s (n_orig=%d)", rname, n_orig)

    vals = pd.DataFrame(value_rows)
    vals.to_csv(TBL / "oinfo_filtering_values.tsv", sep="\t", index=False)
    pd.DataFrame(filt_rows).to_csv(TBL / "region_filtering_stats.tsv", sep="\t", index=False)

    # --- age 回归 ---
    rows = []
    for (rname, arm), g in vals.groupby(["region", "arm"]):
        if arm.startswith("M_down"):
            slopes = []
            for rep, gr in g.groupby("rep"):
                sl, pv = age_slope(gr.set_index("sample_id")["hoi"])
                if np.isfinite(sl):
                    slopes.append(sl)
            slopes = np.array(slopes)
            rows.append(dict(region=rname, arm=arm, n=len(slopes),
                             age_slope=np.mean(slopes),
                             slope_sd=np.std(slopes),
                             slope_q025=np.percentile(slopes, 2.5),
                             slope_q975=np.percentile(slopes, 97.5),
                             p_value=np.nan))
        else:
            sl, pv = age_slope(g.set_index("sample_id")["hoi"])
            rows.append(dict(region=rname, arm=arm, n=1, age_slope=sl,
                             slope_sd=np.nan, slope_q025=np.nan,
                             slope_q975=np.nan, p_value=pv))
    reg = pd.DataFrame(rows)
    reg.to_csv(TBL / "oinfo_region_filtering.tsv", sep="\t", index=False)

    # --- 汇总: 过滤臂斜率 vs matched 零分布分位数 ---
    summ = []
    for rname in sorted(reg.region.unique()):
        for base, armname in (("B", "B_excl_immune"), ("C", "C_excl_strong")):
            obs = reg[(reg.region == rname) & (reg.arm == armname)]
            null = reg[(reg.region == rname) & (reg.arm == f"M_down_{base}")]
            if len(obs) == 0 or len(null) == 0:
                continue
            obs_sl = obs.age_slope.iloc[0]
            # 从原始值重算每个 rep 的斜率, 求观测斜率在零分布中的分位数
            rep_slopes = []
            null_vals = vals[(vals.region == rname) & (vals.arm == f"M_down_{base}")]
            for rep, gr in null_vals.groupby("rep"):
                sl, _ = age_slope(gr.set_index("sample_id")["hoi"])
                if np.isfinite(sl):
                    rep_slopes.append(sl)
            rep_slopes = np.array(rep_slopes)
            pct = float((rep_slopes < obs_sl).mean()) if len(rep_slopes) else np.nan
            summ.append(dict(region=rname, arm=armname,
                             slope_excluded=obs_sl,
                             slope_matched_mean=null.age_slope.iloc[0],
                             slope_matched_q025=null.slope_q025.iloc[0],
                             slope_matched_q975=null.slope_q975.iloc[0],
                             p_unadj=obs.p_value.iloc[0],
                             excluded_below_null_pct=pct,
                             n_reps=len(rep_slopes)))
    pd.DataFrame(summ).to_csv(TBL / "oinfo_filtering_summary.tsv", sep="\t", index=False)

    log.info("focus family result:")
    foc = reg[reg.region.str.contains(FOCUS_FAMILY)]
    print(foc.pivot_table(index="region", columns="arm", values="age_slope").to_string())

    # ---------------- Figure 3 ----------------
    meta = pd.read_csv(ROOT / "data/sample_metadata/analysis_samples.tsv", sep="\t")
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))

    # 3A: focus region O-info vs age, 3 arms
    ax = axes[0]
    focus_region = f"{FOCUS_FAMILY}_AND_None"
    cols = {"A_all": "gray", "B_excl_immune": "crimson", "C_excl_strong": "darkorange"}
    for arm, c in cols.items():
        d = vals[(vals.region == focus_region) & (vals.arm == arm)].merge(meta, on="sample_id")
        f = sm.OLS(d.hoi.values.astype(float),
                   sm.add_constant(d.age.values.astype(float))).fit()
        xs = np.linspace(d.age.min() - 2, d.age.max() + 2, 20)
        ax.scatter(d.age, d.hoi, color=c, s=40, edgecolor="k", label=arm, zorder=3)
        ax.plot(xs, f.params[0] + f.params[1] * xs, color=c, lw=1,
                label=f"slope={f.params[1]:.4f}")
    ax.set_xlabel("age (y)")
    ax.set_ylabel("O-information (recomputed)")
    ax.set_title(f"3A: {focus_region.replace('_bins_AND_','/')}", fontsize=8)
    ax.legend(fontsize=6)

    # 3B: slopes across all regions (A vs B vs matched mean)
    ax = axes[1]
    pv = reg.pivot_table(index="region", columns="arm", values="age_slope")
    key = pv[pv.index.str.contains("proximal_ctcf|proximal_not_ctcf")].sort_values("A_all")
    x = np.arange(len(key))
    ax.errorbar(x, key["A_all"], fmt="o", ms=4, color="gray", label="A: all bins")
    ax.errorbar(x, key["B_excl_immune"], fmt="s", ms=4, color="crimson",
                label="B: excl immune")
    if "M_down_B" in key.columns and key["M_down_B"].notna().any():
        ax.plot(x, key["M_down_B"], "^", ms=4, color="pink",
                label="M: matched random")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels([r.replace("_bins_AND_", "\n") for r in key.index],
                       rotation=90, fontsize=5)
    ax.set_ylabel("age slope of O-info")
    ax.set_title("3B: enhancer regions, all vs excl-immune", fontsize=9)
    ax.legend(fontsize=7)

    # 3C: matched downsampling null distribution for focus region
    ax = axes[2]
    fd = reg[(reg.region == focus_region) & (reg.arm.str.startswith("M_down"))]
    nul = vals[(vals.region == focus_region) & (vals.arm.str.startswith("M_down"))]
    rep_sl = []
    for (arm, rep), gr in nul.groupby(["arm", "rep"]):
        sl, _ = age_slope(gr.set_index("sample_id")["hoi"])
        if np.isfinite(sl):
            rep_sl.append(dict(arm=arm, slope=sl))
    rep_sl = pd.DataFrame(rep_sl)
    for arm, c, lab in (("M_down_B", "pink", "matched to B (n)"),
                        ("M_down_C", "moccasin", "matched to C (n)")):
        if len(rep_sl):
            d = rep_sl[rep_sl.arm == arm]
            ax.hist(d.slope, bins=25, alpha=0.7, color=c, label=f"{lab} mean={d.slope.mean():.4f}")
    for arm, c in (("B_excl_immune", "crimson"), ("C_excl_strong", "darkorange"),
                   ("A_all", "gray")):
        o = reg[(reg.region == focus_region) & (reg.arm == arm)]
        if len(o):
            ax.axvline(o.age_slope.iloc[0], color=c, ls="--", lw=1.5, label=f"obs {arm}")
    ax.set_xlabel("age slope (100 matched random subsets)")
    ax.set_title("3C: matched downsampling null", fontsize=9)
    ax.legend(fontsize=6)

    fig.tight_layout()
    fig.savefig(FIG / "fig3_filtering_sensitivity.pdf")
    fig.savefig(FIG / "fig3_filtering_sensitivity.png", dpi=200)
    plt.close(fig)
    log.info("saved fig3")
    return 0


def oinfo_from_states_sub(X, row_idx):
    from hoi_recompute import oinfo_from_states, subsample_indices
    sel = subsample_indices(len(row_idx))
    return oinfo_from_states(X[row_idx[sel]])


if __name__ == "__main__":
    sys.exit(main())
