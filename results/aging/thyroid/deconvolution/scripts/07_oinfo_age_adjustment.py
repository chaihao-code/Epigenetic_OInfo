#!/usr/bin/env python3
"""07_oinfo_age_adjustment.py — Step 8: O-information 年龄效应的免疫比例调整.

对 40 个 genomic context (存储的原始 O-info 值) 分别拟合:
  Model 0: hoi ~ age
  Model 1: hoi ~ age + immune_fraction
报告 beta_age(调整前/后), beta_immune, 95% CI, p, adj R^2, attenuation,
VIF, leave-one-out 敏感性, panel/model 组合间稳定性, BH-FDR。
免疫比例取 Panel A / Model A 主分析, Panel B / Model B 为敏感性。
输出:
  results/tables/oinfo_age_adjustment.tsv
  results/figures/fig2_oinfo_age_adjustment.{pdf,png}
"""
import logging
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
import yaml

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(Path(__file__).resolve().parents[1] / "logs/08_oinfo_adjust.log"),
        logging.StreamHandler(),
    ],
)
log = logging.getLogger("oinfo_adjust")

ROOT = Path(__file__).resolve().parents[1]
TBL = ROOT / "results/tables"
FIG = ROOT / "results/figures"
# 免疫组分按 model 粒度定义 (Model A 细粒度 9 类, Model B 粗粒度 4 类)
IMM_BY_MODEL = {
    "A": ["T_CD4", "T_CD8", "B_cell", "Monocyte", "NK", "Granulocyte"],
    "B": ["Lymphoid", "Myeloid"],
}
FOCUS = "proximal_ctcf_bound_enhancer_purified_bins_AND_None"


def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        val = p[order[i]] * n / rank
        prev = min(prev, val)
        q[order[i]] = prev
    return np.minimum(q, 1.0)


def vif(x1, x2):
    r = np.corrcoef(x1, x2)[0, 1]
    return 1 / (1 - r ** 2) if abs(r) < 1 else np.inf


def fit_region(df, region, panel, model):
    d = df[df.region == region].copy()
    if len(d) < 4:
        return None
    x_age = d.age.values.astype(float)
    x_imm = d.immune_fraction.values.astype(float)
    y = d.hoi.values.astype(float)

    X0 = sm.add_constant(x_age)
    f0 = sm.OLS(y, X0).fit()
    X1 = sm.add_constant(np.column_stack([x_age, x_imm]))
    f1 = sm.OLS(y, X1).fit()

    ci = np.asarray(f1.conf_int())
    ci0 = np.asarray(f0.conf_int())[1]
    ci1_age = ci[1]
    ci1_imm = ci[2]

    b0, b1 = f0.params[1], f1.params[1]
    attenuation = (1 - b1 / b0) if abs(b0) > 1e-9 else np.nan

    # leave-one-out on adjusted model
    loo_slopes = []
    for s in d.sample_id:
        m = d.sample_id != s
        if m.sum() < 4:
            continue
        f = sm.OLS(y[m], sm.add_constant(np.column_stack([x_age[m], x_imm[m]]))).fit()
        loo_slopes.append(f.params[1])
    # bootstrap CI for adjusted age slope (donor resample)
    rng = np.random.default_rng(42)
    boots = []
    idx = np.arange(len(d))
    for _ in range(1000):
        bidx = rng.choice(idx, size=len(idx), replace=True)
        if len(np.unique(x_age[bidx])) < 3:
            continue
        try:
            f = sm.OLS(y[bidx], sm.add_constant(np.column_stack([x_age[bidx], x_imm[bidx]]))).fit()
            boots.append(f.params[1])
        except Exception:
            continue
    boot_ci = np.percentile(boots, [2.5, 97.5]) if len(boots) > 50 else [np.nan, np.nan]

    return dict(
        region=region, panel=panel, model=model, n=len(d),
        beta_age_unadj=b0, ci_lo_unadj=ci0[0], ci_hi_unadj=ci0[1],
        p_unadj=f0.pvalues[1], adj_r2_unadj=f0.rsquared_adj,
        beta_age_adj=b1, ci_lo_adj=ci1_age[0], ci_hi_adj=ci1_age[1],
        p_adj=f1.pvalues[1], adj_r2_adj=f1.rsquared_adj,
        beta_immune=f1.params[2], ci_lo_immune=ci1_imm[0], ci_hi_immune=ci1_imm[1],
        p_immune=f1.pvalues[2],
        attenuation=attenuation if np.isfinite(attenuation) else np.nan,
        vif_age_immune=vif(x_age, x_imm),
        r_age_immune=np.corrcoef(x_age, x_imm)[0, 1],
        loo_slope_min=float(np.min(loo_slopes)) if loo_slopes else np.nan,
        loo_slope_max=float(np.max(loo_slopes)) if loo_slopes else np.nan,
        boot_ci_lo=float(boot_ci[0]), boot_ci_hi=float(boot_ci[1]),
    )


def main():
    oinfo = pd.read_csv(ROOT / "data/processed/oinfo_long.tsv", sep="\t")
    prop = pd.read_csv(TBL / "thyroid_cell_proportions.tsv", sep="\t")
    meta = pd.read_csv(ROOT / "data/sample_metadata/analysis_samples.tsv", sep="\t")

    rows = []
    for panel, model in [("A", "A"), ("B", "B")]:
        imm = (prop[(prop.panel == panel) & (prop.model == model)
                    & prop.cell_type.isin(IMM_BY_MODEL[model])]
               .groupby("sample_id")["fraction"].sum().reset_index()
               .rename(columns={"fraction": "immune_fraction"}))
        df = oinfo.merge(meta, on="sample_id").merge(imm, on="sample_id")
        for region in sorted(df.region.unique()):
            r = fit_region(df, region, panel, model)
            if r:
                rows.append(r)
    out = pd.DataFrame(rows)
    out["fdr_unadj"] = np.nan
    out["fdr_adj"] = np.nan
    for (pn, md), gidx in out.groupby(["panel", "model"]).groups.items():
        out.loc[gidx, "fdr_unadj"] = bh_fdr(out.loc[gidx, "p_unadj"].values)
        out.loc[gidx, "fdr_adj"] = bh_fdr(out.loc[gidx, "p_adj"].values)
    out.to_csv(TBL / "oinfo_age_adjustment.tsv", sep="\t", index=False)

    main_df = out[(out.panel == "A") & (out.model == "A")]
    print(main_df[["region", "beta_age_unadj", "p_unadj", "beta_age_adj", "p_adj",
                   "beta_immune", "p_immune", "attenuation", "fdr_adj"]]
          .sort_values("p_unadj").to_string(index=False))

    # ---------------- Figure 2 ----------------
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    d = main_df.sort_values("beta_age_unadj")

    # 2A slopes unadj vs adj
    ax = axes[0, 0]
    x = np.arange(len(d))
    ax.errorbar(x, d.beta_age_unadj, yerr=[d.beta_age_unadj - d.ci_lo_unadj,
                                           d.ci_hi_unadj - d.beta_age_unadj],
                fmt="o", ms=3, label="unadjusted", color="steelblue")
    ax.errorbar(x, d.beta_age_adj, yerr=[d.beta_age_adj - d.boot_ci_lo,
                                         d.boot_ci_hi - d.beta_age_adj],
                fmt="s", ms=3, label="adjusted", color="crimson")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels([r.replace("_bins_AND_", "\n").replace("_AND_", "\n")
                        for r in d.region], rotation=90, fontsize=4)
    ax.set_ylabel("beta_age (O-info / year)")
    ax.set_title("2A: age slope before/after immune adjustment", fontsize=9)
    ax.legend(fontsize=7)

    # 2B attenuation
    ax = axes[0, 1]
    att = d.attenuation.clip(-2, 2)
    colors = ["crimson" if a < 0 else "steelblue" for a in att]
    ax.bar(x, att, color=colors)
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels([r.replace("_bins_AND_", "\n").replace("_AND_", "\n")
                        for r in d.region], rotation=90, fontsize=4)
    ax.set_ylabel("attenuation = 1 - beta_adj/beta_unadj")
    ax.set_title("2B: attenuation of age effect", fontsize=9)

    # 2C focus region: hoi vs age, colored by immune fraction
    ax = axes[0, 2]
    imm = (prop[(prop.panel == "A") & (prop.model == "A")
                & prop.cell_type.isin(IMM_BY_MODEL["A"])]
           .groupby("sample_id")["fraction"].sum().reset_index()
           .rename(columns={"fraction": "immune_fraction"}))
    fd = oinfo[oinfo.region == FOCUS].merge(meta, on="sample_id").merge(imm, on="sample_id")
    sc = ax.scatter(fd.age, fd.hoi, c=fd.immune_fraction, cmap="viridis",
                    s=70, edgecolor="k", zorder=3)
    plt.colorbar(sc, ax=ax, label="immune fraction")
    f0 = sm.OLS(fd.hoi, sm.add_constant(fd.age.astype(float))).fit()
    xs = np.linspace(fd.age.min() - 2, fd.age.max() + 2, 30)
    ax.plot(xs, f0.params[0] + f0.params[1] * xs, color="crimson", lw=1)
    for _, r in fd.iterrows():
        ax.annotate(r["sample_id"][-4:], (r.age, r.hoi), fontsize=6, alpha=0.7)
    ax.set_xlabel("age (y)")
    ax.set_ylabel("O-information")
    ax.set_title(f"2C: {FOCUS.replace('_bins_AND_','/')}\nr={stats.pearsonr(fd.age, fd.hoi)[0]:.2f}", fontsize=8)

    # 2D hoi vs immune fraction (focus)
    ax = axes[1, 0]
    ax.scatter(fd.immune_fraction, fd.hoi, s=70, color="seagreen", edgecolor="k")
    r_i, p_i = stats.pearsonr(fd.immune_fraction, fd.hoi)
    ax.set_xlabel("estimated immune fraction")
    ax.set_ylabel("O-information")
    ax.set_title(f"2D: O-info vs immune fraction\nr={r_i:.2f}, p={p_i:.3g}", fontsize=9)

    # 2E CI comparison for key regions (focus family)
    ax = axes[1, 1]
    key = d[d.region.str.contains("proximal_ctcf|proximal_not_ctcf")]
    ypos = np.arange(len(key))
    ax.errorbar(key.beta_age_unadj, ypos + 0.15,
                xerr=[key.beta_age_unadj - key.ci_lo_unadj, key.ci_hi_unadj - key.beta_age_unadj],
                fmt="o", ms=4, color="steelblue", label="unadj")
    ax.errorbar(key.beta_age_adj, ypos - 0.15,
                xerr=[key.beta_age_adj - key.boot_ci_lo, key.boot_ci_hi - key.beta_age_adj],
                fmt="s", ms=4, color="crimson", label="adj")
    ax.axvline(0, color="k", lw=0.6)
    ax.set_yticks(ypos)
    ax.set_yticklabels([r.replace("_bins_AND_", "/") for r in key.region], fontsize=6)
    ax.set_xlabel("beta_age")
    ax.set_title("2E: key enhancer regions CIs", fontsize=9)
    ax.legend(fontsize=7)

    # 2F panel/model stability (beta_age_adj across 2 combos)
    ax = axes[1, 2]
    piv = out.pivot_table(index="region", columns=["panel", "model"],
                          values="beta_age_adj")
    ax.scatter(piv[("A", "A")], piv[("B", "B")], s=25, alpha=0.7)
    lim = np.nanmax(np.abs(piv.values)) * 1.1
    ax.plot([-lim, lim], [-lim, lim], "k--", lw=0.8)
    ax.set_xlabel("beta_adj (Panel A / Model A)")
    ax.set_ylabel("beta_adj (Panel B / Model B)")
    ax.set_title("2F: marker-panel stability", fontsize=9)
    r_s = np.corrcoef(piv[("A", "A")], piv[("B", "B")])[0, 1]
    ax.text(0.03, 0.95, f"r = {r_s:.3f}", transform=ax.transAxes, fontsize=9)

    fig.tight_layout()
    fig.savefig(FIG / "fig2_oinfo_age_adjustment.pdf")
    fig.savefig(FIG / "fig2_oinfo_age_adjustment.png", dpi=200)
    plt.close(fig)
    log.info("saved fig2")
    return 0


if __name__ == "__main__":
    sys.exit(main())
