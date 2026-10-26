#!/usr/bin/env python3
"""06_immune_age_association.py — Step 7: 免疫比例与年龄的统计关系.

每 panel/model 组合: total_immune_fraction ~ age
  - Pearson / Spearman 相关
  - OLS 斜率 + 95% CI + p
  - 敏感性: leave-one-out (每样本剔除后重估)
  - 协变量检查: sex (5F/1M, 不纳入主模型, 仅报告分布)
n=6, 强调 effect size 与 CI, 避免过度解释显著性。
输出: results/tables/immune_age_regression.tsv
     results/figures/immune_age_regression.pdf
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
from scipy import stats
import yaml

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(Path(__file__).resolve().parents[1] / "logs/07_immune_age.log"),
        logging.StreamHandler(),
    ],
)
log = logging.getLogger("immune_age")

ROOT = Path(__file__).resolve().parents[1]
CFG = _release_config(yaml.safe_load(open(ROOT / "config.yaml")))
TBL = ROOT / "results/tables"
FIG = ROOT / "results/figures"
# 免疫组分按 model 粒度定义 (Model A 细粒度 9 类, Model B 粗粒度 4 类)
IMM_BY_MODEL = {
    "A": ["T_CD4", "T_CD8", "B_cell", "Monocyte", "NK", "Granulocyte"],
    "B": ["Lymphoid", "Myeloid"],
}


def main():
    prop = pd.read_csv(TBL / "thyroid_cell_proportions.tsv", sep="\t")
    meta = pd.read_csv(ROOT / "data/sample_metadata/analysis_samples.tsv", sep="\t")
    rows = []
    fig, axes = plt.subplots(1, 4, figsize=(13, 3.2))
    for ax_i, (panel, model) in enumerate([("A", "A"), ("A", "B"), ("B", "A"), ("B", "B")]):
        d = prop[(prop.panel == panel) & (prop.model == model)]
        imm = (d[d.cell_type.isin(IMM_BY_MODEL[model])]
               .groupby("sample_id")["fraction"].sum().reset_index()
               .rename(columns={"fraction": "immune_fraction"}))
        df = imm.merge(meta, on="sample_id")
        x, y = df.age.values, df.immune_fraction.values
        n = len(df)
        pear_r, pear_p = stats.pearsonr(x, y)
        spr_r, spr_p = stats.spearmanr(x, y)
        X = sm.add_constant(x)
        fit = sm.OLS(y, X).fit()
        slope, intercept = fit.params[1], fit.params[0]
        ci = np.asarray(fit.conf_int())[1]
        slope_lo, slope_hi = ci[0], ci[1]
        # leave-one-out
        loo = []
        for s in df.sample_id:
            m = df.sample_id != s
            f = sm.OLS(y[m], sm.add_constant(x[m])).fit()
            loo.append(dict(left_out=s, slope=f.params[1], p=f.pvalues[1]))
        loo_df = pd.DataFrame(loo)
        rows.append(dict(panel=panel, model=model, n=n,
                         pearson_r=round(pear_r, 4), pearson_p=round(pear_p, 5),
                         spearman_r=round(spr_r, 4), spearman_p=round(spr_p, 5),
                         slope_per_year=round(slope, 6),
                         slope_ci_lo=round(slope_lo, 6), slope_ci_hi=round(slope_hi, 6),
                         ols_p=round(fit.pvalues[1], 5),
                         adj_r2=round(fit.rsquared_adj, 4),
                         loo_slope_min=round(loo_df.slope.min(), 6),
                         loo_slope_max=round(loo_df.slope.max(), 6),
                         loo_p_max=round(loo_df.p.max(), 5)))
        ax = axes[ax_i]
        ax.scatter(x, y, s=45, color="steelblue", edgecolor="k", zorder=3)
        xs = np.linspace(x.min() - 2, x.max() + 2, 50)
        ax.plot(xs, intercept + slope * xs, color="crimson", lw=1)
        ax.set_xlabel("age (y)")
        ax.set_ylabel("total immune fraction")
        ax.set_title(f"Panel {panel} / Model {model}\nr={pear_r:.2f}, p={pear_p:.3g}",
                     fontsize=9)
        for _, r in df.iterrows():
            ax.annotate(r["sample_id"][-4:], (r.age, r.immune_fraction),
                        fontsize=6, alpha=0.7)
    out = pd.DataFrame(rows)
    out.to_csv(TBL / "immune_age_regression.tsv", sep="\t", index=False)
    print(out.to_string(index=False))
    fig.tight_layout()
    fig.savefig(FIG / "immune_age_regression.pdf")
    fig.savefig(FIG / "immune_age_regression.png", dpi=200)
    plt.close(fig)
    log.info("done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
