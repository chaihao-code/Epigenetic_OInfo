#!/usr/bin/env python3
"""16_fig_lollipop.py — 棒棒糖图: (beta_unadj - beta_adj) / beta_adj, 40 区间.

指标: 调整前后年龄斜率的相对差 (除以调整后斜率)。
注意: beta_adj 接近 0 时该指标不稳定 (|ratio| 可达 10+), 图中截断到 ±3,
被截断的点在边缘用三角标记并单独注明。
输出: results/figures/fig4_slope_ratio_lollipop.{pdf,png}
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TBL = ROOT / "results/tables"
FIG = ROOT / "results/figures"

CM = 1 / 2.54
plt.rcParams.update({
    "font.family": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 6,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size": 2,
    "ytick.major.size": 2,
    "legend.frameon": False,
    "pdf.fonttype": 42,
})

XLIM = 3.0
SHORT = {  # 仅标注重点区间
    "proximal_ctcf_bound_enhancer_purified_bins_AND_CpG_island":
        "prox. enh. CTCF+ / CpG island",
    "proximal_ctcf_bound_enhancer_purified_bins_AND_SINE":
        "prox. enh. CTCF+ / SINE",
    "proximal_ctcf_bound_enhancer_purified_bins_AND_LINE":
        "prox. enh. CTCF+ / LINE",
    "proximal_ctcf_bound_enhancer_purified_bins_AND_LTR":
        "prox. enh. CTCF+ / LTR",
    "proximal_ctcf_bound_enhancer_purified_bins_AND_None":
        "prox. enh. CTCF+ / None",
}


def main():
    m = pd.read_csv(TBL / "oinfo_age_fdr_table.tsv", sep="\t")
    m["ratio"] = (m.beta_age_unadj - m.beta_age_adj) / m.beta_age_adj
    m = m.sort_values("ratio").reset_index(drop=True)
    n = len(m)
    clipped = m.ratio.abs() > XLIM

    fig, ax = plt.subplots(figsize=(5 * CM, 9.5 * CM))
    fig.subplots_adjust(left=0.10, right=0.97, top=0.965, bottom=0.10)

    ypos = np.arange(n)
    is_ctcf = m.region.str.contains("proximal_ctcf_bound").values
    col = np.where(is_ctcf, "#D55E00", "#9C9C9C")

    for i, r in m.iterrows():
        v = np.clip(r.ratio, -XLIM, XLIM)
        ax.plot([0, v], [i, i], color=col[i], lw=0.6,
                alpha=0.9 if is_ctcf[i] else 0.5, zorder=1)
    ax.scatter(m.ratio.clip(-XLIM, XLIM), ypos, s=7, c=col, zorder=2)
    # 被截断的点用三角标在边缘
    if clipped.any():
        ax.scatter(np.sign(m.loc[clipped, "ratio"]) * XLIM,
                   ypos[clipped.values], s=13, marker=">", c="#D55E00",
                   zorder=3,
                   label=f"clipped (|ratio| > {XLIM:.0f}, n={int(clipped.sum())})")

    # 重点区间名字放在左侧空白区, 细引线连到点 (这些点都紧贴 0)
    for i, r in m.iterrows():
        if r.region in SHORT:
            ax.annotate(SHORT[r.region], (0, i),
                        xytext=(-XLIM * 0.96, i), textcoords="data",
                        fontsize=3.8, va="center", ha="left", color="#B34700",
                        arrowprops=dict(arrowstyle="-", lw=0.3,
                                        color="#BBBBBB", shrinkA=0, shrinkB=2))

    ax.axvline(0, color="#333333", lw=0.5)
    ax.set_xlim(-XLIM, XLIM)
    ax.set_ylim(-1, n)
    ax.set_yticks([])
    ax.set_xlabel(r"$(\beta_{unadj}-\beta_{adj})/\beta_{adj}$", fontsize=5.5)
    ax.tick_params(labelsize=4.5)
    if clipped.any():
        ax.legend(fontsize=4.2, loc="lower right", handletextpad=0.2)
    fig.savefig(FIG / "fig4_slope_ratio_lollipop.pdf")
    fig.savefig(FIG / "fig4_slope_ratio_lollipop.png", dpi=400)
    plt.close(fig)
    print("saved fig4_slope_ratio_lollipop | clipped:",
          int(clipped.sum()), "| median ratio: %.3f" % m.ratio.median())
    return 0


if __name__ == "__main__":
    sys.exit(main())
