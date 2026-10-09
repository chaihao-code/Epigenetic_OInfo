#!/usr/bin/env python3
"""15_fig_volcano.py — 火山图: 调整后年龄斜率 vs -log10(nominal p).

x: Model 1 (Ω~age+immune) 的 beta_age (主分析 Panel A / Model A, 40 区间)
y: -log10(p_adj, nominal)
参考线: p = 0.05; BH-FDR = 0.1 (本数据对应 nominal p 阈值 = max{p: q_BH<0.1}).
输出: results/figures/fig3_volcano_slope_p.{pdf,png}
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

FOCUS = "proximal_ctcf_bound_enhancer_purified_bins_AND_None"


def main():
    m = pd.read_csv(TBL / "oinfo_age_fdr_table.tsv", sep="\t")
    x = m.beta_age_adj.values
    y = -np.log10(m.p_adj.values)

    p_fdr_thr = m.loc[m.q_bh_adj < 0.1, "p_adj"].max()  # BH-FDR=0.1 对应的 p 阈值
    fig, ax = plt.subplots(figsize=(5 * CM, 5 * CM))
    fig.subplots_adjust(left=0.155, right=0.975, top=0.965, bottom=0.155)

    is_ctcf_prox = m.region.str.contains("proximal_ctcf_bound")
    ax.scatter(x[~is_ctcf_prox], y[~is_ctcf_prox], s=4, c="#9C9C9C",
               alpha=0.75, lw=0, label="other (n=32)")
    ax.scatter(x[is_ctcf_prox], y[is_ctcf_prox], s=6, c="#D55E00",
               alpha=0.95, lw=0, label="prox. enh. CTCF+ (n=5)")

    ax.axhline(-np.log10(0.05), color="#0072B2", lw=0.6, ls=(0, (4, 3)))
    ax.axhline(-np.log10(p_fdr_thr), color="#009E73", lw=0.6, ls=(0, (1, 1.2)))
    ax.set_ylim(0, max(y) * 1.22)
    ax.set_xlim(min(x) * 1.12, max(0.012, max(x) * 1.3))
    ax.text(ax.get_xlim()[0] * 0.30 + ax.get_xlim()[1] * 0.70,
            -np.log10(0.05) + 0.06,
            "p = 0.05", ha="center", va="bottom", fontsize=4.2,
            color="#0072B2")
    ax.text(ax.get_xlim()[0] * 0.40 + ax.get_xlim()[1] * 0.60,
            -np.log10(p_fdr_thr) + 0.06,
            f"FDR = 0.1 (p ≤ {p_fdr_thr:.4f})", ha="center", va="bottom",
            fontsize=4.2, color="#009E73")

    # 标注重点区间
    lab = {
        "proximal_ctcf_bound_enhancer_purified_bins_AND_CpG_island":
            ("prox. enh. CTCF+ / CpG island", (-5, -12), "right"),
        FOCUS: ("prox. enh. CTCF+ / None", (8, -2), "left"),
        "proximal_ctcf_bound_enhancer_purified_bins_AND_LTR":
            ("prox. enh. CTCF+ / LTR", (4, -10), "left"),
        "proximal_ctcf_bound_enhancer_purified_bins_AND_SINE":
            ("prox. enh. CTCF+ / SINE", (5, 1), "left"),
        "proximal_ctcf_bound_enhancer_purified_bins_AND_LINE":
            ("prox. enh. CTCF+ / LINE", (5, -6), "left"),
    }
    for region, (txt, off, ha) in lab.items():
        r = m[m.region == region].iloc[0]
        xi, yi = r.beta_age_adj, -np.log10(r.p_adj)
        ax.annotate(txt, (xi, yi), textcoords="offset points", xytext=off,
                    fontsize=4.0, color="#333333", ha=ha,
                    arrowprops=dict(arrowstyle="-", lw=0.3, color="#888888",
                                    shrinkA=0, shrinkB=1))
    ax.axvline(0, color="#BBBBBB", lw=0.5)

    ax.tick_params(labelsize=4.5)
    ax.set_xlabel(r"adjusted $\beta_{age}$ (Ω / year)", fontsize=5)
    ax.set_ylabel(r"$-\log_{10}$(p), adjusted", fontsize=5)
    ax.legend(fontsize=4.2, loc="upper left", handletextpad=0.2,
              borderaxespad=0.1)
    fig.savefig(FIG / "fig3_volcano_slope_p.pdf")
    fig.savefig(FIG / "fig3_volcano_slope_p.png", dpi=400)
    plt.close(fig)
    print("saved fig3_volcano_slope_p | FDR=0.1 p threshold:", p_fdr_thr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
