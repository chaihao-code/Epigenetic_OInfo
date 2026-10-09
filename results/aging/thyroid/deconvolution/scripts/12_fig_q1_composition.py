#!/usr/bin/env python3
"""12_fig_q1_composition.py — 图: 细胞组成 (a) + panel/model 稳健性散点 (b,c).

规格: 总宽 18 cm; 子图宽比 a:b:c = 1.2:1:1; (b)(c) 为方形 (高=宽);
Arial 字体 (缺失时回退 Liberation Sans, 与 Arial 度量兼容), 5 pt.

(a) 主分析 Panel A / Model A 的 9 类细胞组成堆积柱状图 (按年龄排序, 各柱和为 1);
    配色: 高饱和突出甲状腺上皮 / 内皮 / 最大免疫组分 (CD4 T), 其余低饱和.
(b) 总免疫比例: 主分析 (A/A) vs 三种敏感性组合.
(c) 甲状腺上皮比例: 同上.
输出: results/figures/figQ1_composition.{pdf,png}
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
TBL = ROOT / "results/tables"
FIG = ROOT / "results/figures"

# ---------------- 规格 ----------------
CM = 1 / 2.54
TOTAL_W_CM = 18.0
W_B_CM = TOTAL_W_CM / 3.2          # (b)/(c) 宽 5.625 cm -> 方形高同
H_CM = W_B_CM
FIG_W = TOTAL_W_CM * CM
FIG_H = (H_CM + 1.15) * CM        # 底部留给 x 轴标签 + 图例

plt.rcParams.update({
    "font.family": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 5,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size": 2,
    "ytick.major.size": 2,
    "legend.frameon": False,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})

CT_ORDER = ["Thyroid_epithelial", "Endothelial", "T_CD4", "T_CD8",
            "B_cell", "Monocyte", "NK", "Granulocyte", "Fibroblast"]
CT_LABEL = {"Thyroid_epithelial": "Thyroid epithelium", "Endothelial": "Endothelial",
            "T_CD4": "CD4 T", "T_CD8": "CD8 T", "B_cell": "B cell",
            "Monocyte": "Monocyte", "NK": "NK", "Granulocyte": "Granulocyte",
            "Fibroblast": "Fibroblast"}
# 高饱和突出三者: 甲状腺上皮 / 内皮 / 最大免疫组分 (CD4 T); 其余低饱和
CT_COLOR = {"Thyroid_epithelial": "#009E73", "Endothelial": "#0072B2",
            "T_CD4": "#D55E00", "T_CD8": "#E9A03B", "B_cell": "#C08AAD",
            "Monocyte": "#9C9C9C", "NK": "#C4B49A", "Granulocyte": "#BFCFD6",
            "Fibroblast": "#A9C6B8"}

# 敏感性组合 (新命名: 主分析 = Panel A / Model A)
SENS = [("B", "A", "o", "Panel B / Model A (UXM markers)"),
        ("A", "B", "s", "Panel A / Model B (coarse)"),
        ("B", "B", "^", "Panel B / Model B (UXM + coarse)")]
IMM_FINE = ["T_CD4", "T_CD8", "B_cell", "Monocyte", "NK", "Granulocyte"]
IMM_COARSE = ["Lymphoid", "Myeloid"]


def main():
    cfg = yaml.safe_load(open(ROOT / "config.yaml"))
    meta = pd.read_csv(ROOT / "data/sample_metadata/analysis_samples.tsv", sep="\t")
    prop = pd.read_csv(TBL / "thyroid_cell_proportions.tsv", sep="\t")
    order = meta.sort_values("age").sample_id.tolist()
    ages = dict(zip(meta.sample_id, meta.age))

    fig = plt.figure(figsize=(FIG_W, FIG_H))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.2, 1, 1], wspace=0.42,
                          left=0.065, right=0.995, top=0.90, bottom=0.30)
    axa, axb, axc = [fig.add_subplot(gs[i]) for i in range(3)]

    # ---------- (a) stacked composition, Panel A / Model A ----------
    d = prop[(prop.panel == "A") & (prop.model == "A")]
    pv = d.pivot_table(index="cell_type", columns="sample_id", values="fraction")
    bottoms = np.zeros(len(order))
    for ct in CT_ORDER:
        v = pv.loc[ct, order].values * 100
        axa.bar(range(len(order)), v, bottom=bottoms, color=CT_COLOR[ct],
                edgecolor="white", lw=0.25, width=0.72, label=CT_LABEL[ct])
        for i, (b, vi) in enumerate(zip(bottoms, v)):
            if vi >= 5:
                axa.text(i, b + vi / 2, f"{vi:.0f}", ha="center", va="center",
                         fontsize=4.2, color="white",
                         fontweight="bold" if ct == "Thyroid_epithelial" else "normal")
        bottoms += v
    axa.set_ylim(0, 104)
    axa.set_ylabel("Estimated fraction (%)")
    axa.set_xticks(range(len(order)))
    axa.set_xticklabels([f"{s.replace('CEMT','')}\n{ages[s]:.0f} y" for s in order], fontsize=4.5)
    axa.set_xlabel("Sample (sorted by age)")
    axa.legend(ncol=3, fontsize=4.2, loc="upper center",
               bbox_to_anchor=(0.5, -0.30), handlelength=1.1,
               columnspacing=0.9, borderaxespad=0)
    axa.set_title("Panel A / Model A (9 cell types)", fontsize=5, pad=3)

    # ---------- (b) total immune fraction agreement ----------
    main = d[d.cell_type.isin(IMM_FINE)].groupby("sample_id").fraction.sum()
    for pn, md, mk, lab in SENS:
        types = IMM_FINE if md == "A" else IMM_COARSE
        v = (prop[(prop.panel == pn) & (prop.model == md)
                  & prop.cell_type.isin(types)]
             .groupby("sample_id").fraction.sum())
        axb.scatter(main * 100, v.reindex(main.index) * 100, s=9, marker=mk,
                    label=lab, edgecolor="k", lw=0.3)
    lim = 12
    axb.plot([0, lim], [0, lim], ls=(0, (3, 2)), color="gray", lw=0.5)
    axb.set_xlim(-0.4, lim)
    axb.set_ylim(-0.4, lim)
    axb.set_aspect("equal")
    axb.set_xlabel("Total immune fraction —\nmain (Panel A / Model A, %)")
    axb.set_ylabel("Total immune fraction —\nalternative (%)")
    axb.legend(fontsize=4.2, loc="upper left", handletextpad=0.3)
    axb.set_title("Robustness: immune fraction", fontsize=5, pad=3)

    # ---------- (c) thyroid fraction agreement ----------
    thy_main = d[d.cell_type == "Thyroid_epithelial"].set_index("sample_id").fraction
    for pn, md, mk, lab in SENS:
        ct = "Thyroid_epithelial" if md == "A" else "Thyroid"
        v = prop[(prop.panel == pn) & (prop.model == md)
                 & (prop.cell_type == ct)].set_index("sample_id").fraction
        axc.scatter(thy_main * 100, v.reindex(thy_main.index) * 100, s=9,
                    marker=mk, label=lab, edgecolor="k", lw=0.3)
    axc.plot([70, 100], [70, 100], ls=(0, (3, 2)), color="gray", lw=0.5)
    axc.set_xlim(70, 100)
    axc.set_ylim(70, 100)
    axc.set_aspect("equal")
    axc.set_xlabel("Thyroid fraction —\nmain (Panel A / Model A, %)")
    axc.set_ylabel("Thyroid fraction —\nalternative (%)")
    axc.set_title("Robustness: thyroid fraction", fontsize=5, pad=3)

    # panel letters
    for ax, lab in ((axa, "a"), (axb, "b"), (axc, "c")):
        ax.text(-0.16, 1.06, lab, transform=ax.transAxes, fontsize=6.5,
                fontweight="bold", va="bottom", ha="left")

    fig.savefig(FIG / "figQ1_composition.pdf")
    fig.savefig(FIG / "figQ1_composition.png", dpi=400)
    plt.close(fig)
    print("saved figQ1_composition")
    return 0


if __name__ == "__main__":
    sys.exit(main())
